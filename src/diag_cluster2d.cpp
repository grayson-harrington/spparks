/* -----------------------q-----------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "output.h"
#include "memory.h"
#include "app.h"
#include "error.h"
#include "timer.h"
#include "diag_cluster2d.h"
#include "app_lattice2d.h"
#include "comm_lattice2d.h"
using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

DiagCluster2d::DiagCluster2d(SPK *spk, int narg, char **arg) : Diag(spk,narg,arg)
{
  fp = NULL;
  fpdump = NULL;
  clustlist = NULL;
  ncluster = 0;
  idump = 0;
  dump_style = STANDARD;

  if (narg < 3 || narg > 4) error->all("Illegal diag_style cluster2d command");

  int iarg = 2;
  if (me == 0) {
    fp = fopen(arg[iarg],"w");
    if (!fp) error->one("Cannot open diag_style cluster2d output file");
  }
  iarg++;

  if (narg == 4) {
    dump_style = STANDARD;
    idump = 1;
    if (me == 0) {
      fpdump = fopen(arg[iarg],"w");
      if (!fpdump) error->one("Cannot open diag_style cluster2d dump file");
    }
  }
}

/* ---------------------------------------------------------------------- */

DiagCluster2d::~DiagCluster2d()
{
  memory->destroy_2d_T_array(cluster_ids,nxlo,nylo);
  free_clustlist();
  if (me == 0 ) {
    if (fp) fclose(fp);
    if (fpdump) fclose(fpdump);
  }
}

/* ---------------------------------------------------------------------- */

void DiagCluster2d::init(double time)
{
  applattice2d = (AppLattice2d *) app;
  nx_global = applattice2d->nx_global;
  ny_global = applattice2d->ny_global;
  nx_procs = applattice2d->nx_procs;
  ny_procs = applattice2d->ny_procs;
  delghost = applattice2d->delghost;
  nx_local = applattice2d->nx_local;
  ny_local = applattice2d->ny_local;
  nx_offset = applattice2d->nx_offset;
  ny_offset = applattice2d->ny_offset;
  nxlo = applattice2d->nxlo;
  nylo = applattice2d->nylo;
  nxhi = applattice2d->nxhi;
  nyhi = applattice2d->nyhi;

  memory->create_2d_T_array(cluster_ids,nxlo,nxhi,nylo,nyhi,
			    "diagcluster2d:cluster");

  write_header();
  analyze_clusters(time);
  diag_time = time + diag_delta;
}


/* ---------------------------------------------------------------------- */

void DiagCluster2d::compute(double time, int done)
{
  if ((diag_delta > 0.0 && time >= diag_time) || done) {
    analyze_clusters(time);
    diag_time += diag_delta;
  }
}

/* ---------------------------------------------------------------------- */

void DiagCluster2d::analyze_clusters(double time)
{
  if (me == 0) {
    fprintf(fp,"\n\n--------------------------------------------------\n");
    fprintf(fp,"Time = %f \n",time);
  }
  generate_clusters();
  if (idump) dump_clusters(time);
  free_clustlist();
}
/* ---------------------------------------------------------------------- */

void DiagCluster2d::write_header()
{
  if (me == 0) {
    fprintf(fp,"Clustering Analysis for 2D Lattice (diag_style cluster2d) \n");
    fprintf(fp,"App = %s \n",style);
    fprintf(fp,"nx_global = %d \n",nx_global);
    fprintf(fp,"ny_global = %d \n",ny_global);
    fprintf(fp,"nprocs = %d \n",nprocs);
  }
}


/* ----------------------------------------------------------------------
   Perform cluster analysis using a definition
   of connectivity provided by the child application
------------------------------------------------------------------------- */

void DiagCluster2d::generate_clusters()
{

  // Psuedocode
  //
  // work on copy of spin array since spin values are changed
  // all id values = 0 initially

  // loop n over all owned sites {
  //   if (id[n] = 0) {
  //     area = 0
  //   } else continue

  //   push(n)
  //   while (stack not empty) {
  //     i = pop()
  //     loop over m neighbor pixels of i:
  //       if (spin[m] == spin[n] and not outside domain) push(m)
  //   }

  // void push(m)
  //   id[m] = id
  //   stack[nstack++] = m

  // int pop()
  //   return stack[nstack--]
  //   area++

  // Set ghost site ids to -1
  // Use four equal sized rectangles arranged in a ring
  for (int j = 1-delghost; j <= ny_local; j++) {
    for (int i = 1-delghost; i <= 0; i++) {
      cluster_ids[i][j] = -1;
    }
  }
  for (int j = ny_local+1; j <= ny_local+delghost; j++) {
    for (int i = 1-delghost; i <= nx_local; i++) {
      cluster_ids[i][j] = -1;
    }
  }
  for (int j = 1; j <= ny_local+delghost; j++) {
    for (int i = nx_local+1; i <= nx_local+delghost; i++) {
      cluster_ids[i][j] = -1;
    }
  }
  for (int j = 1-delghost; j <= 0; j++) {
    for (int i = 1; i <= nx_local+delghost; i++) {
      cluster_ids[i][j] = -1;
    }
  }

  // Set local site ids to zero 
  for (int j = 1; j <= ny_local; j++) {
    for (int i = 1; i <= nx_local; i++) {
      cluster_ids[i][j] = 0;
    }
  }

  int vol,volsum,voltot,nclustertot,ii,jj,id;

  ncluster = 0;
  volsum = 0;

  // loop over all owned sites
  for (int i = 1; i <= nx_local; i++) {
    for (int j = 1; j <= ny_local; j++) {

      // If already visited, skip
      if (cluster_ids[i][j] != 0) {
	continue;
      }

      // Push first site onto stack
      id = ncluster+1;
      vol = 0;
      add_cluster(id,vol,0,NULL);
      cluststack.push(i);
      cluststack.push(j);
      cluster_ids[i][j] = id;

      while (cluststack.size()) {
	// First top then pop
	jj = cluststack.top();
	cluststack.pop();
	ii = cluststack.top();
	cluststack.pop();
	vol++;
	applattice2d->push_connected_neighbors(ii,jj,cluster_ids,ncluster,&cluststack);
      }
      clustlist[ncluster-1].volume = vol;
      volsum+=vol;
    }
  }

  int idoffset;
  MPI_Allreduce(&volsum,&voltot,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&ncluster,&nclustertot,1,MPI_INT,MPI_SUM,world);
  MPI_Scan(&ncluster,&idoffset,1,MPI_INT,MPI_SUM,world);
  idoffset = idoffset-ncluster+1;
  for (int i = 0; i < ncluster; i++) {
    clustlist[i].global_id = i+idoffset;
  }

  // change site ids to global ids
  for (int i = 1; i <= nx_local; i++) {
    for (int j = 1; j <= ny_local; j++) {
      cluster_ids[i][j] = clustlist[cluster_ids[i][j]-1].global_id;
    }
  }

  applattice2d->comm->all(cluster_ids);

  // loop over all owned sites adjacent to boundary
  for (int i = 1; i <= nx_local; i++) {
    for (int j = 1; j <= ny_local; j++) {

      applattice2d->connected_ghosts(i,j,cluster_ids,clustlist,idoffset);

    }
  }

  // pack my clusters into buffer
  int me_size,m,maxbuf;
  int* ibufclust;
  int tmp,nrecv;
  MPI_Status status;
  MPI_Request request;
  int nn;
  
  me_size = 0;
  for (int i = 0; i < ncluster; i++) {
    me_size += 3+clustlist[i].nneigh;
  }
  if (me == 0) me_size = 0;

  MPI_Allreduce(&me_size,&maxbuf,1,MPI_INT,MPI_MAX,world);

  ibufclust = new int[maxbuf];

  if (me != 0) {
    m = 0;
    for (int i = 0; i < ncluster; i++) {
      ibufclust[m++] = clustlist[i].global_id;
      ibufclust[m++] = clustlist[i].volume;
      ibufclust[m++] = clustlist[i].nneigh;
      for (int j = 0; j < clustlist[i].nneigh; j++) {
	ibufclust[m++] = clustlist[i].neighlist[j];
      }
    }
    
    if (me_size != m) {
      error->one("Mismatch in counting for ibufclust");
    }

  }

  // proc 0 pings each proc, receives it's data, adds it to list
  // all other procs wait for ping, send their data to proc 0

  if (me == 0) {
    for (int iproc = 1; iproc < nprocs; iproc++) {
      MPI_Irecv(ibufclust,maxbuf,MPI_INT,iproc,0,world,&request);
      MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
      MPI_Wait(&request,&status);
      MPI_Get_count(&status,MPI_INT,&nrecv);
      
      m = 0;
      while (m < nrecv) {
	id = ibufclust[m++];
	vol = ibufclust[m++];
	nn = ibufclust[m++];
	add_cluster(id,vol,nn,&ibufclust[m]);
	m+=nn;
	volsum+=vol;
      }
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(ibufclust,me_size,MPI_INT,0,0,world);
  }

  delete [] ibufclust;

  // Perform cluster analysis on the clusters

  if (me == 0) {
    int* neighs;
    int jneigh,ncluster_reduced;

    volsum = 0;
    ncluster_reduced = 0;

    // loop over all clusters
    for (int i = 0; i < ncluster; i++) {
      
      // If already visited, skip
      if (clustlist[i].volume == 0) {
	continue;
      }
      
      // Push first cluster onto stack
      id = clustlist[i].global_id;
      vol = 0;
      ncluster_reduced++;
      
      cluststack.push(i);
      vol+=clustlist[i].volume;
      clustlist[i].volume = 0;
      
      while (cluststack.size()) {
	// First top then pop
	ii = cluststack.top();
	cluststack.pop();
	
	neighs = clustlist[ii].neighlist;
	for (int j = 0; j < clustlist[ii].nneigh; j++) {
	  jneigh = neighs[j]-idoffset;
	  if (clustlist[jneigh].volume != 0) {
	    cluststack.push(jneigh);
	    vol+=clustlist[jneigh].volume;
	    clustlist[jneigh].global_id = id;
	    clustlist[jneigh].volume = 0;
	  }
	}
      }
      clustlist[i].volume = vol;
      volsum+=vol;
    }
    
    fprintf(fp,"ncluster = %d \nsize = ",ncluster_reduced);
    for (int i = 0; i < ncluster; i++) {
      if (clustlist[i].volume > 0.0) {
	fprintf(fp," %d",clustlist[i].volume);
      }
    }
    fprintf(fp,"\n");
  }
  

}


/* ---------------------------------------------------------------------- */

void DiagCluster2d::add_cluster(int id, int vol, int nn, int* neighs)
{
  // grow cluster array

  ncluster++;
  clustlist = (Cluster *) memory->srealloc(clustlist,ncluster*sizeof(Cluster),
					 "diagcluster2d:clustlist");
  clustlist[ncluster-1] = Cluster(id,vol,nn,neighs);
}

/* ----------------------------------------------------------------------
   dump a snapshot of cluster identities to file pointer fpdump
------------------------------------------------------------------------- */

void DiagCluster2d::dump_clusters(double time)
{
  int nsend,nrecv,nxtmp,nytmp,nxotmp,nyotmp;
  int size_one = 1;
  int* buftmp;
  int maxbuftmp;
  int id;
  int isite;

  if (me == 0) {
    fprintf(fpdump,"ITEM: TIME\n");
    fprintf(fpdump,"%g\n",time);
    fprintf(fpdump,"ITEM: DIMENSIONS\n");
    fprintf(fpdump,"%d\n%d\n",nx_global,ny_global);
    fprintf(fpdump,"ITEM: ELEMENT CLUSTERID\n");
  }

  // set up communication buffer
  // maxbuftmp must equal the maximum number of spins on one domain 
  // plus some extra stuff
  maxbuftmp = (nx_global/nx_procs+1)*(ny_global/ny_procs+1)+4;
  nsend = nx_local*ny_local+4;
  if (maxbuftmp < nsend) 
    error->one("maxbuftmp size too small in DiagCluster2d::dump_clusters()");
  
  buftmp = (int*) memory->smalloc(maxbuftmp*sizeof(int),"diagcluster2d:dump_clusters:buftmp");

  // proc 0 writes interactive dump header

  // proc 0 writes timestep header

  int m = 0;

  // pack local layout info into buffer

  buftmp[m++] = nx_local;
  buftmp[m++] = ny_local;
  buftmp[m++] = nx_offset;
  buftmp[m++] = ny_offset;

  // pack my lattice values into buffer
  // Violating normal ordering to satisfy output convention

  for (int j = 1; j <= ny_local; j++) {
    for (int i = 1; i <= nx_local; i++) {
      buftmp[m++] = cluster_ids[i][j];
    }
  }

  // proc 0 pings each proc, receives it's data, writes to file
  // all other procs wait for ping, send their data to proc 0

  int tmp;
  MPI_Status status;
  MPI_Request request;

  if (me == 0) {
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
	MPI_Irecv(buftmp,maxbuftmp,MPI_INT,iproc,0,world,&request);
	MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_INT,&nrecv);
      } else nrecv = nsend;
      
      m = 0;
      nxtmp = buftmp[m++];
      nytmp = buftmp[m++];
      nxotmp = buftmp[m++];
      nyotmp = buftmp[m++];

  // print lattice values
  // isite = global grid cell (1:Nglobal)
  // ordered fast in x, slow in y

      for (int j = 1; j <= nytmp; j++) {
	for (int i = 1; i <= nxtmp; i++) {
	  isite = (nyotmp+j-1)*nx_global + (nxotmp+i-1) + 1;
	  id = clustlist[buftmp[m++]-1].global_id;
	  fprintf(fpdump,"%3d %3d \n",isite,id);
	}
      }
    }
    
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(buftmp,nsend,MPI_INT,0,0,world);
  }

  memory->sfree(buftmp);
}


/* ----------------------------------------------------------------------
   dump a snapshot of cluster identities to the screen in 2D layout
   all the cluster identities for each processor domain are printed out
------------------------------------------------------------------------- */

void DiagCluster2d::dump_clusters_detailed()
{
  int nsend,nrecv,nxtmp,nytmp,nztmp,nxhtmp,nyhtmp,nzhtmp,nxotmp,nyotmp,nzotmp;
  int size_one = 1;
  int* buftmp;
  int maxbuftmp;
  int id;

  // set up communication buffer
  // maxbuftmp must equal the maximum number of spins on one domain 
  // plus some extra stuff
  maxbuftmp = ((nx_global-1)/nx_procs+1+2*delghost)*((ny_global-1)/ny_procs+1+2*delghost)+9;
  nsend = (nx_local+2*delghost)*(ny_local+2*delghost)+9;
  if (maxbuftmp < nsend) 
    error->one("maxbuftmp size too small in DiagCluster2d::dump_clusters()");
  
  buftmp = (int*) memory->smalloc(maxbuftmp*sizeof(int),"diagcluster2d:dump_clusters:buftmp");

  // proc 0 writes interactive dump header

  if (me == 0) {
    fprintf(fpdump,"*** Cluster Dump ***\n");
    fprintf(fpdump,"nx_global = %d ny_global = %d\n",nx_global,ny_global);
  }

  int m = 0;

  // pack local layout info into buffer

  buftmp[m++] = nx_local;
  buftmp[m++] = ny_local;
  buftmp[m++] = 0;
  // Need to delete these two
  buftmp[m++] = 0;
  buftmp[m++] = 0;
  buftmp[m++] = 0;
  buftmp[m++] = nx_offset;
  buftmp[m++] = ny_offset;
  buftmp[m++] = 0;

  // pack my lattice values into buffer
  // Need to violate normal ordering in order to simplify output

  for (int j = 1-delghost; j <= ny_local+delghost; j++) {
    for (int i = 1-delghost; i <= nx_local+delghost; i++) {
      buftmp[m++] = cluster_ids[i][j];
    }
  }

  // proc 0 pings each proc, receives it's data, writes to file
  // all other procs wait for ping, send their data to proc 0

  int tmp;
  MPI_Status status;
  MPI_Request request;

  if (me == 0) {
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
	MPI_Irecv(buftmp,maxbuftmp,MPI_INT,iproc,0,world,&request);
	MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_INT,&nrecv);
      } else nrecv = nsend;
      
      m = 0;
      nxtmp = buftmp[m++];
      nytmp = buftmp[m++];
      nztmp = buftmp[m++];
      nxhtmp = buftmp[m++];
      nyhtmp = buftmp[m++];
      nzhtmp = buftmp[m++];
      nxotmp = buftmp[m++];
      nyotmp = buftmp[m++];
      nzotmp = buftmp[m++];
      fprintf(fpdump,"iproc = %d \n",iproc);
      fprintf(fpdump,"nxlocal = %d \nnylocal = %d \n",nxtmp,nytmp);
      fprintf(fpdump,"nx_offset = %d \nny_offset = %d \n",nxotmp,nyotmp);
      m = nrecv;
      for (int j = nytmp+delghost; j >= 1-delghost; j--) {
	m-=nxtmp+2*delghost;
	for (int i = 1-delghost; i <= nxtmp+delghost; i++) {
	  // work out the correct global id of this site
	  id = clustlist[buftmp[m++]-1].global_id;
	  fprintf(fpdump,"%3d",id);
	}
	fprintf(fpdump,"\n");
	m-=nxtmp+2*delghost;
      }
    }
    
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(buftmp,nsend,MPI_INT,0,0,world);
  }

  memory->sfree(buftmp);
}

/* ---------------------------------------------------------------------- */

void DiagCluster2d::free_clustlist()
{
  // Can not call Cluster destructor, because 
  // that would free memory twice.
  // Instead, need to delete neighlist manually.
  for (int i = 0; i < ncluster; i++) {
    free(clustlist[i].neighlist);
  }
  memory->sfree(clustlist);
  clustlist = NULL;
  ncluster = 0;
}
