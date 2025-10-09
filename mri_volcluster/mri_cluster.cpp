/**
 * @brief Peforms clustering (connected components) on vol or surf; 
 * works on frames (eg, time) as well. 
 *
 *
 */
/*
 * Original Author: Douglas N. Greve
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


/*
  BEGINHELP

  ENDHELP
*/

/* 
  BEGINUSAGE

  ENDUSAGE
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/utsname.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#include "macros.h"
#include "utils.h"
#include "fio.h"
#include "version.h"
#include "cmdargs.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri2.h"
#include "timer.h"
#include "region.h"
#include "resample.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "randomfields.h"
#include "volcluster.h"

#ifdef _OPENMP
#include "romp_support.h"
#endif

MRI *MRIcorify(MRI *seg, double corethresh){
  int segidlist[10000], firstframe[10000], lastframe[10000];

  // For each seg, get the first and last frame it appears in
  int segnomax = -1;
  for(int f=0; f < seg->nframes; f++){
    for(int c=0; c < seg->width; c++){
      for(int r=0; r < seg->height; r++){
	for(int s=0; s < seg->depth; s++){
	  int segid = MRIgetVoxVal(seg,c,r,s,f);
	  if(segid > segnomax) segnomax = segid;
	  if(segidlist[segid] == 0){
	    firstframe[segid] = f;
	    lastframe[segid] = f;
	  }
	  segidlist[segid]++;
	  if(firstframe[segid] > f) firstframe[segid] = f;
	  if(lastframe[segid]  < f) lastframe[segid] = f;
	}
      }
    }
  }

  MRI *core = MRIclone(seg,NULL);
  MRIcopyHeader(seg,core);
  MRIcopyPulseParameters(seg,core);
  if(seg->ct) core->ct = CTABdeepCopy(seg->ct);
  if(seg->nframes == 1) return(core);
  printf("Starting corification %d\n",segnomax);
  for(int segid = 1; segid <= segnomax; segid++){
    if(segidlist[segid] == 0) continue;
    int nf = lastframe[segid] - firstframe[segid] + 1;
    printf("segid %4d %5d  %2d %2d %2d ",segid,segidlist[segid],firstframe[segid],lastframe[segid],nf);
    int nhits = 0;
    for(int c=0; c < seg->width; c++){
      for(int r=0; r < seg->height; r++){
	for(int s=0; s < seg->depth; s++){
	  int nmiss=0;
	  for(int f=firstframe[segid]; f <= lastframe[segid]; f++){
	    int val = MRIgetVoxVal(seg,c,r,s,f);
	    if(val == segid) continue;
	    nmiss ++;
	  }// f
	  // Check if the number of misses exceeds threshold
	  if(nmiss > corethresh*nf) continue;
	  nhits++;
	  // If not enough is missing, then fill in all frames (core)
	  for(int f=firstframe[segid]; f <= lastframe[segid]; f++){
	    int val = MRIgetVoxVal(seg,c,r,s,f);
	    if(val == segid) MRIsetVoxVal(core,c,r,s,f,val);
	  }
	}// s
      }// r
    }// c
    printf(" nhits=%d\n",nhits);
  }//seg
  printf("Done corification\n");
  return(core);
}
#if 0
class SpatTempCluster {
  // test masking, edge and corner
  // clusters - centroid, pointsets, size
public:
  int topo=1; //1 = volume, 2 = surface
  int nbrtype = 0; //1=face, 2=+edge, 3=+corner
  MRI *binmask=NULL;
  MRI *cnomap =NULL;
  MRIS *surf=NULL;
  int debug = 0;
  std::vector<std::vector<int>> voxlist; //c,r,s,f
  std::vector<std::vector<int>> GetNearestNeighbors(std::vector<int> vox){
    std::vector<std::vector<int>> nbrs;
    if(topo==1){ // volume
      for(int dc=-1; dc <= +1; dc++){
	for(int dr=-1; dr <= +1; dr++){
	  for(int ds=-1; ds <= +1; ds++){
	    if(abs(dc)+abs(dr)+abs(ds) > nbrtype) continue;
	    if(vox[0]+dc < 0 || vox[0]+dc >= binmask->width) continue;
	    if(vox[1]+dr < 0 || vox[1]+dr >= binmask->height) continue;
	    if(vox[2]+ds < 0 || vox[2]+ds >= binmask->depth) continue;
	    for(int df=-1; df <= +1; df++){
	      if(abs(dc)+abs(dr)+abs(ds)+abs(df)==0) continue; // not self
	      // this line forces a face topology across frames
	      if(abs(df) > 0 && (abs(dc)+abs(dr)+abs(ds) > 0) ) continue;
	      if(vox[3]+df < 0 || vox[3]+df >= binmask->nframes) continue;
	      std::vector<int> dcrsf = {vox[0]+dc,vox[1]+dr,vox[2]+ds,vox[3]+df};
	      nbrs.push_back(dcrsf);
	    }
	  }
	}
      }
    } else { // surface 
      int vno = vox[0];
      // Doing frames separately forces face topology across frame
      for(int df=-1; df <= +1; df++){
	if(vox[3]+df < 0 || vox[3]+df >= binmask->nframes) continue;
	std::vector<int> dcrsf = {vox[0],0,0,vox[3]+df};
	nbrs.push_back(dcrsf);
      }
      VERTEX_TOPOLOGY *vtop = &surf->vertices_topology[vno];
      int vnum = vtop->vtotal;
      for(int i = 0; i < vnum; i++) {
	int vnbno = vtop->v[i];
	std::vector<int> dcrsf = {vnbno,0,0,vox[3]};
	nbrs.push_back(dcrsf);
      }
    }
    return(nbrs);
  }
  class Cluster {
  public: 
    int cno = 0;
    int nmembers = 0;
    std::vector<std::vector<int>> crst;
  };
  std::vector<Cluster> ClusterList;
  int GrowOne(std::vector<int> vox, int cno)
  {
    int m;
    m = MRIgetVoxVal(this->binmask,vox[0],vox[1],vox[2],vox[3]);
    if(!m) return(0); // not a set voxel
    m = MRIgetVoxVal(this->cnomap,vox[0],vox[1],vox[2],vox[3]);
    if(m) return(0); // already in a cluster
    // To get here, the vox must be active in the binmask and not in a cluster already
    MRIsetVoxVal(this->cnomap,vox[0],vox[1],vox[2],vox[3], cno+1);
    this->ClusterList[cno].crst.push_back(vox);
    int nhits = 1;
    std::vector<std::vector<int>> nbrs = this->GetNearestNeighbors(vox);
    for(int n=0; n < nbrs.size(); n++){
      nhits += this->GrowOne(nbrs[n],cno);
    }
    this->ClusterList[cno].nmembers += nhits;
    return(nhits);
  }
  int Clusterize(void){
    if(this->cnomap) MRIfree(&this->cnomap);
    this->cnomap = MRIclone(this->binmask,NULL);
    this->ClusterList.clear();
    int pointno = -1;
    int cno = -1;
    while(1){
      pointno ++;
      if(pointno >= this->voxlist.size()) break;
      int done = 0;
      std::vector<int> vox;
      while(1){ // find the next point that has not been marked
	vox = voxlist[pointno];
	int val = MRIgetVoxVal(this->cnomap,vox[0],vox[1],vox[2],vox[3]);
	if(val) { // point already found
	  pointno++;
	  if(pointno >= this->voxlist.size()) {
	    done = 1;
	    break;
	  }
	}
	else break; // found an unmarked point
      }
      if(done) break;
      cno++;
      Cluster cl;
      cl.cno = cno;
      if(debug) printf("Adding cno=%d pointno=%d  %d %d %d %d =====\n",cno,pointno,vox[0],vox[1],vox[2],vox[3]);
      this->ClusterList.push_back(cl);
      this->GrowOne(vox,cno);
    }
    printf("Found %d clusters\n",(int)this->ClusterList.size());
    printf("Adding ctab\n");
    this->cnomap->ct = CTABalloc(this->ClusterList.size()+1);
    CTABunique(this->cnomap->ct, 100); //100 = number of tries

    return(this->ClusterList.size());
  }
  int PrintClusterSum(FILE *fp){
    for(int n=0; n < this->ClusterList.size(); n++){
      Cluster cl = ClusterList[n];
      fprintf(fp,"%2d %5d   %3d %3d %3d  %3d\n",n+1,(int)cl.crst.size(),
	      cl.crst[0][0],cl.crst[0][1],cl.crst[0][2],cl.crst[0][3]);
    }
    fflush(fp);
    return(0);
  }
  int GetBinMask(MRI *ov, double thmin, double thmax, int sign, MRI *mask)
  {
    if(this->binmask) MRIfree(&this->binmask);
    this->binmask = MRIallocSequence(ov->width,ov->height,ov->depth,MRI_INT,ov->nframes);
    MRIcopyHeader(ov,this->binmask);
    this->voxlist.clear();
    //not thread safe
    for(int c=0; c < ov->width; c++){
      for(int r=0; r < ov->height; r++){
	for(int s=0; s < ov->depth; s++){
	  if(mask && MRIgetVoxVal(mask,c,r,s,0)<0.5) continue;
	  for(int f=0; f < ov->nframes; f++){
	    double val = MRIgetVoxVal(ov,c,r,s,f);
	    if(sign ==  0) val = fabs(val);
	    if(sign == -1) val = -val;
	    if(val < thmin || val > thmax) continue;
	    std::vector<int> vox = {c,r,s,f};
	    voxlist.push_back(vox);
	    MRIsetVoxVal(binmask,c,r,s,f,1);
	  }
	}
      }
    }
    printf("thmin=%g thmax=%g sign=%d nhits=%d\n",thmin,thmax,sign,(int)voxlist.size());
    return(voxlist.size());
  }

};
#endif


static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

const char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;
char *ovfile = NULL;
double thmin = -std::numeric_limits<double>::infinity();
double thmax = +std::numeric_limits<double>::infinity();
int thsign = 0; // abs
char *maskfile = NULL;
char *surffile = NULL;
char *outdir = NULL;
int threads = 1;
char *SUBJECTS_DIR=NULL;
SpatTempCluster stc;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int nargs, err=0;
  char fname[1000];
  Timer timer, mytimer;

  nargs = handleVersionOption(argc, argv, "mri_gtmpvc");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  if (argc == 0) usage_exit();
  parse_commandline(argc, argv);
  check_options();
  if (checkoptsonly) return(0);
  dump_options(stdout);
  setenv("FS_COPY_HEADER_CTAB","1",1);

#ifdef _OPENMP
  printf("%d avail.processors, using %d\n",omp_get_num_procs(),omp_get_max_threads());
#endif

  MRI *ov = MRIread(ovfile);
  if(!ov) exit(1);

  stc.topo = 1;
  if(surffile) {
    stc.surf = MRISread(surffile);
    if(!stc.surf) exit(1);
    stc.topo = 2;
  }
  MRI *mask = NULL;
  if(maskfile){
    mask = MRIread(maskfile);
    if(!mask) exit(1);
  }

  printf("Creating output directory %s\n",outdir);
  err = mkdir(outdir,0777);
  if(err != 0 && errno != EEXIST) {
    printf("ERROR: creating directory %s\n",outdir);
    perror(NULL);
    return(1);
  }

  stc.GetBinMask(ov,thmin,thmax,thsign,mask);
  stc.Clusterize();

  sprintf(fname,"%s/ocn.mgz",outdir);
  MRIwrite(stc.cnomap,fname);
  sprintf(fname,"%s/ctab.ocn",outdir);
  CTABwriteFileASCII(stc.cnomap->ct,fname);

  sprintf(fname,"%s/binmask.nii.gz",outdir);
  MRIwrite(stc.binmask,fname);
  sprintf(fname,"%s/clusters.dat",outdir);
  FILE *fp = fopen(fname,"w");
  stc.PrintClusterSum(fp);
  fclose(fp);

  printf("#VMPC# mri_cluster VmPeak  %d\n",GetVmPeak());
  printf("mris_cluster-runtime %5.2f min\n",timer.minutes());
  printf("mris_cluster done\n");
  return(0);
  exit(0);

} // end of main

/*--------------------------------------------------------------------*/
/*---------------------------------------------------------------*/
/*---------------------------------------------------------------*/
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;

  if (argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while (nargc > 0) {

    option = pargv[0];
    if(debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if(!strcasecmp(option, "--help"))  print_help() ;
    else if(!strcasecmp(option, "--version")) print_version() ;
    else if(!strcasecmp(option, "--debug"))   debug = 1; 
    else if(!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if(!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;
    else if(!strcasecmp(option, "--o")) {
      if(nargc < 1) CMDargNErr(option,1);
      outdir = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--i")) {
      if(nargc < 1) CMDargNErr(option,1);
      ovfile = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--thmin")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&thmin);
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--thmax")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&thmax);
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--abs")) thsign =  0;
    else if(!strcasecmp(option, "--pos")) thsign = +1;
    else if(!strcasecmp(option, "--neg")) thsign = -1;
    else if(!strcasecmp(option, "--surf")) {
      if(nargc < 1) CMDargNErr(option,1);
      surffile = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--mask")) {
      if(nargc < 1) CMDargNErr(option,1);
      maskfile = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--face"))   stc.nbrtype = 1;
    else if(!strcasecmp(option, "--edge"))   stc.nbrtype = 2;
    else if(!strcasecmp(option, "--corner")) stc.nbrtype = 3;
    else if(!strcasecmp(option, "--threads")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&threads);
      #ifdef _OPENMP
      omp_set_num_threads(threads);
      #endif
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--max-threads")){
      threads = 1;
      #ifdef _OPENMP
      threads = omp_get_max_threads();
      omp_set_num_threads(threads);
      #endif
    } 
    else if(!strcasecmp(option, "--max-threads-1") || !strcasecmp(option, "--max-threads-minus-1")){
      threads = 1;
      #ifdef _OPENMP
      threads = omp_get_max_threads()-1;
      if(threads < 0) threads = 1;
      omp_set_num_threads(threads);
      #endif
    } 
    else if(!strcasecmp(option, "--corify")) {
      // standalone invol thresh outvol
      if(nargc < 3) CMDargNErr(option,2);
      MRI *invol = MRIread(pargv[0]);
      if(!invol) exit(1);
      double thresh;
      sscanf(pargv[1],"%lf",&thresh);
      printf("corify thresh %g\n",thresh);
      if(thresh < 0 || thresh > 1) {
	printf("ERROR: threshold must be between 0 and 1\n");
	exit(1);
      }
      MRI *core = MRIcorify(invol,thresh);
      if(!core) exit(1);
      int err = MRIwrite(core,pargv[2]);
      exit(err);
      nargsused = 2;
    }
    else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if(CMDsingleDash(option))
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/*---------------------------------------------------------------*/
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/*---------------------------------------------------------------*/
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --o outdir\n");
  printf("   --i overlay\n");
  printf("   --mask mask\n");
  printf("   --surf surffile : when overlay is a surface\n");
  printf("   --thmin thmin (default is -infinity)\n");
  printf("   --thmax thmax (default is +infinity)\n");
  printf("   --abs, --pos, --neg\n");
  printf("   --face, --edge, --cornder : neighbor definition (vol only)\n");
  printf("   --corify tsegvol thresh corevol (thresh 0-1, 0 more strict)\n");
  #ifdef _OPENMP
  printf("   --threads N : use N threads (with Open MP)\n");
  printf("   --max-threads : use the maximum allowable number of threads for this computer\n");
  printf("   --max-threads-minus-1 : use one less than the maximum allowable number of threads for this computer\n");
  #endif
  //printf("   --sd SUBJECTS_DIR\n");
  //printf("   --gdiag diagno : set diagnostic level\n");
  //printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
  printf("\n");
}
/*---------------------------------------------------------------*/
static void print_help(void) {
  print_usage() ;
  exit(1) ;
}
/*---------------------------------------------------------------*/
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/*---------------------------------------------------------------*/
static void check_options(void) 
{
  if(outdir == NULL){
    printf("ERROR: must spec outdir\n");
    exit(1);
  }
  if(ovfile == NULL){
    printf("ERROR: must spec input overlay\n");
    exit(1);
  }
  if(stc.nbrtype && surffile){
    printf("ERROR: cannot use --face, --edge, or --corner with surfaces\n");
    exit(1);
  } else stc.nbrtype = 1; // use face by default
  
  stc.debug = debug;

  return;
}
/*---------------------------------------------------------------*/
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"setenv SUBJECTS_DIR %s\n",SUBJECTS_DIR);
  fprintf(fp,"cd %s\n",cwd);
  fprintf(fp,"%s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  fprintf(fp,"overlay   %s\n",ovfile);
  if(maskfile) fprintf(fp,"mask  %s\n",maskfile);
  if(surffile) fprintf(fp,"mask  %s\n",surffile);
  fprintf(fp,"thmin  %8.4lf\n",thmin);
  fprintf(fp,"thmax  %8.4lf\n",thmax);
  fprintf(fp,"thsign  %d\n",thsign);
  fprintf(fp,"nbrtype %d\n",stc.nbrtype);
  return;
}
