/**
 * @brief GLM fit for multi-variate data
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

/*
  ToDo:
  1. Spec contrast to perm if desired (rather than all)
  2. Spec whether to use residual or not
  3. Allow arrays to threshold values and signs
  4. Done, sort of. Change/Add SpatTempClust to compute size as mm2 or mm3 tricky
  5. Done. Mem leak, not sure how serious or if there is something that can be done
  6. Run separately after running simple GLM
  7. Done. Handle surf and vol better
  8. Done. multiple threads
  9. Done: ctab - hanges when trying to find ctab
 10. Write out point set
 11. Keep cluster list around so it can be written out with perm p-value
 12. Sort list by size, allow thresholding of list
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
#include "version.h"
#include "cmdargs.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri2.h"
#include "timer.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "volcluster.h"
#include "fsglm.h"
#include "fsgdf.h"
#include "fmriutils.h"
#include "pdf.h"

class MVGLMPERM {
public:
  GLMMAT *glm=NULL;
  MATRIX *ymat=NULL;
  MRI *templt=NULL;
  MRI *mask=NULL;
  int permsign=0, permshuffle=0;
  int UseResid=0;
  std::vector<double> threshlist;
  std::vector<int> signlist;
};

#ifdef _OPENMP
#include "romp_support.h"
#endif

MATRIX *fMRIarrayToMatrix(std::vector<MRI*> mriArray, MRI *mask, int major=1);

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
std::vector<char*> yflist;
char *maskfile=NULL;
char *surffile = NULL;

char *fsgdfile=NULL;
int DoOSGM = 0;
char *outdir=NULL;
std::vector<double> threshlist;
int threads = 1;
char *SUBJECTS_DIR=NULL;
int prunemask = 1;
double prune_thr = FLT_MIN;
int nperm=0;
unsigned long seed = 0;
double thmin = 1.3;
double thmax = +std::numeric_limits<double>::infinity();
int thsign = 0; // abs
int arraymajor=1;
char *ymatfile = NULL;
int cpermno = 0;
int nbrtype = 3;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int nargs, err=0;
  char fname[1000];
  Timer timer, mytimer;
  char logfile[1000];

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

  MRIS *surf = NULL;
  if(surffile){
    printf("Reading surf  %s\n",surffile);
    surf = MRISread(surffile);
    if(!surf) exit(1);
  }

  MRI *mask = NULL;
  if(maskfile){
    printf("Reading mask  %s\n",maskfile);
    mask = MRIread(maskfile);
    if(!mask) exit(1);
  }

  std::vector<MRI *> yarray;
  for(int n=0; n < yflist.size(); n++){
    printf("Reading %d %s\n",n,yflist[n]);
    MRI *y = MRIread(yflist[n]);
    if(y==NULL) exit(1);
    yarray.push_back(y);
  }
  int nvariates;
  if(arraymajor == 1) nvariates = yarray.size(); // number of time points
  else                nvariates = yarray[0]->nframes;
  printf("nvariates = %d\n",nvariates);

  printf("Creating output directory %s\n",outdir);
  err = mkdir(outdir,0777);
  if(err != 0 && errno != EEXIST) {
    printf("ERROR: creating directory %s\n",outdir);
    perror(NULL);
    return(1);
  }

  sprintf(logfile,"%s/mri_mvglmfit.log",outdir);
  FILE *logfp = fopen(logfile,"w");
  dump_options(logfp);

  sprintf(fname,"%s/seed.txt",outdir);
  FILE *fp = fopen(fname,"w");
  fprintf(fp,"%lu\n",seed);
  fclose(fp);

  MATRIX *ymat = fMRIarrayToMatrix(yarray, mask, arraymajor);
  int nvox = ymat->cols;
  printf("%f\n",ymat->rptr[1][1]);
  if(ymatfile) MatrixWriteTxt(ymatfile, ymat);

  GLMMAT *glm = (GLMMAT *) calloc(sizeof(GLMMAT),1);

  // Create design
  if(fsgdfile){
    char gd2mtx_method[1000];
    sprintf(gd2mtx_method,"DODS");
    FSGD *fsgd = gdfRead(fsgdfile,gd2mtx_method,0);
    if(fsgd==NULL) exit(1);
    glm->X = gdfMatrix(fsgd,fsgd->gd2mtx_method,NULL);
    if(glm->X==NULL) exit(1);
    if(glm->X->rows != ymat->rows){
      printf("ERROR: fsgd %s has %d rows, expecting %d\n",fsgdfile,glm->X->rows,ymat->rows);
      exit(1);
    }
    if(fsgd->nContrasts == 0){
      printf("ERROR: fsgd %s has no contrasts\n",fsgdfile);
      exit(1);
    }
    if(fsgd->nContrasts > 1){
      printf("INFO: fsgd %s has %d contrasts, permutation testing number %d\n",fsgdfile,fsgd->nContrasts,cpermno);
    }
    glm->C[0] = MatrixCopy(fsgd->C[cpermno],NULL);
    glm->Cname[0] = strcpyalloc(fsgd->ContrastName[cpermno]);
    glm->ncontrasts = 1;
  } 
  else {
    if(DoOSGM == 1){
      glm->X = MatrixConstVal(1.0,ymat->rows,1,NULL);
      glm->ncontrasts = 1;
      glm->Cname[0] = strcpyalloc("osgm");
      glm->C[0] = MatrixConstVal(1.0,1,1,NULL);
    } else {
      glm->X = MatrixAlloc(ymat->rows,2,MATRIX_REAL);
      for(int n=0; n < ymat->rows; n++){
	if(n < ymat->rows/2) glm->X->rptr[n+1][1] = 1;
	else                 glm->X->rptr[n+1][2] = 1;
      }
      glm->ncontrasts = 2;
      glm->Cname[0] = strcpyalloc("tsgd");
      glm->C[0] = MatrixAlloc(1,2,MATRIX_REAL);
      glm->C[0]->rptr[1][1] = +1;
      glm->C[0]->rptr[1][2] = -1;
      glm->Cname[1] = strcpyalloc("mean");
      glm->C[1] = MatrixAlloc(1,2,MATRIX_REAL);
      glm->C[1]->rptr[1][1] = 1.0;
      glm->C[1]->rptr[1][2] = 1.0;
    }
  }
  for(int c=0; c < glm->ncontrasts; c++){
    printf("contrast %d %s\n",c,glm->Cname[c]);
    char condir[1000];
    sprintf(condir,"%s/%s",outdir,glm->Cname[c]);
    err = mkdir(condir,0777);
    if(err != 0 && errno != EEXIST) {
      printf("ERROR: creating directory %s\n",outdir);
      perror(NULL);
      return(1);
    }
  }

  GLMcMatrices(glm);
  GLMallocY(glm);
  GLMxMatrices(glm);

  // Allocate matrix arrays to hold results for each contrast
  std::vector<MATRIX *> pmatarray;
  std::vector<MATRIX *> gammaarray;
  for(int c=0; c < glm->ncontrasts; c++){
    MATRIX *gamma = MatrixAlloc(1,nvox,MATRIX_REAL); // univariate only
    gammaarray.push_back(gamma);
    MATRIX *pmat = MatrixAlloc(1,nvox,MATRIX_REAL);
    pmatarray.push_back(pmat);
  }

  // Do the fit at each voxel
  MATRIX *emat = MatrixAlloc(ymat->rows,ymat->cols,MATRIX_REAL);
  for(int vox = 0; vox < nvox; vox++){
    for(int f=0; f < ymat->rows; f++) glm->y->rptr[f+1][1] = ymat->rptr[f+1][vox+1];
    GLMfit(glm);
    GLMtest(glm);
    for(int f=0; f < ymat->rows; f++) emat->rptr[f+1][vox+1] = glm->eres->rptr[f+1][1];
    // pack the gamma and pvalue into the arrays
    for(int c=0; c < glm->ncontrasts; c++){
      gammaarray[c]->rptr[1][vox+1] = glm->gamma[c]->rptr[1][1];//univarate only
      pmatarray[c]->rptr[1][vox+1]  = glm->p[c];
    }
  }

  // Convert to MRI, cluster, Write out results
  MRI *pvol  = MRIcloneBySpace(yarray[0], MRI_FLOAT,nvariates);
  MRI *gamma = MRIcloneBySpace(yarray[0], MRI_FLOAT,nvariates);
  std::vector<std::vector<double>> clustersizes;
  for(int c=0; c < glm->ncontrasts; c++){
    printf("contrast %d %s\n",c,glm->Cname[c]);
    char condir[1000];
    sprintf(condir,"%s/%s",outdir,glm->Cname[c]);
    err = fMRIfromMatrix(gammaarray[c], gamma, mask);
    if(err) exit(1);

    // Write out gamma
    sprintf(fname,"%s/gamma.nii.gz",condir);
    MRIwrite(gamma,fname);
    int vno = 0;
    printf("%d  gamma %g %g\n",vno,MRIgetVoxVal(gamma,vno,0,0,0),MRIgetVoxVal(gamma,vno,0,0,1));

    // Write out sig (convert p-value to log10() and sign) and create link
    err = fMRIfromMatrix(pmatarray[c], pvol, mask);
    if(err) exit(1);
    MRI *sigvol = MRIlog10(pvol,mask,NULL,1);
    printf("%d  %g %g   %g %g\n",vno,MRIgetVoxVal(gamma,vno,0,0,0),MRIgetVoxVal(sigvol,vno,0,0,0),
	                             MRIgetVoxVal(gamma,vno,0,0,1),MRIgetVoxVal(sigvol,vno,0,0,1));
    MRIsetSign(sigvol,gamma,-1);
    printf("%d  %g %g   %g %g\n",vno,MRIgetVoxVal(gamma,vno,0,0,0),MRIgetVoxVal(sigvol,vno,0,0,0),
	                             MRIgetVoxVal(gamma,vno,0,0,1),MRIgetVoxVal(sigvol,vno,0,0,1));
    sprintf(fname,"%s/sig.nii.gz",condir);
    MRIwrite(sigvol,fname);
    char link[2000];
    sprintf(link,"%s.sig.nii.gz",glm->Cname[c]);
    makelocallink(fname,link,1);

    // Cluster
    SpatTempCluster stc;
    if(surf){
      stc.surf = surf;
      stc.topo = 2;
    }
    else {
      stc.topo = 1;
      stc.nbrtype = nbrtype;
    }
    printf("th %g %g %d\n",thmin,thmax,thsign);
    stc.GetBinMask(sigvol,thmin,thmax,thsign,mask);
    stc.Clusterize();
    stc.SortClusters(sigvol);
    sprintf(fname,"%s/ocn.nii.gz",condir);
    MRIwrite(stc.cnomap,fname);
    sprintf(fname,"%s/ctab.ocn",condir);
    CTABwriteFileASCII(stc.cnomap->ct,fname);
    sprintf(fname,"%s/clusters.dat",condir);
    stc.WriteClusterSum(fname,sigvol);
    sprintf(fname,"%s/clusters.json",condir);
    stc.WritePointSet(fname,sigvol);
    double maxcluster = stc.MaxClusterSize();
    sprintf(fname,"%s/max.cluster.dat",condir);
    fp = fopen(fname,"w");
    fprintf(fp,"%lf\n",maxcluster);
    fclose(fp);
    clustersizes.push_back(stc.GetClusterSizes());
    MRIfree(&sigvol);
  } // contrast
  MRIfree(&gamma);

  // Permutation loop =========================================================
  if(nperm > 0){
    printf("\n\nStarting perm loop %d\n",nperm);

    // Delete/empty perm.cluster-size.dat
    for(int c=0; c < glm->ncontrasts; c++){
      sprintf(fname,"%s/%s/perm.cluster-size.dat",outdir,glm->Cname[c]);
      FILE *fp = fopen(fname,"w");
      fclose(fp);
    }

    std::vector<std::vector<double>> permcsize(glm->ncontrasts,std::vector<double>(nperm,0));
#ifdef HAVE_OPENMP
  #pragma omp parallel for 
#endif
    for(int n = 0; n < nperm; n++){
      printf("%d t=%lf --------------\n",n,mytimer.seconds());fflush(stdout);
      //printf("VMPC1 %d\n",GetVmPeak()); fflush(stdout);
      GLMMAT *glmn = GLMdeepCopy(glm);
      MatrixRandPermRows(glmn->X, 2, seed+n); 
      GLMxMatrices(glmn);
      GLMcMatrices(glmn);
      // Allocate sigmat array
      std::vector<MATRIX *> sigmatarray;
      std::vector<MATRIX *> gammaarray, pmatarray;
      for(int c=0; c < glmn->ncontrasts; c++){
	MATRIX *sigmat = MatrixAlloc(1,nvox,MATRIX_REAL);
	sigmatarray.push_back(sigmat);
	MATRIX *pmat = MatrixAlloc(1,nvox,MATRIX_REAL);
	pmatarray.push_back(pmat);
	MATRIX *gamma = MatrixAlloc(1,nvox,MATRIX_REAL); // univariate only
	gammaarray.push_back(gamma);
      }
      // Do the fit at each voxel (use ymat or emat)
      for(int vox = 0; vox < nvox; vox++){
	for(int f=0; f < ymat->rows; f++) glmn->y->rptr[f+1][1] = ymat->rptr[f+1][vox+1];
	GLMfit(glmn);
	GLMtest(glmn);
	// pack the sig into an array for each contrast. 

	// Note: Currently, perm is giving clusters that are a little
	// too large and causing the stats to be conserv. Need to find
	// places where the perm is different than the test stat. test
	// stat uses MRIlog10() and MRIsetSign() (diff here). Maybe
	// the perm is not random enough (unlikely). It got worse as
	// the ico increased. Something specifically to do with
	// surfaces since vol *seems* ok?
	for(int c=0; c < glmn->ncontrasts; c++){
	  sigmatarray[c]->rptr[1][vox+1] = -log10(glmn->p[c]);
	  if(glmn->gamma[c]->rptr[1][1] < 0) sigmatarray[c]->rptr[1][vox+1] *= -1;
	  pmatarray[c]->rptr[1][vox+1]  = glmn->p[c];
	  gammaarray[c]->rptr[1][vox+1] = glmn->gamma[c]->rptr[1][1];//univarate only
	} //contrast
      }// vox


      MRI *pvol = MRIcloneBySpace(yarray[0], MRI_FLOAT,nvariates);
      MRI *sig = MRIcloneBySpace(yarray[0], MRI_FLOAT,nvariates);
      MRI *gamma = MRIcloneBySpace(yarray[0], MRI_FLOAT,nvariates);
      //MRI *pvoln  = MRIcloneBySpace(yarray[0], MRI_FLOAT,nvariates);
      fMRIfromMatrix(gammaarray[0], gamma, mask);
      fMRIfromMatrix(pmatarray[0], pvol, mask);//0=hack
      sig = MRIlog10(pvol,mask,NULL,1);
      MRIsetSign(sig,gamma,-1);
      //MRIwrite(sig,"junkp.nii.gz"); exit(1);
      MRIfree(&gamma);
      MRIfree(&pvol);

      for(int c=0; c < glmn->ncontrasts; c++){
	// Map the sig array into a volume
	//--err = fMRIfromMatrix(sigmatarray[c], sig, mask);
	//--if(err) exit(1);
	//MRIwrite(sig,"junkp.nii.gz");exit(1);
	// Cluster
	SpatTempCluster stcn;
	stcn.GetCtab=0;
	if(surf){
	  stcn.surf = surf;
	  stcn.topo = 2;
	}
	else {
	  stcn.topo = 1;
	  stcn.nbrtype = nbrtype;
	}
	stcn.GetBinMask(sig,thmin,thmax,thsign,mask);
	//printf("  th %g %g %d %d\n",thmin,thmax,thsign,nbrtype);
	stcn.Clusterize();
	double maxcluster = stcn.MaxClusterSize();
	permcsize[c][n] = maxcluster;
	printf("permn = %d maxc=%g\n",n,maxcluster);
	char condir[1000];
	sprintf(condir,"%s/%s",outdir,glm->Cname[c]);
	//sprintf(fname,"%s/mask.%04d.nii.gz",condir,n);
	//MRIwrite(mask,fname);
	sprintf(fname,"%s/perm.cluster-size.dat",condir);
	FILE *fp = fopen(fname,"a");
	fprintf(fp,"%d %g\n",n,maxcluster);fflush(fp);
	fclose(fp);
	//sprintf(fname,"%s/binmask.%04d.nii.gz",condir,n);
	//MRIwrite(stcn.binmask,fname);
	if(stcn.binmask) MRIfree(&stcn.binmask);
	if(n==0){
	  sprintf(fname,"%s/perm.sig.%04d.nii.gz",condir,n);
	  MRIwrite(sig,fname);
	  sprintf(fname,"%s/perm.clusters.%04d.dat",condir,n);
	  printf("fname=%s\n",fname);
	  FILE *fpn = fopen(fname,"w");
	  stcn.PrintClusterSum(fpn);
	  fclose(fpn);
	  sprintf(fname,"%s/ocn.%04d.nii.gz",condir,n);
	  MRIwrite(stcn.cnomap,fname);
	  //exit(1);
	}
      } // contrast
      for(int c=0; c < glmn->ncontrasts; c++) MatrixFree(&sigmatarray[c]);
      for(int c=0; c < glmn->ncontrasts; c++) MatrixFree(&pmatarray[c]);
      for(int c=0; c < glmn->ncontrasts; c++) MatrixFree(&gammaarray[c]);
      GLMfree(&glmn);
      MRIfree(&sig);
    } // for(nperm) ============================
    printf("done perm ============================\n\n");

    for(int c=0; c < glm->ncontrasts; c++){
      char condir[1000];
      sprintf(condir,"%s/%s",outdir,glm->Cname[c]);

      //Read in previously generated sig and redo clustering
      SpatTempCluster stc;      
      stc.surf = surf;
      stc.topo = 1; if(surf) stc.topo = 2;
      stc.nbrtype = nbrtype;
      sprintf(fname,"%s/sig.nii.gz",condir);
      MRI *sigvol0 = MRIread(fname);
      stc.GetBinMask(sigvol0,thmin,thmax,thsign,mask);
      stc.Clusterize();
      stc.SortClusters(sigvol0);

      std::vector<double> csizelist = clustersizes[c];
      std::vector<double> permcsizelist = permcsize[c];
      sprintf(fname,"%s/cluster.pval.dat",condir);
      FILE *fp = fopen(fname,"w");
      fprintf(fp,"Cno CSize p-value\n");
      double pmin = 1;
      for(int cno=0; cno < csizelist.size(); cno++){
	// Take average of > and >= to match volcluster
	int count=0,count2=0;
	for(int n=0; n < nperm; n++){
	  if(permcsizelist[n] > csizelist[cno])  count++; 
	  if(permcsizelist[n] >= csizelist[cno]) count2++;
	}
	double cpval = (double)(count+count2)/(2*nperm);
	cpval = (double)(count2+1)/(nperm+1); // Tom Nichols
	fprintf(fp,"%3d %8.1lf %20.18lf\n",cno+1,csizelist[cno],cpval);
	fflush(fp);
	if(cpval < pmin) pmin = cpval;
	stc.ClusterList[cno].pvalue = cpval;
      }
      fclose(fp);
      sprintf(fname,"%s/cluster.pval.min.dat",condir);
      fp = fopen(fname,"w");
      fprintf(fp,"%20.18lf\n",pmin);
      fclose(fp);
      sprintf(fname,"%s/clusters.perm.dat",condir);
      printf("writing %s\n",fname);
      stc.WriteClusterSum(fname,sigvol0);
      sprintf(fname,"%s/clusters.perm.json",condir);
      stc.WritePointSet(fname,sigvol0);
      MRIfree(&sigvol0);
    }

  } // if(nperm)


  printf("#VMPC# mri_mvglmfit VmPeak  %d\n",GetVmPeak());
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
    else if(!strcasecmp(option, "--o") || !strcasecmp(option, "--glmdir")) {
      if(nargc < 1) CMDargNErr(option,1);
      outdir = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--i") || !strcasecmp(option, "--y")) {
      if(nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while(CMDnthIsArg(nargc, pargv, nargsused)){
	yflist.push_back(pargv[nargsused]);
        nargsused++;
      }
      arraymajor=1;
    }
    else if(!strcasecmp(option, "--i2") || !strcasecmp(option, "--y2")) {
      if(nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while(CMDnthIsArg(nargc, pargv, nargsused)){
	yflist.push_back(pargv[nargsused]);
        nargsused++;
      }
      arraymajor=2;
    }
    else if(!strcasecmp(option, "--ymatfile")) {
      if(nargc < 1) CMDargNErr(option,1);
      ymatfile = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--th")) {
      if(nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while(CMDnthIsArg(nargc, pargv, nargsused)){
	double thresh;
	sscanf(pargv[nargsused],"%lf",&thresh);
	threshlist.push_back(thresh);
	thmin = thresh;
        nargsused ++;
      }
    }
    else if(!strcasecmp(option, "--fsgd")) {
      if(nargc < 1) CMDargNErr(option,1);
      fsgdfile = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--osgm")) DoOSGM = 1;
    else if(!strcasecmp(option, "--tsgd")) DoOSGM = 2;
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
    else if(!strcasecmp(option, "--nperm")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nperm);
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--cpermno")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&cpermno);
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--face")) nbrtype = 1;
    else if(!strcasecmp(option, "--edge")) nbrtype = 2;
    else if(!strcasecmp(option, "--corner")) nbrtype = 3;
    else if(!strcasecmp(option, "--seed")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lu",&seed);
      nargsused = 1;
    }
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
  printf("   --y  var1 <var2...> : each file is a stack of subjects from a time point\n");
  printf("   --y2 var1 <var2...> : each file is a stack of time points from a subject\n");
  printf("   --mask mask\n");
  printf("   --surf surffile : when input is a surface\n");
  printf("   --fsgd fsgdfile\n");
  printf("   --cpermno contrastno : test using contrastno from fsgd file\n");
  printf("   --osgm, --tsgm\n");
  printf("   --nperm nperm\n");
  printf("   --th threshold : cluster-forming threshold\n");
  printf("   --face : volume neighborhood defined by adjancent faces only\n");
  printf("   --edge : volume neighborhood defined by adjancent edge and faces only\n");
  printf("   --corner : volume neighborhood defined by adjancent edge and faces and corners (default)\n");
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
  if(yflist.size() == 0){
    printf("ERROR: must spec at least one input\n");
    exit(1);
  }
  if(seed == 0) seed = PDFtodSeed();

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
  fprintf(fp,"input   %s\n",yflist[0]);
  fprintf(fp,"seed   %lu\n",seed);
  fprintf(fp,"DoOSGM %d\n",DoOSGM);
  fprintf(fp,"arraymajor %d\n",arraymajor);
  fprintf(fp,"nperm %d\n",nperm);
  fprintf(fp,"thmin %lf\n",thmin);
  fprintf(fp,"thmax %lf\n",thmax);
  fprintf(fp,"thsign %d\n",thsign);
  fprintf(fp,"nbrtype %d\n",nbrtype);
  fprintf(fp,"threads %d\n",threads);
  if(maskfile) fprintf(fp,"mask  %s\n",maskfile);
  if(surffile) fprintf(fp,"mask  %s\n",surffile);
  return;
}

MATRIX *fMRIarrayToMatrix(std::vector<MRI*> mriArray, MRI *mask, int major)
{
  int narray = mriArray.size(), nvariates, nsubjects;
  MRI *v0 = mriArray[0];

  if(major == 1){
    nvariates = narray;
    nsubjects = v0->nframes;
  } else {
    nvariates = v0->nframes;
    nsubjects = narray;
  }

  int nvox=0;
  if(mask){
    if(MRIdimMismatch(v0, mask, 0)){
      printf("MRIarray2Mat(): mask dimension mismatch\n");
      return(NULL);
    }
    for(int s=0; s < mask->depth; s++){
      for(int r=0; r < mask->height; r++){
	for(int c=0; c < mask->width; c++){
	  if(MRIgetVoxVal(mask,c,r,s,0)<0.5) continue;
	  nvox++;
	}
      }
    }
    printf("nmask = %d\n",nvox);
  } else {
    nvox = v0->width * v0->height * v0->depth;
  }
  printf("Volume size %d %d %d %d major=%d\n", 
	 v0->width,v0->height,v0->depth,v0->nframes,major);
  nvox *= nvariates;
  printf("major %d nvarites = %d nvox = %d\n",major,nvariates,nvox);

  MATRIX *m = MatrixAlloc(nsubjects,nvox,MATRIX_REAL);
  printf("matrix size: %d %d\n",m->rows,m->cols);

  // Create the matrix where the rows are the dependent variable
  // (frame/subject) and the columns are the spatial variables and the
  // array variable.  Eg, if the array variable is time, there will be
  // nrows*ncols*nslices*narray columns in the matrix arranged with
  // narray is the slowest variable. If comparing to matlab, permute
  // the first and second variables when reading the volume in with
  // MRIread().  Order compatible with fMRIfromMatrix() Can't
  // parallelize. This allows for each MRI input file to be a subject
  // stack from a given timepoint as would be output from
  // isxconcat-sess for that time point.
  int k=0;
  for(int q=0; q < nvariates; q++){ // variate dim
    for(int s=0; s < v0->depth; s++){
      for(int r=0; r < v0->height; r++){
	for(int c=0; c < v0->width; c++){
	  if(mask && MRIgetVoxVal(mask,c,r,s,0)<0.5) continue;
	  for(int f=0; f < nsubjects; f++){// subject/dep variable
	    double val;
	    if(major == 1){
	      MRI *v = mriArray[q];
	      val = MRIgetVoxVal(v,c,r,s,f);
	    } else {
	      MRI *v = mriArray[f];
	      val = MRIgetVoxVal(v,c,r,s,q);
	    }
	    m->rptr[f+1][k+1] = val;
	  }
	  k++;
	}
      }
    }
  }

  return(m);
}

#if 0
// Smooth across both space and array. mriArray [ncols, nrows, nslices, nsubj, nvar]
std::vector<MRI*> MRIarrayGaussianSmooth(std::vector<MRI*>mriArray,std::vector<double>gstd,MRI *mask, double arrayRes)
{
  int nvariates = mriArray.size();
  std::vector<MRI*> mriArraySm;
  //MRImaskedGaussianSmooth(InVals, mritmp, ingstd, InVals);
  //MRIgaussianSmoothNI(InVals, ingstdc, ingstdr, ingstds, InVals);
  // Can't do non-iso smoothing with mask?
  return(mriArraySm);
}
#endif


