/**
 * @brief Uses the 4 surfaces of a scan to construct a mask volume
 *
 * Uses the 4 surfaces of a scan to construct a mask volume showing the
 * position of each voxel with respect to the surfaces - GM, WM, LH or RH.
 *
 * Uses the MRISOBBTree algorithm
 */
/*
 * Original Author: Krish Subramaniam
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

// STL
#include <string>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <vector>

#include "MRISOBBTree.h"
#include "MRISdistancefield.h"
#include "fastmarching.h"
#include "cmd_line_interface.h"

// FS


#include "fsenv.h"
#include "mrisurf.h"
#include "mri.h"
#include "error.h"
#include "cma.h"
#include "diag.h"
#include "macros.h"
#include "timer.h"
#include "gca.h"
#include "version.h"

#include "romp_support.h"
#undef private

;
const char *Progname;

typedef Math::Point<int> Pointd;

// static function declarations
// forward declaration
struct IoParams;
using namespace std;

class IoError : public std::exception
{
public:
  IoError(const std::string& what)
    : m_what(what)
  {}
  virtual ~IoError() throw()
  {}
  using std::exception::what;
  virtual const char* what()
  {
    return m_what.c_str();
  }
private:
  std::string m_what;
};

/*
  loads all the input files
  the return value is the output path for the volume files
*/
std::string
operator/(const std::string& a,
          const std::string& b)
{
  std::string s(a);
  s += "/" + b;
  return s;
}

std::string
LoadInputFiles(const IoParams& params,
               MRI*& mriInput,
               MRIS*& surfLeftWhite,
               MRIS*& surfLeftPial,
               MRIS*& surfRightWhite,
               MRIS*& surfRightPial);

MRI* ComputeSurfaceDistanceFunction
(MRIS* mris, //input surface
 MRI* mriInOut, //output MRI structure
 float resolution,
 int indexSurf,
 bool savedistance,
 std::string outputPath, std::string outRoot);

MRI* CreateHemiMask(MRI* dpial, MRI* dwhite,
                    const unsigned char lblWhite,
                    const unsigned char lblRibbon,
                    const unsigned char lblBackground);

MRI* CombineMasks(MRI* maskOne,
                  MRI* maskTwo,
                  const unsigned char lblBackground);

//  return a binary mask
MRI* FilterLabel(MRI* initialMask,
                 const unsigned char lbl);

std::string InitVolmask(const IoParams& params, MRI*& mriTemplate, std::vector<MRI*>& surfDistVector, std::vector<MRIS*>& inSurfVector);
void Cleanup(std::vector<MRI*>& surfDistVector, std::vector<MRIS*>& inSurfVector);
void SaveSurfaceDistance2(std::vector<MRI*>& surfDistVector, std::string outputPath, std::string outRoot);
void SaveSurfaceDistances(MRI *surfDist, int threadid, const char *hemi, const char *surfname, std::string outputPath, std::string outRoot);
void SaveHemispheresRibbon(const char *hemi, MRI *mriTemplate, MRI *maskHemi, const unsigned char labelRibbon, std::string outputPath, std::string outRoot);

/*
  IO structure

  The user can either specify an input subject (thus implicitly
  setting the names of the in-out files)
  or use an advanced mode and manually specify the files
*/
struct IoParams
{
  typedef std::string StringType;

  StringType  subject;
  StringType  subjectsDir;

  StringType  templateMri;
  StringType  surfWhiteRoot;
  StringType  surfPialRoot;

  StringType  asegName;

  StringType  outRoot;

  unsigned char labelLeftWhite;
  unsigned char labelLeftRibbon;
  unsigned char labelRightWhite;
  unsigned char labelRightRibbon;
  unsigned char labelBackground;
  int DoLH, DoRH;
  bool bLHOnly, bRHOnly;
  bool bParallel;
  bool bVerbose;

  float capValue;

  bool bSaveDistance;
  bool bSaveRibbon;
  bool bEditAseg;

  IoParams();
  void parse(int ac, char* av[]);
};

int
main(int ac, char* av[])
{
  // first, handle stuff for --version and --all-info args
  int nargs = 0, msec;
  Timer then ;

  then.reset() ;
  nargs = handleVersionOption(ac, av, "mris_volmask");
  if (nargs && ac - nargs == 1)
    exit (0);
  ac -= nargs;

  // parse command-line
  IoParams params;
  try{
    params.parse(ac,av);
  }
  catch (std::exception& excp){
    std::cerr << " Exception caught while parsing the command-line\n"<< excp.what() << std::endl;
    exit(1);
  }

  if (params.bVerbose) Gdiag |= DIAG_VERBOSE;
  if(params.bLHOnly){params.DoLH=1;params.DoRH=0;}
  if(params.bRHOnly){params.DoLH=0;params.DoRH=1;}
  //printf("lhrhonly %d %d %d %d\n",params.bLHOnly,params.bRHOnly,params.DoLH,params.DoRH);
  int idx_offset = 0, nSurfDists = 2;
  if (!params.DoLH)
    idx_offset = 2;
  if (params.DoLH && params.DoRH)
    nSurfDists = 4;
  if (Gdiag & DIAG_VERBOSE)
    printf("[DEBUG] nSurfDists = %d, vector idx_offset = %d\n", nSurfDists, idx_offset);

  MRI *mriTemplate = NULL;
  std::vector<MRI*> surfDistVector;  
  std::vector<MRIS*> inSurfVector;

  //setenv("FS_MGZIO_USEVOXELBUFREAD", "1", 1);
  //setenv("FS_MGZIO_USEVOXELBUFWRITE", "1", 1);
  std::string outputPath = InitVolmask(params, mriTemplate, surfDistVector, inSurfVector);

  if (params.bEditAseg)
  {
    std::string subjDir = params.subjectsDir / params.subject;
    std::string pathSurf(subjDir / "surf"),
        pathMriOutput = outputPath / "aseg.ribbon.mgz";

    if (params.bRHOnly == 0)
    {
      printf("inserting LH into aseg...\n") ;
      insert_ribbon_into_aseg(mriTemplate, mriTemplate,
			      inSurfVector[0], inSurfVector[1], LEFT_HEMISPHERE) ;
    }
    if (params.bLHOnly == 0)
    {
      printf("inserting RH into aseg...\n") ;
      insert_ribbon_into_aseg(mriTemplate, mriTemplate,
			      inSurfVector[2], inSurfVector[3], RIGHT_HEMISPHERE);
    }
    printf("writing output to %s\n",
           (const_cast<char*>( pathMriOutput.c_str() )));
    mriTemplate->ct = CTABreadDefault();
    MRIwrite(mriTemplate,  (const_cast<char*>( pathMriOutput.c_str() ))) ;
    msec = then.milliseconds() ;
    fprintf(stderr, "mris_volmask took %2.2f seconds\n", (float)msec/1000.0f);   ///(1000.0f*60.0f));
    exit(0) ;
  }

#ifdef _OPENMP
  if (params.bParallel){
    int nthreads = nSurfDists;
    printf("Computing surface distance in parallel (%d threads)\n", nthreads);
    omp_set_nested(1); // Enable nested parallelism
    omp_set_num_threads(nthreads);
  }
  else{
    printf("Computing surface distance serially\n");
    omp_set_num_threads(1);
  }
#endif

  Timer then1;
  if (Gdiag & DIAG_VERBOSE)
    then1.reset();

  #ifdef HAVE_OPENMP
  #pragma omp parallel for 
  #endif
  for (int n = 0; n < nSurfDists; n++)
      ComputeSurfaceDistanceFunction(inSurfVector[n+idx_offset],
				     surfDistVector[n+idx_offset],
				     params.capValue,
				     n+idx_offset,
                                     params.bSaveDistance,
                                     outputPath, params.outRoot);
  if (Gdiag & DIAG_VERBOSE)
    fprintf(stderr, "[TIMER] %d ComputeSurfaceDistanceFunction() (%2.2f seconds)\n", nSurfDists, (float)then1.milliseconds()/1000.0f);

  MRI *maskLeftHemi = NULL, *maskRightHemi=NULL;
  #ifdef HAVE_OPENMP
  #pragma omp parallel for 
  #endif
  for (int n = 0; n < nSurfDists/2; n++) {
    if (n == 0 && params.DoLH) {
      maskLeftHemi = CreateHemiMask(surfDistVector[1], surfDistVector[0],
				      params.labelLeftWhite,
				      params.labelLeftRibbon,
				      params.labelBackground);
      if (params.bSaveRibbon)
        SaveHemispheresRibbon("lh", mriTemplate, maskLeftHemi, params.labelLeftRibbon, outputPath, params.outRoot);
    }
    if (n == 1 && params.DoRH) {
      maskRightHemi = CreateHemiMask(surfDistVector[3], surfDistVector[2],
				       params.labelRightWhite,
				       params.labelRightRibbon,
                                       params.labelBackground);
      if (params.bSaveRibbon)
        SaveHemispheresRibbon("rh", mriTemplate, maskRightHemi, params.labelRightRibbon, outputPath, params.outRoot);	 
    }
  }

  /*  finally combine the two created masks -- need to resolve overlap  */

  if (Gdiag & DIAG_VERBOSE)
    then1.reset();
  MRI* finalMask = NULL;
  if(params.DoLH && params.DoRH)
   finalMask= CombineMasks(maskLeftHemi, maskRightHemi,params.labelBackground);
  else if(params.DoLH) finalMask = maskLeftHemi;
  else if(params.DoRH) finalMask = maskRightHemi;
  MRIcopyHeader( mriTemplate, finalMask);
  finalMask->ct = CTABreadDefault();
  if (Gdiag & DIAG_VERBOSE)
    fprintf(stderr, "[TIMER] finalMask (%2.2f msec)\n", (float)then1.milliseconds());

  // write final mask
  std::cout << "writing cortical ribbon mask " << const_cast<char*>( (outputPath / (params.outRoot +".mgz")).c_str() ) << endl;
  MRIwrite( finalMask,const_cast<char*>( (outputPath / (params.outRoot +".mgz")).c_str() ));
  // sanity-check: make sure location 0,0,0 is background (not brain)
  if ( MRIgetVoxVal(finalMask,0,0,0,0) != 0 ) {
    cerr << "ERROR: ribbon has non-zero value at location 0,0,0" << endl;
    exit(1);
  }

  // cleanup
  MRIfree(&mriTemplate);
  Cleanup(surfDistVector, inSurfVector);
  if (maskLeftHemi != NULL)
    MRIfree(&maskLeftHemi);
  if (maskRightHemi != NULL)
    MRIfree(&maskRightHemi);
  
  msec = then.milliseconds() ;
  fprintf(stderr, "mris_volmask took %2.2f seconds\n", (float)msec/1000.0f);   ///(1000.0f*60.0f));
  return 0;
}


IoParams::IoParams()
{
  labelBackground = 0;

  labelLeftWhite = 20;
  labelRightWhite = 120;

  labelLeftRibbon = 10;
  labelRightRibbon = 110;

  capValue = 3;
  bSaveDistance = false;
  bEditAseg = false ;
  bSaveRibbon = false;
  bLHOnly = false;
  bRHOnly = false;
  bParallel = false;
  DoLH = 1;
  DoRH = 1;

  outRoot = "ribbon";
  asegName = "aseg";
  surfWhiteRoot = "white";
  surfPialRoot = "pial";

  char *sd = FSENVgetSUBJECTS_DIR();
  if(NULL == sd) subjectsDir = "";
  else           subjectsDir = sd;
}

void
IoParams::parse(int ac, char* av[])
{
  std::string sl = "left_";
  std::string sr = "right_";
  std::string srib = "ribbon";
  std::string sw = "white";
  std::string ssurf = "surf_";
  std::string slbl = "label_";

  int iLeftWhite(labelLeftWhite),
      iLeftRibbon(labelLeftRibbon),
      iRightWhite(labelRightWhite),
      iRightRibbon(labelRightRibbon),
      iBackground(labelBackground);

  std::string strUse =  "surface root name (i.e. <subject>/surf/$hemi.<NAME>";
  CCmdLineInterface interface(av[0]);
  bool showHelp(false);

  interface.AddOptionBool
  ( "help", &showHelp, "display help message");
  interface.AddOptionBool
  ( "usage", &showHelp, "display help message");
  interface.AddOptionString
  ( (ssurf+sw).c_str(), &surfWhiteRoot,
    (strUse + " - default value is white").c_str()
  );
  interface.AddOptionString
  ( (ssurf+"pial").c_str(), &surfPialRoot,
    (strUse + " - default value is pial").c_str()
  );
  interface.AddOptionString
  ( "aseg_name", &asegName,
    "default value is aseg, allows name override"
  );
  interface.AddOptionString
  ( "out_root", &outRoot,
    "default value is ribbon - output will then be mri/ribbon.mgz and "
    "mri/lh.ribbon.mgz and mri/rh.ribbon.mgz "
    "(last 2 if -save_ribbon is used)"
  );
  interface.AddOptionInt
  ( (slbl+"background").c_str(), &iBackground,
    "override default value for background label value (0)"
  );
  interface.AddOptionInt
  ( (slbl+sl+sw).c_str(), &iLeftWhite,
    "override default value for left white matter label - 20"
  );
  interface.AddOptionInt
  ( (slbl+sl+srib).c_str(), &iLeftRibbon,
    "override default value for left ribbon label - 10"
  );
  interface.AddOptionInt
  ( (slbl+sr+sw).c_str(), &iRightWhite,
    "override default value for right white matter label - 120"
  );
  interface.AddOptionInt
  ( (slbl+sr+srib).c_str(), &iRightRibbon,
    "override default value for right ribbon label - 110"
  );
  interface.AddOptionFloat
  ( "cap_distance", &capValue,
    (char*)"maximum distance up to which the signed distance function "
    "computation is accurate"
  );
  interface.AddOptionBool
  ( "save_distance", &bSaveDistance,
    "option to save the signed distance function as ?h.dwhite.mgz "
    "?h.dpial.mgz in the mri directory"
  );
  interface.AddOptionBool( "lh-only", &bLHOnly,"only analyze the left hemi");
  interface.AddOptionBool( "rh-only", &bRHOnly,"only analyze the right hemi");
  interface.AddOptionBool( "parallel", &bParallel,"run hemis in parallel");
  interface.AddOptionBool( "verbose",  &bVerbose, "output debug info");
  interface.AddOptionBool
  ( "edit_aseg", &bEditAseg,
    "option to edit the aseg using the ribbons and save to "
    "aseg.ribbon.mgz in the mri directory"
  );
  interface.AddOptionBool
  ( "save_ribbon", &bSaveRibbon,
    "option to save just the ribbon for the hemispheres - "
    "in the format ?h.ribbon.mgz"
  );
  interface.AddOptionString
  ( "sd", &subjectsDir,
    "option to specify SUBJECTS_DIR, default is read from enviro"
  );
  interface.AddIoItem(&subject, " subject (required param!)");

  // if ac == 0, then print complete help
  if ( ac == 1 )
  {
    interface.PrintHelp();
    exit(0);
  }

  interface.Parse(ac,av);
  if ( showHelp )
  {
    exit(0);
  }

  labelLeftWhite = (unsigned char)(iLeftWhite);
  labelLeftRibbon = (unsigned char)(iLeftRibbon);
  labelRightWhite = (unsigned char)(iRightWhite);
  labelRightRibbon = (unsigned char)(iRightRibbon);
  labelBackground = (unsigned char)(iBackground);
}


std::string
LoadInputFiles(const IoParams& params,
               MRI*& mriTemplate,
               MRIS*& surfLeftWhite,
               MRIS*& surfLeftPial,
               MRIS*& surfRightWhite,
               MRIS*& surfRightPial)
{
  // determine the mode of the application and infer the input file names

  // declare path objects
  std::string pathSurfLeftWhite,
      pathSurfLeftPial,
      pathSurfRightWhite,
      pathSurfRightPial,
      pathMriInput,
      pathOutput;

  if ( params.subjectsDir.empty() )
  {
    cerr << "SUBJECTS_DIR not found. Use --sd <dir>, or set SUBJECTS_DIR"
         << endl;
    exit(1);
  }
  else
  {
    cout << "SUBJECTS_DIR is " << params.subjectsDir << endl;
  }

  if ( !params.subject.empty() ) // application is in subject-mode
  {
    std::string subjDir = params.subjectsDir / params.subject;
    std::string pathSurf(subjDir / "surf");
    pathSurfLeftWhite = pathSurf / "lh." +  params.surfWhiteRoot ;
    pathSurfLeftPial = pathSurf / "lh." + params.surfPialRoot;
    pathSurfRightWhite = pathSurf / "rh." +  params.surfWhiteRoot;
    pathSurfRightPial = pathSurf / "rh." + params.surfPialRoot;
    pathMriInput = subjDir / "mri" / params.asegName + ".mgz";

    pathOutput = subjDir / "mri";
  }

  printf("loading input data...\n") ;

  // load actual files now
  if(params.DoLH){

    printf("Loading %s\n",pathSurfLeftWhite.c_str());
    surfLeftWhite = MRISread( const_cast<char*>( pathSurfLeftWhite.c_str() ));
    if ( !surfLeftWhite )
      throw IoError( std::string("failed to read left white surface ") + pathSurfLeftWhite );
    printf("Loading %s\n",pathSurfLeftPial.c_str());
    surfLeftPial  = MRISread( const_cast<char*>( pathSurfLeftPial.c_str() ));
    if ( !surfLeftPial )
      throw IoError( std::string("failed to read left pial surface ")+ pathSurfLeftPial );

    // Must verify that determinant is < 0
    MATRIX *m = MRIgetVoxelToRasXform(&surfLeftWhite->vg) ;
    if(MatrixDeterminant(m) > 0) {
      printf("INFO: lh surf vg vox2ras determinant is > 0 so reversing face order\n");
      MRISreverseFaceOrder(surfLeftWhite);
      MRISreverseFaceOrder(surfLeftPial);
    }
    MatrixFree(&m);
  }

  if(params.DoRH){
    printf("Loading %s\n",pathSurfRightWhite.c_str());
    surfRightWhite = MRISread( const_cast<char*>( pathSurfRightWhite.c_str() ));
    if ( !surfRightWhite )
      throw IoError( std::string("failed to read right white surface ")
		     + pathSurfRightWhite );
    printf("Loading %s\n",pathSurfRightPial.c_str());
    surfRightPial = MRISread( const_cast<char*>( pathSurfRightPial.c_str() ));
    if ( !surfRightPial )
      throw IoError( std::string("failed to read right pial surface ")
		     + pathSurfRightPial );
    // Must verify that determinant is < 0
    MATRIX *m = MRIgetVoxelToRasXform(&surfRightWhite->vg) ;
    if(MatrixDeterminant(m) > 0) {
      printf("INFO: rh surf vg vox2ras determinant is > 0 so reversing face order\n");
      MRISreverseFaceOrder(surfRightWhite);
      MRISreverseFaceOrder(surfRightPial);
    }
    MatrixFree(&m);
  }
    
  printf("Loading %s\n",pathMriInput.c_str());
  mriTemplate = MRIread( const_cast<char*>( pathMriInput.c_str() ));
  if ( !mriTemplate )
    throw IoError( std::string("failed to read template mri ")
                   + pathMriInput );

  return pathOutput;
}



MRI*
ComputeSurfaceDistanceFunction(MRIS* mris,
                               MRI* mri_distfield,
                               float thickness,
			       int indexSurf,
			       bool savedistance,
			       std::string outputPath, std::string outRoot)
{
  int tid = 0;
#ifdef HAVE_OPENMP
  tid = omp_get_thread_num();
#endif
  printf("[thread %d] computing distance to %s %s surface\n", tid, (indexSurf < 2) ? "lh" : "rh", (indexSurf%2 == 0) ? "white" : "pial");

  /******* use mri_distfield directly *******/
  MRI *mri_visited  = MRIcloneDifferentType(mri_distfield, MRI_INT);

  // Convert surface vertices coordinates from tkrRAS to vox space
  Math::ConvertSurfaceRASToVoxel(tid, mris, mri_distfield);

  Timer then1, then2, then3, then4;
  int then3_sum = 0, then4_sum = 0;
  int res0 = 0, resplusone = 0, resminusone = 0, resall = 0;

  if (Gdiag & DIAG_VERBOSE)
    then1.reset();

  // Find the distance field
  MRISDistanceField *distfield = new MRISDistanceField(mris, mri_distfield);
  distfield->SetMaxDistance(thickness);
  distfield->Generate(); //mri_distfield now has the distancefield

  if (Gdiag & DIAG_VERBOSE) {
    fprintf(stderr, "[TIMER] [thread %d] ComputeSurfaceDistanceFunction() MRISDistanceField->Generate() (%2.2f seconds)\n", tid, (float)then1.milliseconds()/1000.0f);

    char distfieldname[256] = {'\0'};
    sprintf(distfieldname, "/autofs/cluster/scratch_wednesday/yh887/bert.copy/mri/%s.%s.distfield.mgz", (indexSurf < 2) ? "lh" : "rh", (indexSurf%2 == 0) ? "white" : "pial");
    MRIwrite(mri_distfield, distfieldname);
    printf("[DEBUG] [thread %d] saved %s %s unsigned distance field: %s\n", tid, (indexSurf < 2) ? "lh" : "rh", (indexSurf%2 == 0) ? "white" : "pial", distfieldname);

    then2.reset();
  }
  
  // Construct the OBB Tree
  MRISOBBTree* OBBTree = new MRISOBBTree(mris);
  OBBTree->ConstructTree();

  if (Gdiag & DIAG_VERBOSE) {
    fprintf(stderr, "[TIMER] [thread %d] ComputeSurfaceDistanceFunction() MRISOBBTree->ConstructTree() (%2.2f seconds)\n", tid, (float)then2.milliseconds()/1000.0f);
    OBBTree->PrintSurfBoundingBox(tid, (indexSurf < 2) ? "lh" : "rh", (indexSurf%2 == 0) ? "white" : "pial");

    then2.reset();
  }
	
  std::queue<Pointd* > ptsqueue;
  // iterate through all the volume points
  // and apply sign
  for(int i=0; i< mri_distfield->width; i++)
  {
    //std::cerr << i <<" ";
    for(int j=0; j< mri_distfield->height; j++)
    {
      for(int k=0; k< mri_distfield->depth; k++)
      {
        if (Gdiag & DIAG_VERBOSE)
	  then3.reset();

        if ( MRIIvox(mri_visited, i, j, k ))
        {  
          continue;
        }
	int res = OBBTree->PointInclusionTest(i, j, k);
	
	if (Gdiag & DIAG_VERBOSE) {
	  resall++;
	  if (res == -1)     resminusone++;
	  else if (res == 1) resplusone++;
	  else if (res == 0) res0++;
	  else printf("[ERROR] ****** res is not -1, 1, or 0\n");

          then3_sum += then3.milliseconds();
	  then4.reset();
	}

        Pointd *pt = new Pointd;
        pt->v[0] = i;
        pt->v[1] = j;
        pt->v[2] = k;
        ptsqueue.push( pt );

        // First serve all the points in the queue before going to the next voxel
        while ( !ptsqueue.empty() )
        {
          // serve the front and pop it
          Pointd *p = ptsqueue.front();
          const int x = p->v[0];
          const int y = p->v[1];
          const int z = p->v[2];
          delete p;
          ptsqueue.pop();
	  
          if ( MRIIvox(mri_visited, x, y, z) )
          {
            continue;
          }
          MRIIvox(mri_visited, x, y, z) =  res;
          const float dist = MRIFvox(mri_distfield, x, y, z);
          MRIFvox(mri_distfield, x, y, z) =  dist*res;

          // mark its 6 neighbors if distance > 1 ( triangle inequality )
          if ( dist > 1 )
          {
            // left neighbor in x
            if ( x>0 && !MRIIvox(mri_visited, x-1, y, z))
            {
              Pointd *ptemp = new Pointd;
              ptemp->v[0]   = x - 1;
              ptemp->v[1]   = y;
              ptemp->v[2]   = z;
              ptsqueue.push(ptemp);
            }
            // bottom neighbor in y
            if ( y>0 && !MRIIvox(mri_visited, x, y-1, z))
            {
              Pointd *ptemp = new Pointd;
              ptemp->v[0]   = x;
              ptemp->v[1]   = y - 1;
              ptemp->v[2]   = z;
              ptsqueue.push(ptemp);
            }
            // front neighbor in z
            if ( z>0 && !MRIIvox(mri_visited, x, y, z-1))
            {
              Pointd *ptemp = new Pointd;
              ptemp->v[0]   = x;
              ptemp->v[1]   = y;
              ptemp->v[2]   = z - 1;
              ptsqueue.push(ptemp);
            }
            // right neighbor in x
            if ( x<mri_visited->width-1 && !MRIIvox(mri_visited, x+1, y, z))
            {
              Pointd *ptemp = new Pointd;
              ptemp->v[0]   = x + 1;
              ptemp->v[1]   = y;
              ptemp->v[2]   = z;
              ptsqueue.push(ptemp);
            }
            // top neighbor in y
            if ( y<mri_visited->height-1 && !MRIIvox(mri_visited, x, y+1, z))
            {
              Pointd *ptemp = new Pointd;
              ptemp->v[0]   = x;
              ptemp->v[1]   = y + 1;
              ptemp->v[2]   = z;
              ptsqueue.push(ptemp);
            }
            // back neighbor in z
            if ( z<mri_visited->depth-1 && !MRIIvox(mri_visited, x, y, z+1))
            {
              Pointd *ptemp = new Pointd;
              ptemp->v[0]   = x;
              ptemp->v[1]   = y;
              ptemp->v[2]   = z + 1;
              ptsqueue.push(ptemp);
            }
          }
        }

	if (Gdiag & DIAG_VERBOSE)
	  then4_sum += then4.milliseconds();
      }
    }
  }

  if (Gdiag & DIAG_VERBOSE) {
    fprintf(stderr, "[TIMER] [thread %d] ComputeSurfaceDistanceFunction() (for loops) apply sign to all the volume points (%2.2f seconds)\n",
                    tid, (float)then2.milliseconds()/1000.0f);

    fprintf(stderr, "[TIMER] [thread %d] ComputeSurfaceDistanceFunction() MRISOBBTree->PointInclusionTest() (%2.2f seconds) "
                    "points tested: %d, outside: %d, inside: %d, zeors: %d\n",
                    tid, (float)then3_sum/1000.0f, resall, resminusone, resplusone, res0);
    fprintf(stderr, "[TIMER] [thread %d] ComputeSurfaceDistanceFunction() (while loop) apply sign to voxels connected with dist > 1 (%2.2f seconds)\n", tid, (float)then4_sum/1000.0f);
  }

  if (savedistance)
    SaveSurfaceDistances(mri_distfield, tid, (indexSurf < 2) ? "lh" : "rh", (indexSurf%2 == 0) ? "white" : "pial", outputPath, outRoot);
    
  MRIfree(&mri_visited);
  delete OBBTree;
  delete distfield;
  return(mri_distfield);
}


// combine pial and white surface distance to create a mask for the hemi
// Must be outside of white and inside pial. Creates labels for WM and Ribbon.
MRI* CreateHemiMask(MRI* dpial,
                    MRI* dwhite,
                    const unsigned char lblWhite,
                    const unsigned char lblRibbon,
                    const unsigned char lblBackground)
{
  int tid = 0;
#ifdef HAVE_OPENMP
  tid = omp_get_thread_num();
#endif

  Timer then1;
  if (Gdiag & DIAG_VERBOSE)
    then1.reset();
  
  // allocate return volume
  MRI* mri = MRIalloc(dpial->width,
                      dpial->height,
                      dpial->depth,
                      MRI_UCHAR);

  #ifdef HAVE_OPENMP
  #pragma omp parallel for 
  #endif
  for (int z = 0; z < dpial->depth; ++z)
    for (int y = 0; y < dpial->height; ++y)
      for (int x = 0; x < dpial->width; ++x)
      {
        if ( MRIFvox(dwhite,x,y,z) > 0 )
        {
          MRIsetVoxVal(mri, x,y,z,0, lblWhite);
        }
        else if ( MRIFvox(dpial,x,y,z) > 0 )
        {
          MRIsetVoxVal(mri, x,y,z,0, lblRibbon);
        }
        else
        {
          MRIsetVoxVal(mri,x,y,z,0, lblBackground);
        }
      } // next x,y,z

  if (Gdiag & DIAG_VERBOSE)
    fprintf(stderr, "[TIMER] [thread %d] CreateHemiMask(white=%d, ribbon=%d, background=%d) (%2.2f msec)\n", tid, lblWhite, lblRibbon, lblBackground, (float)then1.milliseconds());
      
  return mri;
}

/*
  implicit assumption the masks do not overlap
  if this is the case, a message will be printed
*/
MRI* CombineMasks(MRI* maskOne,
                  MRI* maskTwo,
                  const unsigned char lblBackground)
{
  MRI* mri = MRIalloc(maskOne->width,
                      maskOne->height,
                      maskOne->depth,
                      MRI_UCHAR);
  unsigned int overlap = 0;
  #ifdef HAVE_OPENMP
  #pragma omp parallel for reduction(+ : overlap)
  #endif  
  for (int z = 0; z < maskOne->depth; ++z)
    for (int y = 0; y < maskOne->height; ++y)
      for (int x = 0; x < maskOne->width; ++x)
      {
        unsigned char voxOne = MRIvox(maskOne,x,y,z);
        unsigned char voxTwo = MRIvox(maskTwo,x,y,z);
        if ( voxOne!=lblBackground && voxTwo!=lblBackground )
        {
          // overlap
          ++overlap;
          MRIsetVoxVal(mri,x,y,z,0,voxOne);
        }
        else if ( voxOne == lblBackground )
        {
          MRIsetVoxVal(mri,x,y,z,0,voxTwo);
        }
        else if ( voxTwo == lblBackground )
        {
          MRIsetVoxVal(mri,x,y,z,0,voxOne);
        }
      } // next x,y,z

  std::cout << " hemi masks overlap voxels = " << overlap << std::endl;
  return mri;
}


MRI* FilterLabel(MRI* initialMask,
                 const unsigned char lbl)
{
  MRI* mri = MRIalloc( initialMask->width,
                       initialMask->height,
                       initialMask->depth,
                       MRI_UCHAR);

  #ifdef HAVE_OPENMP
  #pragma omp parallel for
  #endif   
  for (int z = 0; z< initialMask->depth; ++z)
    for (int y = 0; y < initialMask->height; ++y)
      for (int x = 0; x < initialMask->width; ++x)
      {
        if ( MRIvox(initialMask,x,y,z) == lbl )
        {
          MRIsetVoxVal(mri, x,y,z,0, (unsigned char)(1) );
        }
        else
        {
          MRIsetVoxVal(mri, x,y,z,0, (unsigned char)(0) );
        }
      } // next x,y,z

  return mri;
}


std::string InitVolmask(const IoParams& params, MRI*& mriTemplate, std::vector<MRI*>& surfDistVector, std::vector<MRIS*>& inSurfVector)
{
  // process input files
  // will also resolve the paths depending on the mode of the application
  // (namely if the subject option has been specified or not)
  std::string outputPath;

  try
  {
    MRIS *surfLeftWhite = NULL, *surfLeftPial = NULL, *surfRightPial = NULL, *surfRightWhite = NULL;
    outputPath = LoadInputFiles(params,
                                mriTemplate,
                                surfLeftWhite,
                                surfLeftPial,
                                surfRightWhite,
                                surfRightPial);

    inSurfVector.push_back(surfLeftWhite);
    inSurfVector.push_back(surfLeftPial);
    inSurfVector.push_back(surfRightWhite);
    inSurfVector.push_back(surfRightPial);
  }
  catch (std::exception& e)
  {
    std::cerr << " Exception caught while processing input files \n"
              << e.what() << std::endl;
    exit(1);
  }

  MRI *dLeftWhite = NULL, *dLeftPial = NULL, *dRightWhite = NULL, *dRightPial = NULL;
  if (params.DoLH) {
      // allocate distance
      dLeftWhite = MRIalloc( mriTemplate->width,
				  mriTemplate->height,
				  mriTemplate->depth,
				  MRI_FLOAT );
      MRIcopyHeader(mriTemplate, dLeftWhite);
  
      dLeftPial = MRIalloc( mriTemplate->width,
				 mriTemplate->height,
				 mriTemplate->depth,
				 MRI_FLOAT);
      MRIcopyHeader(mriTemplate,dLeftPial);
  }
  if (params.DoRH) {
      // allocate distance
      dRightWhite = MRIalloc( mriTemplate->width,
				   mriTemplate->height,
			 	   mriTemplate->depth,
				   MRI_FLOAT );
      MRIcopyHeader(mriTemplate, dRightWhite);

      dRightPial = MRIalloc( mriTemplate->width,
				  mriTemplate->height,
				  mriTemplate->depth,
				  MRI_FLOAT);
      MRIcopyHeader(mriTemplate, dRightPial);
  }

  surfDistVector.push_back(dLeftWhite);
  surfDistVector.push_back(dLeftPial);
  surfDistVector.push_back(dRightWhite);
  surfDistVector.push_back(dRightPial);

  return outputPath;
}


void Cleanup(std::vector<MRI*>& surfDistVector, std::vector<MRIS*>& inSurfVector)
{
  for (int n = 0; n < surfDistVector.size(); n++) {
    if (surfDistVector[n] != NULL)
      MRIfree(&surfDistVector[n]);
    // ??? MRISfree() fails, why ???
    //if (inSurfVector[n] != NULL)
    //  MRISfree(&inSurfVector[n]);
  }  
}


void SaveSurfaceDistances(std::vector<MRI*>& surfDistVector, std::string outputPath, std::string outRoot)
{
  for (int n = 0; n < surfDistVector.size(); n++) {
    MRI *dist = surfDistVector[n];
    if (dist != NULL) {
      const char *hemi = (n < 2) ? "lh" : "rh";
      const char *surfname = (n%2 == 0) ? "white" : "pial";

      char dist_mgz[256] = {'\0'};
      sprintf(dist_mgz, "%s/%s.d%s.%s.mgz", outputPath.c_str(), hemi, surfname, outRoot.c_str());
      MRIwrite(dist, dist_mgz);
      printf("saved %s %s signed distance field (%s)\n", hemi, surfname, dist_mgz);      
    }
  }
}


void SaveSurfaceDistances(MRI *surfDist, int threadid, const char *hemi, const char *surfname, std::string outputPath, std::string outRoot)
{
  char dist_mgz[256] = {'\0'};
  sprintf(dist_mgz, "%s/%s.d%s.%s.mgz", outputPath.c_str(), hemi, surfname, outRoot.c_str());
  MRIwrite(surfDist, dist_mgz);
  printf("[thread %d] saved %s %s signed distance field (%s)\n", threadid, hemi, surfname, dist_mgz);      
}


void SaveHemispheresRibbon(const char *hemi, MRI *mriTemplate, MRI *maskHemi, const unsigned char labelRibbon, std::string outputPath, std::string outRoot)
{
  int tid = 0;
#ifdef HAVE_OPENMP
  tid = omp_get_thread_num();
#endif

  // filter the mask of the hemi
  MRI* ribbon = FilterLabel(maskHemi, labelRibbon);
  MRIcopyHeader( mriTemplate, ribbon);
  ribbon->ct = CTABreadDefault();

  char ribbon_mgz[256] = {'\0'};
  sprintf(ribbon_mgz, "%s/%s.%s.mgz", outputPath.c_str(), hemi, outRoot.c_str());
  MRIwrite( ribbon, ribbon_mgz);   //const_cast<char*>( (outputPath / hemi + "." + outRoot + ".mgz").c_str()  ));
  printf("[thread %d] saved %s cortical ribbon mask (%s)\n", tid, hemi, ribbon_mgz);  

  // sanity-check: make sure location 0,0,0 is background (not brain)
  if( MRIgetVoxVal(ribbon,0,0,0,0) != 0 )    {
    fprintf(stderr, "[ERROR] [thread %d] %s ribbon has non-zero value at location 0,0,0\n", tid, hemi);
    exit(1);
  }
  MRIfree(&ribbon);
}
