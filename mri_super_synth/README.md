This is a U-Net trained to make a set of useful predictions from any 3D brain image
(in vivo, ex vivo, single hemispheres, etc) using a common backbone/ It predicts:
- Segmentation:
- Registration to MNI atlas
- Joint super-resolution and synthesis of 1mm isotropic T1w, T2w, and FLAIR scans. 

The code relies on:
"A Modality-agnostic Multi-task Foundation Model for Human Brain Imaging"
Liu et al. (under revision)

##  Usage: 

The entry point / main script is mri_super_synth. There are two way of running the code:

A. For a single scan: just provide input file with --i, output directory with --o, and type of volume with --mode.

B. For a set of scans: you need to prepare a CSV file, where each row has 3 columns separated with commas:
- Column 1: input file
- Column 2: output directory
- Column 3: mode (must be invivo, exvivo, cerebrum, left-hemi, or right-hemi)

Please note that there is no leading/header row in the CSV file. The first row already corresponds to an input volume.
Tip: you can comment out a line by starting it with #

Important note: for 32 vs 64-bit reasons, inference is tiled on the GPU but not on the CPU, so results are expected to be 
slightly different on the 2 platforms. 
You can use --force_tiling option on the CPU to force tiling and get the same results as on the GPU

The command line options are:

  --i [IMAGE_OR_CSV_FILE]
                        Input image to segment - mode A - or CSV file with list of scans - mode B (required argument)

  --o [OUTPUT_DIRECTORY]
                        Directory where outputs will be written (ignored in mode B)

  --mode [MODE]
                        Type of input. Must be invivo, exvivo, cerebrum, left-hemi, or right-hemi (ignored in mode B)

  --threads [THREADS]     
                        Number of cores to be used. You can use -1 to use all available cores. Default is -1 (optional)

  --device [DEV]     
                        Device used for computations (cpu or cuda). The default is to use cuda if a GPU is available (optional)
                        
  --force_tiling     
                        Use this flag to force tiling on CPU and get the same results as on GPU, as explained above (optional)





## Prerequisites:

The first time you run the method, it will prompt you to download the machine learning model files, which are not distributed with the code.
