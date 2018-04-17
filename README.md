# cam-raftint-mv-camera
acquisition code for capturing images from MatrixVision BlueFOX3 cameras, for guiding raft integration via edge detection
This code needs to be compiled against vendor distribution source tree residing in 
/opt/mvIMPACT_acquire on lsst-vw02.slac.stanford.edu. specifically this source code should be placed in 
a directory created by mvIMPACT_acquire/apps/mknewappl.sh <app_name> where <app_name> may be QuadCameraScope. 
Copy the sole source file here (QuadCameraScope.cpp) to replace the file with the same name within 
the directory created by the sh script. Modify Makefile.inc therein to contain cfitsio header and library locations.
then type "make x86_64" from within mvIMPACT_acquire to recompile this as an application.

commandline arguments are being added. not much help on this yet, here at least..
