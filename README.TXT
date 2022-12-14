** DODRI source code files for surface mesh building based on 3-dimensional image stacks

* Authors: Ana C. Chang-Gonzalez and Wonmuk Hwang

* Contact: Wonmuk Hwang, hwm@tamu.edu (see COPYRIGHT.TXT)

* How to cite:

Title: Building a three-dimensional model of early-stage zebrafish embryo brain

Authors: Ana C. Chang-Gonzalez, Holly C. Gibbs, Arne C. Lekven, Alvin T. Yeh, and Wonmuk Hwang

Journal: Biophysical Reports, Volume 1, Issue 1. 100003

Year: 2021


----------------------------------------------------------------------
----------------------------------------------------------------------

[Source code contents]

  ./inc/  : program header files 
  ./src/  : *.cpp source files, includes ./io/ folder with IO-related files
  ./obj/  : compiled object files  
  ./bin/  : dodri executable
  ./test/ : test files and input datasets 


[How to compile]

* Prerequisites:
  - C++ compiler: Although we compiled the code using the GNU g++, since the code uses only C++ STL, any C++ compiler should work.
  
  - GraphicsMagick++: http://www.graphicsmagick.org/Magick++/
    For reading/writing image files.

  - This code has been compiled and tested on the Ubuntu 20.04LTS operating system.

* To compile:
     make clean
     make


[To execute/test the code]

  ------------
  NOTE:
  
  0) A basic knowledge of the Linux operating system is desired since
  commands in the README.txt file in each folder are meant to be entered
  into a Linux terminal window (BASH environment). In Mac OS, a similar
  terminal can be used. In Microsoft Windows, a Linux emulation program
  can be used.
  
  1) Visualization was done mostly using the VMD software.
  
  ------------

Go to one of the sub-folders in ./test, and follow instructions in README.txt.
Each sub-folder is briefly explained below, in relation to our reference paper.

 ./test/wt_fixed: Fixed WT embryo used for the "Methods" section (Figs. 1-6 of the paper). (SUGGESTION: start with this folder) 

 ./zf_time_lapse : Live time-lapse embryos (Figs. 7, 8, 10)

 ./fetal_brain    : Human Fetal Brain (Fig. 9; 33 weeks)  

- For visualization of the output files, following programs can be used. They can be installed free of charge (see README.txt in the above sub-folders):

    	VMD (https://www.ks.uiuc.edu/Research/vmd/)
    	Chimera (https://www.cgl.ucsf.edu/chimera/) 
    	MeshLab (https://www.meshlab.net/)


[Image data]

- ./test/wt_fixed/ and ./test/zf_time_lapse/ (folders named `in'):
  Provided by Holly Gibbs, Arne Lekven and Alvin Yeh.

- For ./test/fetal_brain: From
  brain-development.org/brain-atlases/fetal-brain-atlases/


