* ./test/fetal_brain: Human Fetal Brain (Fig. 9, 33 weeks)

*  Dataset from: brain-development.org/brain-atlases/fetal-brain-atlases/
         Serag, A., Aljabar, P., Ball, G., Counsell, S. J.,
         Boardman, J. P., Rutherford, M. A., Edwards, A. D.,
         Hajnal, J. V., and Rueckert, D. (2012) Construction
         of a consistent high-definition spatio-temporal atlas
         of the developing brain using adaptive kernel regres-
         sion. Neuroimage, 59, 2255–2265.
	 
   --> Note: Images from download on brain-development.org are in *.nii.gz
       format. We use python3 med2image command to extract *.jpg images from
       cortex.nii.gz and ventricles.nii.gz. Respective folders in ./in/
       contain *.jpg cortex and ventricle images
       
* Input script overview:
       fetal_cortex.inp	    : Build cortex model 
       fetal_ventricle.inp  : Build ventricle model

* List of images to process (read by the above input scripts): fetal_*.dat

* ./in: Input image files.
* ./out: Folder where outputs are stored.

----------------------------------------------------------------------
----------------------------------------------------------------------
[ Procedure ]

* Make symbolic link of the executable:
ln -s ../../bin/dodri

./dodri fetal_cortex.inp
	--> output: ffc.cor, ffc.psf
./dodri fetal_ventricle.inp
	--> output: ffvt.cor, ffvt.psf

--> to view outputs (Fig. 9, 33 weeks, slanted view)
vmd -e view.vmd
 [Graphics -> Representations : toggle views on/off to see selected slices]  
