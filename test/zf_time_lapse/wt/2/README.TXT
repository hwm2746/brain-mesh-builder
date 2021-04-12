* ./test/zf_time_lapse/wt/2: Time-lapse WT embryo at 27.5 hpf
  			     (Figs. 7, 8, 10) 
* Input script overview:

       zft2_merge.inp      : Merge ./in/auto and ./in/gfp images to enhance MHB
       			     (not needed for fixed embryo)
       zft2_img.inp	   : Read merged images
       			     Includes stage 1 pixel-processing steps 
       zft2_boc.inp	   : Read binary image from stage 1
       			     Includes stages 2, 3 
       zft2_img_vent.inp   : Extract ventricle binary image (stage 1)
       zft2_boc_vent.inp   : Read binary from zft2_img_vent.inp
       			     Generate ventricle BOC
       zft2_measure.inp    : Collect measurements from ventricle BOC
       			     Outputs w_dist*.txt for each slice
       zft2_gfp.inp	   : Assign weighed beads to GFP channel
       
* List of images to process (read by the above input scripts): wt2*.dat

* ./in: Input image files.
* ./out: Folder where outputs are stored.

----------------------------------------------------------------------
----------------------------------------------------------------------
[ Procedure ]

* Make symbolic link of the executable:
ln -s ../../../../bin/dodri

------------------------------
**  Enhance MHB by merging autofluorecence and GFP images
   (this is not done for wt_fixed)

./dodri zft2_merge.inp
    --> output: ./in/merge/imrg*.tiff  (for reading to input scripts) 

------------------------------
** Full embryo

  --> STAGE 1
    ./dodri zft2_img.inp
    	--> output: ./binary/iwa_*.tif
    
  --> STAGES 2,3
    ./dodri zft2_boc.inp
    	--> output: ./out/fwa.{cor,psf}

  --> TO VIEW
      BOC (Fig. 7B): vmd -e view_fwa.vmd

------------------------------
** Ventricle 

  --> STAGE 1
    ./dodri zft2_img_vent.inp 
        --> output: ./binary/vent_mhb_*.tif
	    	    ./out/vent_mhb.mrc (for viewing with Chimera)

  --> STAGES 2,3
    ./dodri zft2_boc_vent.inp
        --> output: ./out/fwv.{cor,psf}
      	      	    ./out/fwv.dat, ./out/bwv.dat
	- note: uncomment first "write" to compare with initial BOC   

  --> Measure ventricle 
    ./dodri zft2_measure.inp
        --> output: ./out/fwv_orie.{cor,psf}   (oriented BOC)
	    	     w_dist*.txt 	 (measurements per slice)
    mv w_dist* distance/.
    	- If no measurements are made for a slice, w_dist*.txt for that slice is not created 
	- Written output includes MHB position, position of LOCATEd
    	  ventricle widths, and measurements. Positions are printed to copy
    	  and paste in VMD for checking

  --> TO VIEW
      Electron density map    : chimera ./out/vent_mhb.mrc
      BOC (like Fig. 8C)      : vmd -e view_fwv_orie.vmd  
      Width [in vmd terminal; *==slice number]:  source distance/w_dist*.txt   
      	     - To view all slices (like Fig.S4C):
	       bash dist.sh
	       [in vmd] source distance/w_dist_all.txt
 	     - To view MHB positions: In the printed output of
 	     './dodri zft2_measure.inp', load into the VMD console
    	      any of the "draw sphere ... radius 3" folowing "# MHB position:"
     	     - To view measured widths: In the printed output of
 	     './dodri zft2_measure.inp', load into the VMD console
	      any of the "draw line ... width 4" following "#  LOCATE ventricles, MHB"
	     - To view measured widths: load any of the "draw line"
      	       positions under the output "#  LOCATE ventricles, MHB"	      

------------------------------
** GFP

  --> Build GFP beads
    ./dodri zft2_gfp.inp
        --> output:  bwg.{cor,psf}

  --> TO VIEW:
      Overlay with auto (Fig. 10A):    vmd -e view_fwa_gfp.vmd  
