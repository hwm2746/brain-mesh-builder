* ./test/zf_time_lapse/ace/1: Time-lapse ace embryo at 20.5 hpf
  			      (Figs 7, 8 of the paper)
			      
* Input script overview:

  zface1_merge.inp     : Merge ./in/auto and ./in/gfp images to enhance MHB
  		       	 (not needed for fixed embryo)
  zface1_img.inp       : Read merged images
        	       	 Includes stage 1 pixel-processing steps 
  zface1_boc.inp       : Read binary image from stage 1
  		       	 Includes stages 2, 3			 
  zface1_img_vent.inp  : Extract ventricle binary image (stage 1)
  zface1_boc_vent.inp  : Read binary from zface1_img_vent.inp
  		       	 Generate ventricle BOC
  zface1_measure.inp   : Collect measurements from ventricle BOC
       			 Outputs w_dist*.txt for each slice
  zface1_gfp.inp       : Assign weighed beads to GFP channel  

* List of images to process (read by the above input scripts): ace1*.dat

* ./in: Input image files.
* ./out: Folder where outputs are stored.

----------------------------------------------------------------------
----------------------------------------------------------------------
[ Procedure ]

* Make symbolic link of the executable:
ln -s ../../../../bin/dodri

------------------------------
** Enhance MHB by merging autofluorecence and GFP images
   (this is not done for wt_fixed)

./dodri zface1_merge.inp
    --> output: ./in/merge/imrg_*.tiff (for reading)

------------------------------
** Full embryo

  --> STAGE 1
    ./dodri zface1_img.inp
        --> output: ./binary/iaa_*.tif

  --> STAGES 2,3
    ./dodri zface1_boc.inp
    	--> output: ./out/faa.{cor,psf}

  --> TO VIEW
      BOC (Fig. 7D): vmd -e view_faa.vmd

------------------------------
** Ventricle

  --> STAGE 1
    ./dodri zface1_img_vent.inp
        --> output: ./binary/vent_mhb_*.tif
	    	    ./out/vent_mhb.mrc (for viewing with Chimera)

  --> STAGES 2,3
    ./dodri zface1_boc_vent.inp
        --> output: ./out/fav.{cor,psf}
	    	    ./out/fav.dat, ./out/bav.dat

  --> Measure ventricle
    ./dodri zface1_measure.inp
        --> output: ./out/fwv_orie.{cor,psf}   (oriented BOC)
	    	    w_dist*.txt 	 (measurements per slice)
    mv w_dist* distance/.
   	- If no measurements are made for a slice, w_dist*.txt for that slice is not created 
	- Written output includes MHB position, position of LOCATEd
    	  ventricle widths, and measurements. Positions are printed to copy
    	  and paste in VMD for checking

  --> TO VIEW 
      Electron density map    : chimera ./out/vent_mhb.mrc
      BOC (like Fig. 8C)      : vmd -e view_fav_orie.vmd  
      Width: [in vmd terminal; *==slice number]  source distance/w_dist*.txt   
      	     - To view all slices (like Fig.S4C):
	       bash dist.sh # output: distance/w_dist_all.txt
	       [in vmd] source distance/w_dist_all.txt
 	     - To view MHB positions: In the printed output of
 	     './dodri zface1_measure.inp', load into the VMD console
    	      any of the "draw sphere ... radius 3" folowing "# MHB position:"
     	     - To view measured widths: In the printed output of
 	     './dodri zface1_measure.inp', load into the VMD console
	      any of the "draw line ... width 4" following "#  LOCATE ventricles, MHB"	       
	       
	     - To view MHB positions: load any of the "draw sphere"
      	       positions under the output "# MHB position:"
	     - To view measured widths: load any of the "draw line"
      	       positions under the output "#  LOCATE ventricles, MHB"

------------------------------
** GFP

  --> Build GFP beads
    ./dodri zface1_gfp.inp
        --> output: bag.{cor,psf}

  --> TO VIEW
      Overlay with auto (like Fig. 10B, earlier time point):
      	    vmd -e view_faa_gfp.vmd
