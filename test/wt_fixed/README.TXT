* ./test/wt_fixed: Fixed WT embryo used for the "Methods" section (Figs. 1-6 of the paper).

* Input script overview:

       zf_img.inp      : Read input stack
       		         Includes stage 1 pixel-processing steps
       zf_boc.inp      : Read binary image from stage 1
       		         Includes stages 2 and 3 steps
       zf_mesh.inp     : Read data file from zf_boc.inp
       		         Writes .txt file for surface and normal drawing in VMD
		     	  	.stl file for surface rendering in MeshLab
       zf_img_vent.inp : Extract ventricle binary image (stage 1)
       zf_boc_vent.inp : Read binary from zf_img_vent.inp
       		       	 Generate ventricle BOC
       zf_measure.inp  : Collect measurements from ventricle BOC
       		       	 Outputs w_dist*.txt for each slice

* List of images to process (read by the above input scripts): wt_*.dat

* ./in: Input image files.
* ./out: Folder where outputs are stored.

----------------------------------------------------------------------
----------------------------------------------------------------------
[ Procedure ]

* Make symbolic link of the executable:
ln -s ../../bin/dodri

------------------------------
** Full embryo

  --> STAGE 1:
    ./dodri zf_img.inp
    	--> output: ./binary/wt_*.tif

  --> STAGES 2,3:
    ./dodri zf_boc.inp
        --> output: out/fwt.cor, out/fwt.psf (for viewing)
	    	    out/fwt.dat, out/bwt.dat (for reading)
        - note: uncomment "write" commands in zf_boc.inp to view progressive BOCs
	
  --> Generate files to visualize the mesh and normal vectors:
    ./dodri zf_mesh.inp
        --> output: out/fwt_quad.txt, out/fwt_quad_norm.txt
		    out/fwt.stl
		    
 --> TO VIEW
     BOC (Fig. 1G):		   vmd -e view_fwt.vmd
     With mesh [in vmd terminal]:  source out/fwt_quad.txt
 	   	  --> may need to click "=" on keyboard 
     	  	     	 	   source out/fwt_quad_norm.txt
     For MeshLab:  meshlab ./out/fwt.stl 				   

------------------------------
** Ventricle

  --> STAGE 1
    ./dodri zf_img_vent.inp
        --> output: binary/vent_mhb_*.tif   : MB, HB ventricles only
		    binary/all_*.tif        : MB, HB, FB merged with input
 	- note: uncomment "write" commands to view individual regions

  --> STAGES 2,3
    ./dodri zf_boc_vent.inp
        --> output: out/fvt.cor, out/fvt.psf (for viewing)
	      	    out/fvt.dat, out/bvt.dat (for reading)
		    
  --> Measure ventricle
    ./dodri zf_measure.inp
        --> output: out/fvt_orie.{cor,psf}  (oriented BOC)
	    	    w_dist*.txt             (measurements per slice)
    mv w_dist* distance/.
	- If no measurements are made for a slice, w_dist*.txt for that slice is not created 
	- Written output includes MHB position, position of LOCATEd
    	  ventricle widths, and measurements. Positions are printed to copy
    	  and paste in VMD for checking

  --> TO VIEW
      BOC:   vmd -e view_fvt_orie.vmd
      Width [in vmd terminal; *==slice number]:  source distance/w_dist*.txt   
      	    - To view all slices (like Fig.S4C):
	      bash dist.sh # output: distance/w_dist_all.txt
	      [in vmd] source distance/w_dist_all.txt
	    - To view MHB positions: In the printed output of
	    './dodri zf_measure.inp', load into the VMD console
   	      any of the "draw sphere ... radius 3" folowing "# MHB position:"
    	    - To view measured widths: In the printed output of
	    './dodri zf_measure.inp', load into the VMD console
	      any of the "draw line ... width 4" following "#  LOCATE ventricles, MHB"	       

	       

