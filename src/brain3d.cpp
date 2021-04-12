#include "dodri.h"

/************************************************************************/
//Constructors
/************************************************************************/

brain3d::brain3d(string tag0)
{ /* Constructor for brain3d class without pre-loaded objects.
     This lets individual functions call necessry tags from 
     img3d, bead3d, fnet3d, and surface objects. */
  
  tag=tag0;
  cout<<"  brain3d object= "<<tag0<<" initialized"<<endl;
}

/************************************************************************/
// Member functions
/************************************************************************/

void brain3d::analyzeZfVent(fnet3d *fn3d,double arr[],int nslice,
			    int dz,int *id_vec,string ax0)  
{ /* Measure ventricle values for zf. 
     Get: mesencephalic width, rhombencephalic width, MHB th
          respective positions (in VMD coors), calc MHB angle
	  
	   arr,nslice,ax0: for fnet3d::measureWidth3D
	   dz: number of slices to consider for locating ventricle positions in 3D 
	   set_mhb: set first position for mhb loc: [istk,pos]  */ 
  
  // axis_pos[istk][pt][cmx,cmy] -> cmz not needed b/c same as zcoord 
  // end_pos{0,1}[istk][pt][cmx,cmy] -> pt order should match b/en axis and end
  // range_loc[istk][2]: mhb opening range on axis_pos (stacks where mhb is not present; else =-1;
  // mhb_loc[istk]: mhb estimated location on axis_pos (stacks after mhb located; else ==-1)
  // wide_loc[istk][2]: loc of maximum width measured; iterator 0=before, 1=after mhb_loc
  // axis_width[istk][pt] -> thickness (width) at axis_pos
  // vent{0,1}_pos and mhb_pos -> position of found vent and mhb 
  
  int i,j,istk,istk1;
  int iloc,iloc1,iloc2,iloc3,mhb_check,flag,skip;
  double rdum[3],rdum1[3],rdum2[3],range_mdpt[2];
  double scale,min_width,tol; 
  double dist,dist1,dist2,max_dist=-1,max_dist2=-1,min_dist=-1,ax_dist_avg,xdist,ydist;
  int range_loc[fn3d->nStack][2], wide_loc[fn3d->nStack][2];
  int mhb_loc[fn3d->nStack];

  double vent0_pos[2][3], vent1_pos[2][3]; // 2 ventricles, each with two sides and 3 coors
  double mhb_pos[2][3]; 
  
  vector<double> axis_dist, copy_vec,*axis_width;
  vector<array<double,2>> *axis_pos,*end_pos0,*end_pos1;
  axis_width=new vector<double>[fn3d->nStack]; 
  axis_pos=new vector<array<double,2>>[fn3d->nStack];
  end_pos0=new vector<array<double,2>>[fn3d->nStack];
  end_pos1=new vector<array<double,2>>[fn3d->nStack];
  
  cout<<"  anal_zf_vent: fnet3d="<<fn3d->tag<<" dz="<<dz<<endl; 

  vector<fnet3d>::iterator it_fnet3d;  
  it_fnet3d=findFnet3d(&vfnet3d,fn3d->tag);
  it_fnet3d->measureWidth3D(9,arr,nslice,"",id_vec,axis_pos,axis_width,
			    end_pos0,end_pos1,ax0); //3 ref pts (N=3*3)

  // initialize
  for (istk=id_vec[0];istk<=id_vec[1];istk++) {
    range_loc[istk][0]=range_loc[istk][1]=-1;
    wide_loc[istk][0]=wide_loc[istk][1]=-1;    
    mhb_loc[istk]=-1;
  }
  
  /*** Locate MHB position ***/
  mhb_check=0;       
  for (istk=id_vec[0]+1;istk<=id_vec[1];istk++) {
    if ((int)axis_pos[istk].size()<=2) {mhb_loc[istk]=-2; continue; } // -2== nothing there 
    axis_dist.clear(); 
    for (i=0;i<(int)axis_pos[istk].size()-1;i++) {
      for (j=0;j<2;j++) {
	rdum[j]=axis_pos[istk][i][j];
	rdum1[j]=axis_pos[istk][i+1][j];
      }
      dist=getDistPos(rdum,rdum1,0); //return sqrt(dist)
      axis_dist.push_back(dist);
    }
    ax_dist_avg =0; //calc avg dist    
    ax_dist_avg = (accumulate(axis_dist.begin(),axis_dist.end(),0.))/(int)axis_dist.size();
    tol=1./ax_dist_avg; //tolerance for axis distance difference

    // for the largest dist, save max dist positions in range_loc if > avg+tol; else=-1
    iloc=-1;
    iloc=distance(axis_dist.begin(),max_element(axis_dist.begin(),axis_dist.end()));
    range_loc[istk][0]=iloc; range_loc[istk][1]=iloc+1;
    
    skip=0; 
    scale=sqrt((int)axis_pos[istk].size()-1);
    if (axis_dist[iloc]>scale*(ax_dist_avg+tol)) { // found region of large opening
      if (mhb_loc[istk-1]==-1 ) {
	iloc=range_loc[istk-1][0]; iloc1=range_loc[istk-1][1]; 
	for (j=0;j<2;j++) {
	  rdum1[j]=axis_pos[istk-1][iloc][j];
	  rdum2[j]=axis_pos[istk-1][iloc1][j];
	}
	dist=getDistPos(rdum1,rdum2,0);
	scale=(sqrt((int)axis_pos[istk-1].size()-1))/2.; 
	if (dist>scale*(ax_dist_avg+tol)) { // mhb not here, but start checking for mhb 
	  mhb_check=1; skip=1;
	}
      }
    } //     if (axis_dist[iloc]>scale*(ax_dist_avg+tol)) { // found region of large opening

    // check: before mhb checked, verify that mhb_loc does exist in this slice
    if (!skip && mhb_loc[istk-1]==-1) {
      iloc=range_loc[istk][0]; 
      copy_vec=axis_dist;
      copy_vec.erase(copy_vec.begin()+iloc); //remove max dist value
      rdum[0]=axis_dist[iloc]; 
      rdum1[0]=*max_element(copy_vec.begin(),copy_vec.end());
      if (rdum[0]-rdum1[0]>ax_dist_avg+tol) skip=1; // no mhb in this slice
    }
    
    if (mhb_check && !skip) {    // potential first mhb_loc
      if (mhb_loc[istk-1]==-1 && range_loc[istk-1][0]!=-1) {
	// mhb_loc== axis point nearest to range_loc[istk-1] mdpt
	iloc=range_loc[istk-1][0]; iloc1=range_loc[istk-1][1]; 
	for (j=0;j<2;j++)
	  range_mdpt[j]=(axis_pos[istk-1][iloc][j]+axis_pos[istk-1][iloc1][j])/2.;
	min_dist=1./TOL; 
	for (i=1;i<((int)axis_pos[istk].size())-1;i++) {
 	  for (j=0;j<2;j++) rdum[j]=axis_pos[istk][i][j];
	  dist=getDistPos(rdum,range_mdpt,0);
	  if (dist<min_dist) { //save min_dist, mhb_pos
	    min_dist=dist;
	    mhb_loc[istk]=i;
	  }
	} // 	for (i=1;i<((int)axis_pos[istk].size())-1;i++) {

	// check: if diff b/en istk range_loc < ax_dist_avg+tol,
	//  likely not first mhb_loc b/c mes/rho vent opening not present 
	if (mhb_loc[istk]>=0) { 
	  iloc=range_loc[istk][0]; iloc1=range_loc[istk][1]; 
	  for (j=0;j<2;j++) {
	    rdum1[j]=axis_pos[istk][iloc][j];
	    rdum2[j]=axis_pos[istk][iloc1][j];
	  }
	  dist=getDistPos(rdum1,rdum2,0);
	  scale=sqrt((int)axis_pos[istk].size()-1);
	  if (dist>scale*(ax_dist_avg+tol)) mhb_loc[istk]=-1; //not first mhb_loc
	}
      } //       if (mhb_loc[istk-1]==-1 && range_loc[istk-1][0]!=-1) {
      else if (mhb_loc[istk-1]>=0 && range_loc[istk-1][0]!=-1) { 
	iloc=range_loc[istk-1][0]; iloc1=range_loc[istk-1][1];
	for (j=0;j<2;j++) {
	  rdum1[j]=axis_pos[istk-1][iloc][j];
	  rdum2[j]=axis_pos[istk-1][iloc1][j];
	}
	min_width=1./TOL;
	for (i=0;i<(int)axis_pos[istk].size();i++) {
	  for (j=0;j<2;j++) rdum[j]=axis_pos[istk][i][j];
	  dist1=getDistPos(rdum,rdum1,0);
	  dist2=getDistPos(rdum,rdum2,0);
	  //find most narrow axis point near istk-1 range_loc
	  if (dist1<ax_dist_avg+tol || dist2<ax_dist_avg+tol) {
	    if (axis_width[istk][i]<min_width) {
	      min_width=axis_width[istk][i];
	      mhb_loc[istk]=i; 
	    }
	  }
	}
      } //       else if (mhb_loc[istk-1]>=0 && range_loc[istk-1][0]!=-1) { 
    } //     if (mhb_check && !skip) {    // potential first mhb_loc 
    
    /* If mhb_loc was set, verify location is valid. */
    
    // check #1: dist b/en prospective mhb and neighbors < ax_dist_avg	
    if (mhb_loc[istk]>=0) { 
      flag=0; 
      iloc=mhb_loc[istk]; 	  
      for (j=0;j<2;j++) {
	rdum[j]=axis_pos[istk][iloc][j];
	rdum1[j]=axis_pos[istk][iloc+1][j];
	rdum2[j]=axis_pos[istk][iloc-1][j];	    
      }
      scale=(sqrt((int)axis_pos[istk].size()-1))/2.; // dec tolerance range 
      dist=getDistPos(rdum,rdum1,0);
      if (dist>scale*(ax_dist_avg+tol)) flag=1; //mhb_loc[istk]=-1; //not mhb
      dist=getDistPos(rdum,rdum2,0);
      if (dist>scale*(ax_dist_avg+tol)) flag=1;// mhb_loc[istk]=-1; //not mhb

      // if condition met above, check dist from mhb_loc[istk-1] to verify 
      if (flag) {
	if (mhb_loc[istk-1]>=0) { 
	  iloc1=mhb_loc[istk-1];
	  for (j=0;j<2;j++) rdum1[j]=axis_pos[istk-1][iloc1][j];
	  dist=getDistPos(rdum,rdum1,0); //rdum from above
	  if (dist>scale*(ax_dist_avg+tol)) mhb_loc[istk]=-1; //not mhb}
	}
	else  mhb_loc[istk]=-1; 
      }
    }
    
    // mhb_loc verified; update range_loc[istk] 
    if (mhb_loc[istk]>=0) { 
      // take care of ends: if mhb_loc is at end, move +/-1 s.t. range_loc is valid
      if (mhb_loc[istk]==0)  mhb_loc[istk]=1;
      else if (mhb_loc[istk]==(int)axis_pos[istk].size()-1)
	mhb_loc[istk]=mhb_loc[istk]-1;
      
      range_loc[istk][0]=mhb_loc[istk]-1;
      range_loc[istk][1]=mhb_loc[istk]+1;
    }
    
    // mes/rho: save widest pos above/below mhb_loc (if exists)
    if (mhb_loc[istk]>=0) {
      if (mhb_loc[istk]-1>=0) 
	wide_loc[istk][0]=distance(axis_width[istk].begin(),
				   max_element(axis_width[istk].begin(),
					       axis_width[istk].begin()+(mhb_loc[istk]-1)));
      if (mhb_loc[istk]+1<=(int)axis_width[istk].size()-1) 
	wide_loc[istk][1]=distance(axis_width[istk].begin(),
				   max_element(axis_width[istk].begin()+(mhb_loc[istk]+1),
					       axis_width[istk].end()));
    }
  } //   for (istk=id_vec[0]+1;istk<=id_vec[1];istk++) {
  
  // mes/rho axis positions (wide_loc) where range_loc==-1 (only one vent present and no mhb_loc)
  for (istk=id_vec[0]+1;istk<id_vec[1];istk++) {
    if (wide_loc[istk][0]==-1 || wide_loc[istk][1]==-1) { // if either missing
      // skip if one already defined from a previous check 
      if (wide_loc[istk][0]!=-1 || wide_loc[istk][1]!=-1) continue;
      if (axis_width[istk].empty()) continue; // no width values at all in this slice
      iloc=distance(axis_width[istk].begin(),max_element(axis_width[istk].begin(),
  							 axis_width[istk].end()));
      for (j=0;j<2;j++) 
  	rdum[j]=axis_pos[istk][iloc][j];
      iloc3=-1;
      for (istk1=istk+1;istk1<=id_vec[1];istk1++) {  // get non -1 wide_loc instance
  	iloc1=wide_loc[istk1][0]; iloc2=wide_loc[istk1][1];
  	if (iloc1==-1 || iloc2==-1) continue;
  	for (j=0;j<2;j++) {
  	  rdum1[j]=axis_pos[istk1][iloc1][j];
  	  rdum2[j]=axis_pos[istk1][iloc2][j];
  	}
  	dist1=getDistPos(rdum,rdum1,0) ;
  	dist2=getDistPos(rdum,rdum2,0) ;
  	iloc3=(dist1<dist2)?(0):(1);  // save vent option (0 or 1) nearest to iloc (iloc1 or iloc2)
  	break; // only need to compare one vent instance    
      } //    for (istk1=istk+1;istk1<=id_vec[1];istk1++) {  // get non -1 wide_loc instance
      // save max dist (iloc) to correct vent (0 or 1)
      wide_loc[istk][iloc3]=iloc;       
    }
  } 
  
  //print
  cout<<"  # MHB position:"<<endl;
  for (istk=id_vec[0]+1;istk<=id_vec[1];istk++) {
    iloc=mhb_loc[istk];
    if (iloc>=0) { // -1==no mhb, -2==no data
      cout<<"    #istk:"<<istk<<"    ";       
      cout<<"  draw sphere {"<<axis_pos[istk][iloc][0]<<" "<<
	axis_pos[istk][iloc][1]<<" "<<fn3d->fnet_stack[istk].zCoord<<"} radius 3"<<endl;
    }
  }
  
  /**********************************/
  /*** Measure:mes_w, rho_w, mhb_th ***/
  //mhb_th: save max_width from mhb_loc widths

  cout<<" # LOCATE ventricles, MHB"<<endl; 
  
  // Locate position of 1st and 2nd widest pts 
  iloc1=-1; iloc2=-1; 
  //initiate max_dist and max_dist2
  for (i=0;i<2;i++) {
    for (istk=id_vec[0]+1;istk<=id_vec[1];istk++) { // from rho/mes 2D axis, find 2 thickest
      iloc=wide_loc[istk][i]; iloc1=wide_loc[istk+1][i];
      if (iloc==-1 || iloc1==-1) continue; 
      if (axis_width[istk][iloc]>axis_width[istk+1][iloc1]) {
  	max_dist=axis_width[istk][iloc];
  	max_dist2=axis_width[istk+1][iloc1];
      }
      else {
  	max_dist2=axis_width[istk][iloc];
  	max_dist=axis_width[istk+1][iloc1];
      }
      break; // just need to find initial max_dist and max_dist2 comparison 
    }
    
    for (istk=id_vec[0]+1;istk<=id_vec[1];istk++) { // from rho/mes 2D axis, find 2 thickest 
      iloc=wide_loc[istk][i];
      if (iloc==-1) continue;
      if (axis_width[istk][iloc]>max_dist) {
  	max_dist2=max_dist; // second=largest
  	iloc2=iloc1; 
  	max_dist=axis_width[istk][iloc]; // largest update
  	iloc1=iloc; istk1=istk;
      }
      else if (axis_width[istk][iloc]>max_dist2 && axis_width[istk][iloc]!=max_dist) {
  	max_dist2=axis_width[istk][iloc]; iloc2=iloc;
      }
    }
    
    // get plane where range_loc line is normal
    if (iloc1==-1) error_dodri("brain3d::analyzeZfVent: no reference line specified."); 
  } //   for (i=0;i<2;i++) {
  
  // {rho,mes,mhb}_w: locate rho/mes positions for 3D distance measure 
  //  criteria: maximize x-dist, minimize y-dist, within 0.5*slope from slope[i]  
  for (i=0;i<3;i++) { // for each section: wide_loc[istk][i]
    max_dist=TOL;
    for (istk=id_vec[0]+dz;istk<=id_vec[1]-dz;istk++) { // from rho/mes 2D axis, find 2 thickest
      if (i==0 || i==1) iloc=wide_loc[istk][i];
      if (i==2)  iloc=mhb_loc[istk];
      if (iloc<=0) continue;       // skip undefined slices
      for (istk1=istk-dz;istk1<=istk+dz;istk1++) {
	if (i==0 || i==1) iloc1=wide_loc[istk1][i];
	if (i==2)  iloc1=mhb_loc[istk1];
	if (iloc1<=0) continue;
	xdist=abs(end_pos0[istk][iloc][0]-end_pos1[istk1][iloc1][0]); // maximum x-dist
	ydist=abs(end_pos0[istk][iloc][1]-end_pos1[istk1][iloc1][1]); // minimum y-dist
	if (xdist>max_dist) { //save this point combo
	  for (j=0;j<2;j++) rdum[j]=end_pos0[istk][iloc][j];
	  rdum[2]=fn3d->fnet_stack[istk].zCoord;
	  for (j=0;j<2;j++) rdum1[j]=end_pos1[istk1][iloc1][j];
	  rdum1[2]=fn3d->fnet_stack[istk1].zCoord;
	  max_dist=xdist;
	}
	xdist=abs(end_pos1[istk][iloc][0]-end_pos0[istk1][iloc1][0]); // maximum x-dist - other end
	if (xdist>max_dist) { //save this point combo
	  for (j=0;j<2;j++) rdum[j]=end_pos0[istk][iloc][j];
	  rdum[2]=fn3d->fnet_stack[istk].zCoord;
	  for (j=0;j<2;j++) rdum1[j]=end_pos1[istk1][iloc1][j];
	  rdum1[2]=fn3d->fnet_stack[istk1].zCoord;
	  max_dist=xdist;	
	}
	if (ydist<min_dist) { 
	  for (j=0;j<2;j++) rdum2[j]=end_pos0[istk][iloc][j];
	  rdum2[2]=fn3d->fnet_stack[istk].zCoord;
	  min_dist=ydist;
	}
      }
    } //     for (istk=id_vec[0]+dz;istk<=id_vec[1]-dz;istk++) {
    cout<<"   draw line {"<<rdum[0]<<" "<<rdum[1]<<" "<<rdum[2]<<"} {"<<
      rdum1[0]<<" "<<rdum1[1]<<" "<<rdum1[2]<<"} width 4"<<endl;
    
    // save vent and mhb positions
    if (i==0)  {
      for (j=0;j<3;j++) vent0_pos[0][j]=rdum[j];
      for (j=0;j<3;j++) vent0_pos[1][j]=rdum1[j];      
    }
    if (i==1) {
      for (j=0;j<3;j++) vent1_pos[0][j]=rdum[j];
      for (j=0;j<3;j++) vent1_pos[1][j]=rdum1[j];      
    }
    if (i==2) {
      for (j=0;j<3;j++) mhb_pos[0][j]=rdum[j];
      for (j=0;j<3;j++) mhb_pos[1][j]=rdum1[j];      
    }
  } //   for (i=0;i<3;i++) { // for each section: wide_loc[istk][i]

  /** PRINT MEASUREMENTS **/ 
  cout<<" # PRINT MEASUREMENTS:"<< endl ;
  dist=getDistPosNdim(3,vent0_pos[0],vent0_pos[1],0); 
  cout<<"   vent_width="<<dist<<endl; 
  dist=getDistPosNdim(3,vent1_pos[0],vent1_pos[1],0); 
  cout<<"   vent_width="<<dist<<endl;
  dist=getDistPosNdim(3,mhb_pos[0],mhb_pos[1],0); 
  cout<<"   mhb_th="<<dist<<endl; 
  for (i=0;i<2;i++) {
    for (j=0;j<3;j++) {
      rdum[j]=vent0_pos[i][j]-mhb_pos[i][j];
      rdum1[j]=vent1_pos[i][j]-mhb_pos[i][j];    
    }
    cout<<"   mhb_bend:"<< getAngle(rdum,rdum1,0)*rad2degree<<endl; 
  }
    
  return; 
}

/************************************************************************/
void brain3d::writeSTL(string surf_tag,string fname)
{ /* Write bead-bond model to STL file. */

  FILE *stl;
  string ofname;
  int i,j,k,istk,kstk,bead,iStack,iQuad;
  unsigned udum;
  double triNorm[3],triVert[3][3],rdum1[3];
  
  // determine filename extension
  ofname = addFileNameExt(fname,".stl");
  cout<<"  Write brain3d to "<<ofname<<endl; 
  stl=fopen(ofname.c_str(),"w");

  fprintf(stl,"solid %s \n","zf");

  //FORMAT: triangles, vertex order follows RHR
  //facet normal n1 n2 n3
  // fprintf(stl,"outer loop");
  //   vertex x0 y0 z0
  //   vertex x1 y1 z1
  //   vertex x2 y2 z2
  // fprintf(stl,"endloop");
  //fprintf(stl,"endfacet");

  //build from surface to get triangles/normals
  vector<surface>::iterator it_surf;
  it_surf=findSurface(&vsurface,surf_tag); 

  for (istk=it_surf->fn3d_id_vec[0];istk<it_surf->fn3d_id_vec[1];istk++) {
    kstk=istk+1;
    for (udum=0;udum<(unsigned)it_surf->quad[istk].size();udum++) {

      /*** triangle 1 of quad ***/      
      iQuad=it_surf->quad[istk][udum].at(0)*2; // to match quad_normal indexing
      for (j=0;j<3;j++)
	triNorm[j]=it_surf->quad_normal.find(iQuad)->second.at(j);  // triangle normal
      for (i=1;i<4;++i) { // triangle vertices
	if (i!=3)  iStack=istk;
	else iStack=kstk;
	bead=it_surf->quad[istk][udum].at(i);
	triVert[i-1][0]=it_surf->fnet_stack[iStack].cm[bead][0];
	triVert[i-1][1]=it_surf->fnet_stack[iStack].cm[bead][1];
	triVert[i-1][2]=it_surf->fnet_stack[iStack].zCoord;
      }
      if (!verifyTriRHR(triVert,triNorm)) {
	for (k=0;k<3;k++) {
	  rdum1[k]=triVert[1][k]; //temp storage
	  triVert[1][k]=triVert[2][k];
	  triVert[2][k]=rdum1[k];
	}
      }
      fprintf(stl,"  facet normal %1.5E %1.5E %1.5E \n",triNorm[0],triNorm[1],triNorm[2]); 
      fprintf (stl,"      outer loop \n");             
      for (i=0;i<3;i++)
	fprintf (stl,"         vertex %1.5E %1.5E %1.5E \n",
		 triVert[i][0],triVert[i][1],triVert[i][2]); 
      fprintf(stl,"      endloop \n");
      fprintf(stl,"   endfacet \n");

      
      /*** triangle 2 of quad ***/
      iQuad=it_surf->quad[istk][udum].at(0)*2+1; // to match quad_normal indexing      
      for (j=0;j<3;j++)
	triNorm[j]=it_surf->quad_normal.find(iQuad)->second.at(j);  // triangle normal
      j=0; //counter
      for (i=1;i<5;++i) { //triangle vertices
	if (i==2) continue;
	if (i==1) iStack=istk;
	else iStack=kstk;
	bead=it_surf->quad[istk][udum].at(i);
	triVert[j][0]=it_surf->fnet_stack[iStack].cm[bead][0];
	triVert[j][1]=it_surf->fnet_stack[iStack].cm[bead][1];
	triVert[j][2]=it_surf->fnet_stack[iStack].zCoord;
	j++;	
      }
      if (!verifyTriRHR(triVert,triNorm)) {
	for (k=0;k<3;k++) {
	  rdum1[k]=triVert[1][k]; //temp storage
	  triVert[1][k]=triVert[2][k];
	  triVert[2][k]=rdum1[k];
	}
      }
      fprintf(stl,"  facet normal %1.5E %1.5E %1.5E \n",triNorm[0],triNorm[1],triNorm[2]);
      fprintf (stl,"      outer loop \n");      
      for (i=0;i<3;i++) 
	fprintf (stl,"         vertex %1.5E %1.5E %1.5E \n",
		 triVert[i][0],triVert[i][1],triVert[i][2]); 
      fprintf(stl,"      endloop \n");
      fprintf(stl,"   endfacet \n");
    }
  }
  
  fprintf(stl,"endsolid %s \n","zf");
  fclose(stl); 
  return;   
}

/************************************************************************/
// Non-member functions
/************************************************************************/

void checkDuplicateBrain3dTag(vector<brain3d> *vbrain3d, string tag)
{ /* In vfnet3d, find iterator for the element with tag. */
  
  vector<brain3d>::iterator it_brain3d;
  for (it_brain3d=(*vbrain3d).begin();it_brain3d!=(*vbrain3d).end();++it_brain3d) 
    if ( (*it_brain3d).tag==tag)
      error_dodri("Duplicate brain3d tag: " + tag);
  return;
}

/************************************************************************/
vector<brain3d>::iterator findBrain3d(vector<brain3d> *vbrain3d,string tag)
{ /* In vbrain3d, find iterator for the element with tag. */ 
  vector<brain3d>::iterator it_brain3d;
  for (it_brain3d=(*vbrain3d).begin();it_brain3d!=(*vbrain3d).end();++it_brain3d) {
    if ( (*it_brain3d).tag==tag) break;
  }
  return it_brain3d;
}
