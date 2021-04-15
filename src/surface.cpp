#include "dodri.h"

/************************************************************************/
// Constructors
/************************************************************************/

surface::surface(fnet3d *fn3d, int *id_vec, string tag0) 
{ /* Constructor using existing fnet3d. 
     id_vec: id range for loaded images. */ 
  
  int i,ndigit;  
  tag=tag0;
  fnet3d_tag=fn3d->tag;
  
  //variables inherited from fnet3d
  for (i=0;i<2;i++) fn3d_id_vec[i]=id_vec[i]; 
  nStack=fn3d->nStack;
  fnet_stack=fn3d->fnet_stack;
  nQuad=0; 
  
  //Initialize 
  bead_quad=new multimap<int,int>[nStack];
  quad = new vector<array<int,5>>[nStack];
  cout<<"  "<<tag<<": "<<nStack<<" planes in the fnet3d stack."<<endl;
  
  ndigit=nDigit(nStack);
  cout<<"  "<<setw(ndigit)<<std::right<<
    tag<<": "<<nStack<<" planes in the surface stack."<<endl;
}

/************************************************************************/
// Member functions
/***********************************************************************/

void surface::addOneQuad(int istk, array<int,5> quadbeads)
{ /* Assign the quad defined by quadbeads. */

  int j,kstk=istk+1;
  int iQuad=quadbeads.at(0);
  if (checkDuplicateQuad(istk,quadbeads))
    error_dodri("  Duplicate quad being assigned");

  quad[istk].push_back(quadbeads);
  for (j=1;j<3;j++) 
    bead_quad[istk].insert( pair<int,int>(quadbeads.at(j),iQuad));
  for (j=3;j<5;j++)
    bead_quad[kstk].insert( pair<int,int>(quadbeads.at(j),iQuad));

  nQuad++;
}

/************************************************************************/
void surface::calcParam(int flag)
{ /* Calculates various parameters for quad.
     0 = calculate centroid (2 centroids for quad)
     1 = calculate unit normal vector   */
  
  int i,j,istk,kdum,ifig,ibead,bead,dum,temp,iQuad;
  array<double,3> normarr,bcmarr;
  double **coor = new double*[3];
  double *bcm=new double[3];
  double *bnorm=new double[3];
  for (i=0;i<3;i++) coor[i]=new double[3];
  
  for (istk=0;istk<nStack-1;istk++) {    //quad parameters
    for (ifig=0;ifig<(int)quad[istk].size();ifig++) {
      iQuad=quad[istk][ifig].at(0);
      for (dum=0;dum<2;dum++) {
	for (temp=1;temp<4;temp++) { 
	  if (dum==0) {ibead=temp;}
	  else if (dum==1 && temp==1) {ibead=1;}
	  else if (dum==1 && temp!=1) {ibead=temp+dum;}
	  if (ibead<3) kdum=istk;
	  else kdum=istk+1;
	  bead=quad[istk][ifig].at(ibead); 
	  coor[0][temp-1]=fnet_stack[kdum].cm[bead][0]; 
	  coor[1][temp-1]=fnet_stack[kdum].cm[bead][1];
	  coor[2][temp-1]=fnet_stack[kdum].zCoord;	 
	}
	if (flag==0) {
	  bcm[0]=bcm[1]=bcm[2]=0.;
	  getCentroid(coor,bcm,3,3);
	  for (j=0;j<3;j++) bcmarr[j]=bcm[j];
	  if (dum==0) quad_cm.insert(pair<int,array<double,3>>(iQuad*2,bcmarr));
	  else if (dum==1) quad_cm.insert(pair<int,array<double,3>>(iQuad*2+1,bcmarr));
	}
	if (flag==1) { //execute surface::quadNormOrie() to orient after
	  bnorm[0]=bnorm[1]=bnorm[2]=0.;
	  getNormal(coor,bnorm);
	  for (j=0;j<3;j++) normarr[j]=bnorm[j];
	  if (dum==0) quad_normal.insert(pair<int,array<double,3>>(iQuad*2,normarr));
	  else if (dum==1) quad_normal.insert(pair<int,array<double,3>>(iQuad*2+1,normarr));
	}
      }
    }
  } //   for (istk=0;istk<nStack-1;istk++) {

  if (flag==1 && nQuad!=0) quadNormOrie(); // orient normals out of surface
  
  //clear  memory
  for (i=0;i<3;i++) {delete[] coor[i];}
  delete[] coor; delete[] bcm; delete[] bnorm;
}

/************************************************************************/
int surface::checkDuplicateQuad(int istk,array<int,5> quadbeads)
{ /* Check if quad defined by quadbeads was already assigned in 
     registerQuad. */

  int bead1=quadbeads.at(1), bead2=quadbeads.at(2);
  int bead3=quadbeads.at(3), bead4=quadbeads.at(4), flag=0;
  unsigned udum;
  for (udum=0;udum<quad[istk].size();udum++) {
    if (quad[istk][udum].at(1)==bead1 && quad[istk][udum].at(2)==bead2) {
      if (quad[istk][udum].at(3)==bead3 && quad[istk][udum].at(4)==bead4)
	++flag;
      if (quad[istk][udum].at(3)==bead4 && quad[istk][udum].at(4)==bead3)
	++flag;
    }
  }
  return flag;
}

/***********************************************************************/
void surface::drawSurf(string fname,string color,string outMode,int *id_vec)
{ /* Output .txt file to visualize surface figures. 
     Color is surface triangle draw color, outMode is image or verbose,
     id_vec is to select stack range.   */
  
  int istk,kstk,i,bead,iStack,j,iQuad;
  unsigned udum;
  FILE *ff,*ff1;
  string ofname,ofname1,fname1;
  double x1, y1, z1;
  double pts_temp[3], normal_temp[3]; 

  vector<fnet3d>::iterator it_fnet3d;
  it_fnet3d=findFnet3d(&vfnet3d,this->fnet3d_tag); 

  // determine istk range if not set
  if (id_vec[0]==-1)  id_vec[0]=fn3d_id_vec[0]; 
  if( id_vec[1]==-1)  id_vec[1]=fn3d_id_vec[1];
  
  // for quadrilaterals
  ofname = addFileNameExt(fname,".txt");
  cout<<"    Write surface for fnet3d quadrilaterals to "<<ofname<<endl;
  cout<<"          color= "<<color<<" outMode="<<outMode<<endl;
  ff=fopen(ofname.c_str(),"w");
  
  fname1= fname+"_norm"; 
  ofname1 = addFileNameExt(fname1,".txt");
  ff1=fopen(ofname1.c_str(),"w"); 
  fprintf(ff1,"source vmd_draw_arrow.tcl\n");
  
  fprintf(ff,"mol load graphics %s \n",fname.c_str()) ;
  fprintf(ff,"draw material BrushedMetal\n") ;  
  fprintf(ff,"draw color %s\n",color.c_str()) ;
  
  for (istk=id_vec[0];istk<id_vec[1];istk++) {
    kstk=istk+1;
    
    for (udum=0;udum<quad[istk].size();udum++) { // draw 2 triangles to make a quad
      fprintf (ff,"draw triangle "); 
      for (i=1;i<4;++i) {
	if (i!=3)  iStack=istk;
	else iStack=kstk;
	bead=quad[istk][udum].at(i);
	if (outMode=="verbose") {
	  x1=fnet_stack[iStack].cm[bead][0];
	  y1=fnet_stack[iStack].cm[bead][1];
	  z1=fnet_stack[iStack].zCoord;
	  fprintf (ff,"{%9.4f %9.4f %9.4f \n} ",x1,y1,z1);
	}
	else { //outMode=="image"
	  y1=((double)(fnet_stack[iStack].rows)-fnet_stack[iStack].cm[bead][0]);
	  x1=fnet_stack[iStack].cm[bead][1];z1=fnet_stack[iStack].zCoord;
	  fprintf (ff,"{%9.4f %9.4f %9.4f} ",x1,y1,z1);
	}
      }
      fprintf (ff,"\n");
      iQuad=quad[istk][udum].at(0)*2; // match quad_normal indexing
      for (j=0;j<3;j++) { // normals
	pts_temp[j]=quad_cm.find(iQuad)->second.at(j);
	normal_temp[j]=quad_normal.find(iQuad)->second.at(j);
      }
      drawVMD(ff1,pts_temp,normal_temp,3.);
      fprintf (ff,"\n");      
      
      fprintf (ff,"draw triangle ");
      for (i=1;i<5;++i) {
	if (i==2) continue;
	if (i==1) iStack=istk;
	else iStack=kstk;
	bead=quad[istk][udum].at(i);
	if (outMode=="verbose") {
	  x1=fnet_stack[iStack].cm[bead][0];
	  y1=fnet_stack[iStack].cm[bead][1];
	  z1=fnet_stack[iStack].zCoord;
	  fprintf (ff,"{%9.4f %9.4f %9.4f \n} ",x1,y1,z1);
	}
	else { //outMode=="image"
	  y1=((double)(fnet_stack[iStack].rows)-fnet_stack[iStack].cm[bead][0]);
	  x1=fnet_stack[iStack].cm[bead][1];z1=fnet_stack[iStack].zCoord;
	  fprintf (ff,"{%9.4f %9.4f %9.4f} ",x1,y1,z1);
	}
      }
      fprintf (ff,"\n");
      iQuad=quad[istk][udum].at(0)*2+1; // match quad_normal indexing
      for (j=0;j<3;j++) { // normals
	pts_temp[j]=quad_cm.find(iQuad)->second.at(j);
	normal_temp[j]=quad_normal.find(iQuad)->second.at(j);
      }
      drawVMD(ff1,pts_temp,normal_temp,3.);      
      fprintf (ff,"\n");
    }
  }
  fclose(ff);  fclose(ff1);
  //for triangles removed
}

/************************************************************************/
void surface::drawVMD(FILE *ff, double d0[3], double d1[3], double a)
{ /* Draw arrow from point d0 to d1. a==magnitude.  */

  int i;
  double npt[3];
  for (i=0;i<3;i++) {
    npt[i]=(a*d1[i]+d0[i]);
  }
  fprintf(ff,"draw arrow {%9.4f %9.4f %9.4f} {%9.4f %9.4f %9.4f}\n",
	  d0[0],d0[1],d0[2],npt[0],npt[1],npt[2]);
}

/***********************************************************************/
void surface::quadNormOrie()
{ /* Orient the normal vectors of each quad to point in same dir as
     2D pixel intensity gradient. */
  
  int istk,stk,k,j,iQuad,ibd;
  double phi,rdum,avgGrad[3],triNorm[3];
  unsigned udum; 

  /* GET 2D pixel intensity gradient for all beads */
  // positive: gradient "points" from high to low intensity
  vector<fnet3d>::iterator it_fnet3d;
  vector<bead3d>::iterator it_bead3d; 
  it_fnet3d=findFnet3d(&vfnet3d,this->fnet3d_tag); 
  it_bead3d=findBead3d(&vbead3d, it_fnet3d->bead3d_tag);  //get bead tag 
  double rCut=1.5, rWindow;
  multimap<int,array<double,2>> *bd_grad;
  bd_grad=new multimap<int,array<double,2>>[nStack]; 
  for (istk=fn3d_id_vec[0];istk<=fn3d_id_vec[1];istk++) {
    rWindow=(double)it_bead3d->bead_stack[istk].beadDiam;
    if (it_bead3d->bead_stack[istk].nBead==0) continue;     
    it_bead3d->bead_stack[istk].gradient(rCut,rWindow,bd_grad[istk]);
  }
  
  for (istk=fn3d_id_vec[0];istk<=fn3d_id_vec[1];istk++) {
    if (it_bead3d->bead_stack[istk].nBead==0) continue; 
    for (udum=0;udum<(unsigned)quad[istk].size();udum++) {
      avgGrad[0]=avgGrad[1]=avgGrad[2]=0; // initialize ; z-dim stays=0    
      for (k=0;k<4;k++) {
	ibd=quad[istk][udum].at(k+1);
	if (k<2) stk=istk; //beads 0,1 are in istk
	else stk=istk+1;   //beads 1,2 are in istk+1
	for (j=0;j<2;j++)  avgGrad[j]+=bd_grad[stk].find(ibd)->second.at(j);
      }
      avgGrad[0]/=3.;  avgGrad[1]/=3.;
      rdum=sqrt(dotprod(avgGrad,avgGrad,2));
      for (j=0;j<2;j++) avgGrad[j]/=rdum;
      
      for (k=0;k<2;k++) {   // two sets of triangle normals per quad
	iQuad=quad[istk][udum].at(0)*2+k; // to match quad_normal indexing
	for (j=0;j<3;j++) triNorm[j]=quad_normal.find(iQuad)->second.at(j); 
	phi=getAngle(avgGrad,triNorm,3); 
	if (phi>PI/2.) { //flip normal orientation
	  for (j=0;j<3;j++)
	    quad_normal.find(iQuad)->second.at(j) = -1.*triNorm[j];
	  //else, keep the same orientation 
	}
      }
    } //     for (udum=0;udum<(int)quad[istk].size();udum++) {
  } //   for (istk=fn3d_id_vec[0];istk<=fn3d_id_vec[1];istk++) {
  
  return;
}

/***********************************************************************/
void surface::registerHexagon()
{ /* Register hexagonal surface : bead-bead-bead-\zbond-\bead-bead-bead
     Located between consecutive slices.

     Code will locate hexagonal gap, add zbond between middle beads,
     and save to bead_quad and quad structures. 
  */

  int i,istk,kstk,iBead,iQuad,nQuadold;
  int jBead0,jBead1,jBeadZ,jBeadZ0,jBeadZ1;
  vector<int> sideBeads; 
  array<int,5> quadbeads;
  array<int,7> hexgaparr;
  vector<array<int,7>> hexgap; 
  pair<multimap<int, int>::iterator, multimap<int, int>::iterator> it_range,
    jt_range,kt_range;
  multimap<int, int>::iterator it,jt,kt;
  
  vector<fnet3d>::iterator it_fnet3d;
  it_fnet3d=findFnet3d(&vfnet3d,this->fnet3d_tag);

  iQuad=nQuad;  nQuadold=nQuad; 
  hexgap.clear();
  
  for (istk=fn3d_id_vec[0];istk<fn3d_id_vec[1]-1;istk++) {
    
    kstk=istk+1; 
    for (iBead=0;iBead<fnet_stack[istk].nBead;iBead++) {
      sideBeads.clear(); 
      
      //search for beads un-z-bonded to next stack
      if (it_fnet3d->bond_a[istk].find(iBead)!=it_fnet3d->bond_a[istk].end())
	continue;  //if bead has a z bond, continue
      
      it_range=fnet_stack[istk].bond.equal_range(iBead); //bonded to iBead
      for (it=it_range.first;it!=it_range.second;it++) 
	sideBeads.push_back(it->second);
      if (sideBeads.size()!=2) continue;       //if more than two, skip for hex shape

      jBead0=sideBeads[0];      jBead1=sideBeads[1];

      //skip if jBead0 or jBead1 do not have z bonds
      if (it_fnet3d->bond_a[istk].find(jBead0)==it_fnet3d->bond_a[istk].end()) continue;
      if (it_fnet3d->bond_a[istk].find(jBead1)==it_fnet3d->bond_a[istk].end()) continue;
      
      jBeadZ0=it_fnet3d->bond_a[istk].find(jBead0)->second;  //assuming one zbond/bead
      jBeadZ1=it_fnet3d->bond_a[istk].find(jBead1)->second;  //assuming one zbond/bead      

      //if jBeadZ0 and jBeadZ1 are bonded to the same bead, match
      jt_range=fnet_stack[kstk].bond.equal_range(jBeadZ0);      
      for (jt=jt_range.first;jt!=jt_range.second;jt++) {
	kt_range=fnet_stack[kstk].bond.equal_range(jBeadZ1);
	for (kt=kt_range.first;kt!=kt_range.second;kt++) {
	  if (jt->second==kt->second) {
	    jBeadZ=jt->second; 
	    //ADD QUAD AND UPDATE Z BOND OUTSIDE OF LOOP
	    //quad order: [jbead0 ibead jbeadz jbeadz0] and [iBead jbead1 jbeadz1 jbeadz]
	    hexgaparr.at(0)=istk;
	    hexgaparr.at(1)=iBead; hexgaparr.at(2)=jBeadZ;
	    hexgaparr.at(3)=jBead0; hexgaparr.at(4)=jBeadZ0;
	    hexgaparr.at(5)=jBead1; hexgaparr.at(6)=jBeadZ1;
	    hexgap.push_back(hexgaparr); 
	    break;
	  }
	} // 	for (kt=kt_range.first;kt!=kt_range.second;kt++) {
      } //       for (jt=jt_range.first;jt!=jt_range.second;jt++) {
    } //     for (iBead=0;iBead<fnet_stack[istk].nBead;iBead++) {
  } // for (istk=fn3d_id_vec[0];istk<fn3d_id_vec[1];istk++) {
  
  //add quad and z bond
  for (i=0;i<(int)hexgap.size();i++) {
    istk=hexgap[i][0];
    iBead=hexgap[i][1]; jBeadZ=hexgap[i][2];
    jBead0=hexgap[i][3]; jBeadZ0=hexgap[i][4];
    jBead1=hexgap[i][5]; jBeadZ1=hexgap[i][6]; 
    
    //quad order: [jbead0 ibead jbeadz jbeadz0] 
    quadbeads.at(0)=iQuad; 
    quadbeads.at(1)=jBead0; quadbeads.at(2)=iBead;
    quadbeads.at(3)=jBeadZ;quadbeads.at(4)=jBeadZ0;
    if (!checkDuplicateQuad(istk,quadbeads)) {
      addOneQuad(istk,quadbeads); iQuad++; 
    }
    
    //and [iBead jbead1 jbeadz1 jbeadz]
    quadbeads.at(0)=iQuad; 
    quadbeads.at(1)=iBead; quadbeads.at(2)=jBead1;
    quadbeads.at(3)=jBeadZ1;quadbeads.at(4)=jBeadZ;
    if (!checkDuplicateQuad(istk,quadbeads)) {
      addOneQuad(istk,quadbeads); iQuad++; 
    }
    
    //add zbond between iBead and jBeadZ
    it_fnet3d->bond_a[istk].insert(pair<int,int>(iBead,jBeadZ));
    it_fnet3d->bond_b[istk+1].insert(pair<int,int>(iBead,jBeadZ));
  } //   for (i=0;i<(int)hexgap.size();i++) {
  
  cout<<"  register_hexagon: "<<nQuad-nQuadold<<" quadrilaterals added"<<endl;
  nQuad=iQuad;  
  return; 
}

/************************************************************************/
void surface::registerMissingSurf()
{ /* Locate beads w/o triangle surf to neighbor slice. 
     Add triangles if necessary. */ 
  
  int i,istk,jstk,ibd0,ibd1,jbd,jbd0,jbd1,iQuad,jFil,nQuadold=nQuad;
  array<int,5> quadbeads;
  vector<int> zEnds,jSeq; 
  vector<fnet3d>::iterator it_fnet3d;
  it_fnet3d=findFnet3d(&vfnet3d,this->fnet3d_tag);
  
  iQuad=nQuad;   
  for (istk=fn3d_id_vec[0];istk<fn3d_id_vec[1];istk++) {
    jstk=istk+1;
    for (jbd=0;jbd<fnet_stack[jstk].nBead;jbd++) {
      if (it_fnet3d->bond_b[jstk].find(jbd)==it_fnet3d->bond_b[jstk].end()) {
	zEnds.clear();	jSeq.clear();
	jFil=fnet_stack[jstk].bead_fil.find(jbd)->second;
	
	it_fnet3d->findNextZ(jstk,jbd,jFil,istk,zEnds);
	if ((int)zEnds.size()!=2) continue ; // to save ENCLOSED region 
	
	// verify and save istk z-bond ends  (use same beads for this sequence)
	if (it_fnet3d->bond_b[jstk].count(zEnds[0])!=1 ||
	    it_fnet3d->bond_b[jstk].count(zEnds[1])!=1) continue;  // skip issues
	ibd0=it_fnet3d->bond_b[jstk].find(zEnds[0])->second; // assume ibd0--ibd1 bond 
	ibd1=it_fnet3d->bond_b[jstk].find(zEnds[1])->second;

	if (it_fnet3d->bond_a[istk].find(ibd0)->second != zEnds[0] ||
	    it_fnet3d->bond_a[istk].find(ibd1)->second != zEnds[1]) continue;  

	if (fnet_stack[istk].bead_fil.find(ibd0)->second !=
	    fnet_stack[istk].bead_fil.find(ibd1)->second) continue; 
	
	it_fnet3d->trackFilSeq(jstk,zEnds[0],zEnds[1],jSeq);

	for (i=0;i<(int)jSeq.size()-1;i++) { 
	  jbd0=jSeq[i]; jbd1=jSeq[i+1]; 
	  quadbeads.at(0)=iQuad;
	  quadbeads.at(1)=ibd0; quadbeads.at(2)=ibd1;
	  quadbeads.at(3)=jbd0; quadbeads.at(4)=jbd1;	
	  if (!checkDuplicateQuad(istk,quadbeads)) {
	    addOneQuad(istk,quadbeads);
	    iQuad++; 
	  }
	} // 	for (i=0;i<(int)jSeq.size()-1;i++) {
	
	for (i=0;i<(int)jSeq.size()-1;i++) { // other side check
	  jbd0=jSeq[i]; jbd1=jSeq[i+1]; 
	  quadbeads.at(0)=iQuad;
	  quadbeads.at(1)=ibd1; quadbeads.at(2)=ibd0;
	  quadbeads.at(3)=jbd0; quadbeads.at(4)=jbd1;	
	  if (!checkDuplicateQuad(istk,quadbeads)) {
	    addOneQuad(istk,quadbeads);
	    iQuad++; 
	  }
	} // 	for (i=0;i<(int)jSeq.size()-1;i++) {
      }
    } //      for (jbd=0;jbd<fnet_stack[jstk].nBead;jbd++) {
  } //    for (istk=fn3d_id_vec[0];istk<fn3d_id_vec[1];istk++) {

  // other side check
  iQuad=nQuad;   
  for (istk=fn3d_id_vec[0]+1;istk<=fn3d_id_vec[1];istk++) {
    jstk=istk-1;
    for (jbd=0;jbd<fnet_stack[jstk].nBead;jbd++) {
      if (it_fnet3d->bond_a[jstk].find(jbd)==it_fnet3d->bond_a[jstk].end()) {
	zEnds.clear();	jSeq.clear();
	jFil=fnet_stack[jstk].bead_fil.find(jbd)->second;
	it_fnet3d->findNextZ(jstk,jbd,jFil,istk,zEnds);
	if ((int)zEnds.size()!=2) continue ; // to save ENCLOSED region
	
	// verify and save istk z-bond ends  (use same beads for this sequence)
	if (it_fnet3d->bond_a[jstk].count(zEnds[0])!=1 ||
	    it_fnet3d->bond_a[jstk].count(zEnds[1])!=1) continue;  // skip issues
	ibd0=it_fnet3d->bond_a[jstk].find(zEnds[0])->second; 
	ibd1=it_fnet3d->bond_a[jstk].find(zEnds[1])->second;
	if (!fnet_stack[istk].isBonded(ibd0,ibd1)) continue; 				       
	if (it_fnet3d->bond_b[istk].find(ibd0)->second != zEnds[0] ||
	    it_fnet3d->bond_b[istk].find(ibd1)->second != zEnds[1]) continue;  
	if (fnet_stack[istk].bead_fil.find(ibd0)->second !=
	    fnet_stack[istk].bead_fil.find(ibd1)->second) 	
	it_fnet3d->trackFilSeq(jstk,zEnds[0],zEnds[1],jSeq);

	for (i=0;i<(int)jSeq.size()-1;i++) { 
	  jbd0=jSeq[i]; jbd1=jSeq[i+1]; 
	  quadbeads.at(0)=iQuad;
	  quadbeads.at(1)=jbd0; quadbeads.at(2)=jbd1;
	  quadbeads.at(3)=ibd0; quadbeads.at(4)=ibd1;	
	  if (!checkDuplicateQuad(jstk,quadbeads)) {
	    addOneQuad(jstk,quadbeads);
	    iQuad++; 
	  }
	} // 	for (i=0;i<(int)jSeq.size()-1;i++) {

	for (i=0;i<(int)jSeq.size()-1;i++) { 
	  jbd0=jSeq[i]; jbd1=jSeq[i+1]; 
	  quadbeads.at(0)=iQuad;
	  quadbeads.at(1)=jbd1; quadbeads.at(2)=jbd0;
	  quadbeads.at(3)=ibd0; quadbeads.at(4)=ibd1;	
	  if (!checkDuplicateQuad(jstk,quadbeads)) {
	    addOneQuad(jstk,quadbeads);
	    iQuad++; 
	  }
	} // 	for (i=0;i<(int)jSeq.size()-1;i++) {
      } //       if (it_fnet3d->bond_a[jstk].find(jbd)==it_fnet3d->bond_a[jstk].end()) {
    } //     for (jbd=0;jbd<fnet_stack[jstk].nBead;jbd++) {
  } //   for (istk=fn3d_id_vec[0]+1;istk<=fn3d_id_vec[1];istk++) {
  
  cout<<"  register_missing: "<<nQuad-nQuadold<<
    " quadrilaterals added"<<endl;   
  return; 
}

/************************************************************************/
void surface::registerPent()
{ /* Register regions that are two quads connected by a filament. Register
     missing filament and create quad.  Geometry: quad1--fil--quad2. 
     Locates pentagon-like, single-slice, gaps.   */
  
  int i,iFil,istk,kstk,nVec,jBead,kBead;
  int iQuad=nQuad,bead1,bead2,bead3,bead4,bead_mid,nQuadold=nQuad;
  unsigned udum;
  int *filSeq=new int[maxVec];
  array<int,5> quadbeads;
  multimap<int,int>::iterator jt,kt,zt,zt0,qt;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> jt_range,
    kt_range,zt_range,zt0_range,qt_range;
  vector<fnet3d>::iterator it_fnet3d;
  it_fnet3d=findFnet3d(&vfnet3d,this->fnet3d_tag); 

  //for bead insert
  vector<bead>::iterator it_bead;
  vector<fnet>::iterator it_fnet;  
  vector<bead3d>::iterator it_bead3d;
  it_bead3d=findBead3d(&vbead3d,it_fnet3d->bead3d_tag);
    
  for (istk=fn3d_id_vec[0];istk<fn3d_id_vec[1]-1;istk++)  {
    if (fnet_stack[istk].nFilament==0) continue; 
    kstk=istk+1;
    for (iFil=1;iFil<=fnet_stack[istk].nFilament;iFil++) {
      nVec=fnet_stack[istk].findFilamentSeq(iFil,filSeq);
      nVec=abs(nVec); // for cyclic
      for (i=1;i<nVec-1;i++) {
      nxtBead:
	//check if 3 beads in istk and 2 in istk+1 
	if (bead_quad[istk].count(filSeq[i])==0) {
	  jBead=filSeq[i-1]; kBead=filSeq[i+1];
	  bead1=jBead; bead2=kBead;bead_mid=filSeq[i]; //save middle bead
	  jt_range=it_fnet3d->bond_a[istk].equal_range(bead1);
	  kt_range=it_fnet3d->bond_a[istk].equal_range(bead2);
	  if (it_fnet3d->bond_a[istk].find(bead1)==it_fnet3d->bond_a[istk].end() ||
	      it_fnet3d->bond_a[istk].find(bead2)==it_fnet3d->bond_a[istk].end())
	    continue; // avoid seg fault
	    
	  for (jt=jt_range.first;jt!=jt_range.second;jt++) {
	    bead4=jt->second;
	    zt_range=fnet_stack[kstk].bond.equal_range(bead4);
	    for (kt=kt_range.first;kt!=kt_range.second;kt++) {
	      bead3=kt->second;
	      for (zt=zt_range.first;zt!=zt_range.second;zt++) {
		if (zt->second!=bead3) continue;
		quadbeads.at(0)=iQuad; //declare quad 1
		quadbeads.at(1)=bead1;quadbeads.at(2)=bead_mid;
		quadbeads.at(3)=bead3;quadbeads.at(4)=bead4;		
		if (!checkDuplicateQuad(istk,quadbeads)) {
		  addOneQuad(istk,quadbeads);
		  bead_quad[istk].insert(pair<int,int>(bead_mid,iQuad));
		  iQuad++; //i++; goto nxtBead;
		}
		
		quadbeads.at(0)=iQuad; //declare quad 2
		quadbeads.at(1)=bead_mid;quadbeads.at(2)=bead2;
		quadbeads.at(3)=bead3;quadbeads.at(4)=bead4;
		if (!checkDuplicateQuad(istk,quadbeads)) {
		  addOneQuad(istk,quadbeads);
		  bead_quad[istk].insert(pair<int,int>(bead_mid,iQuad)); 
		  iQuad++;i++; goto nxtBead;
		}
		
	      } //     for (zt=zt_range.first;zt!=zt_range.second;zt++) {
	    } //     for (kt=kt_range.first;kt!=kt_range.second;kt++) {
	  } //	  for (jt=jt_range.first;jt!=jt_range.second;jt++) {
	} // 	if (bead_quad[istk].count(filSeq[i])==0) {
	else {
	  qt_range=bead_quad[istk].equal_range(filSeq[i]);
	  for (qt=qt_range.first;qt!=qt_range.second;qt++) {
	    for (udum=0;udum<quad[istk].size();udum++) {
	      //not bead1 or bead2 of quad
	      if (quad[istk][udum].at(0)!=qt->second)   { 
		jBead=filSeq[i-1]; kBead=filSeq[i+1];
		bead1=jBead; bead2=kBead; bead_mid=filSeq[i]; //save middle bead
		jt_range=it_fnet3d->bond_a[istk].equal_range(bead1);
		kt_range=it_fnet3d->bond_a[istk].equal_range(bead2); 
		for (jt=jt_range.first;jt!=jt_range.second;jt++) {
		  bead4=jt->second;
		  zt_range=fnet_stack[kstk].bond.equal_range(bead4);
		  for (zt=zt_range.first;zt!=zt_range.second;zt++) {
		    bead3=zt->second;
		    for (kt=kt_range.first;kt!=kt_range.second;kt++) {
		      if (kt->second !=bead3) continue;
		      quadbeads.at(0)=iQuad; //declare quad 1
		      quadbeads.at(1)=bead1;quadbeads.at(2)=bead_mid;
		      quadbeads.at(3)=bead3;quadbeads.at(4)=bead4;
		      if (!checkDuplicateQuad(istk,quadbeads)) {
			addOneQuad(istk,quadbeads);
			bead_quad[istk].insert(pair<int,int>(bead_mid,iQuad));
			iQuad++; 
		      }
		      
		      quadbeads.at(0)=iQuad; //declare quad 2
		      quadbeads.at(1)=bead_mid;quadbeads.at(2)=bead2;
		      quadbeads.at(3)=bead3;quadbeads.at(4)=bead4;
		      if (!checkDuplicateQuad(istk,quadbeads)) {
		      	addOneQuad(istk,quadbeads);
		      	bead_quad[istk].insert(pair<int,int>(bead_mid,iQuad));
		      	iQuad++;i++; goto nxtBead;
		      }
		    } //    for (kt=kt_range.first;kt!=kt_range.second;kt++) {
		  } // 		  for (zt=zt_range.first;zt!=zt_range.second;zt++) {
		} //		for (jt=jt_range.first;jt!=jt_range.second;jt++) {
	      } //	      if (quad[istk][udum].at(0)!=qt->second)   { 
	    } // 	    for (udum=0;udum<quad[istk].size();udum++) {
	  } // 	  for (qt=qt_range.first;qt!=qt_range.second;qt++) {
	} // else {
	//check if 2 beads in istk and 3 in istk+1
	jBead=filSeq[i]; kBead=filSeq[i+1]; 	

	bead1=jBead; bead2=kBead;
	jt_range=it_fnet3d->bond_a[istk].equal_range(bead1);
	kt_range=it_fnet3d->bond_a[istk].equal_range(bead2);
	for (jt=jt_range.first;jt!=jt_range.second;jt++) {
	  bead4=jt->second;
	  for (kt=kt_range.first;kt!=kt_range.second;kt++) {
	    bead3=kt->second;  
	    zt_range=fnet_stack[kstk].bond.equal_range(bead4);
	    zt0_range=fnet_stack[kstk].bond.equal_range(bead3);
	    for (zt=zt_range.first;zt!=zt_range.second;zt++) {
	      for (zt0=zt0_range.first;zt0!=zt0_range.second;zt0++) {
		if (zt->second!=zt0->second) continue;
		bead_mid=zt->second;
		quadbeads.at(0)=iQuad; //declare quad 1
		quadbeads.at(1)=bead1;quadbeads.at(2)=bead2;
		quadbeads.at(3)=bead3;quadbeads.at(4)=bead_mid;
		if (!checkDuplicateQuad(istk,quadbeads)) {
		  addOneQuad(istk,quadbeads);
		  bead_quad[kstk].insert(pair<int,int>(bead_mid,iQuad));
		  iQuad++; 
		}
		
		quadbeads.at(0)=iQuad; //declare quad 2
		quadbeads.at(1)=bead1;quadbeads.at(2)=bead2;
		quadbeads.at(3)=bead_mid;quadbeads.at(4)=bead4;
		if (!checkDuplicateQuad(istk,quadbeads)) {
		  addOneQuad(istk,quadbeads);
		  bead_quad[kstk].insert(pair<int,int>(bead_mid,iQuad));
		  iQuad++;i++; 	  goto nxtBead;
		}
	      } //     for (zt0=zt0_range.first;zt0!=zt0_range.second;zt0++) { 	 
	    } // 	    for (zt=zt_range.first;zt!=zt_range.second;zt++) {
	  } // 	  for (kt=kt_range.first;kt!=kt_range.second;kt++) {
	} // 	for (jt=jt_range.first;jt!=jt_range.second;jt++) {
      } //       for (i=0;i<nVec;i++) {
    } //     for (iFil=1;iFil<=fnet_stack[istk].nFilament;iFil++) {
  } //   for (istk=fn3d_id_vec[0];istk<fn3d_id_vec[1];istk++)  {
  
  cout<<"  register_pentagon: "<<nQuad-nQuadold<<
    " quadrilaterals added"<<endl;
  delete [] filSeq;
}

/************************************************************************/
void surface::registerQuads()
{ /* For an image stack,if two in-plane bonded beads in one image (image i)
     are directly bonded (z-bond) to image i+1, register as quadrilateral. 
     Logic: check beads in plane i that are bonded to plane i+1 (stored in 
     bond_a multimap). If these are directly bonded to a bead in plane i+1 
     that is in-plane bonded to another bead bonded back to plane i 
     (stored in bond_b multimap), and if this bead is = original bead, 
     register as quadrilateral. Bead1 and bead2 are in plane i, bead3 and 
     bead4 in plane i+1.  */

  multimap<int,int>::iterator it,bt,zt,zt0,kt;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> bt_range,
    zt_range,zt0_range,kt_range;
  int istk,kstk,iQuad=nQuad;
  int bead1,bead2,bead3,bead4;
  array<int,5> quadbeads;
  vector<fnet3d>::iterator it_fnet3d;
  it_fnet3d=findFnet3d(&vfnet3d,this->fnet3d_tag);
  if (it_fnet3d==vfnet3d.end())
    error_dodri("Fnet3d tag not found."); 

  //function should be run without prior quad assignment 
  assert(iQuad==0); 
  for (istk=fn3d_id_vec[0];istk<fn3d_id_vec[1];istk++) {
    // index for the next plane: istk is the plane # (ie image in stack)
    kstk=istk+1;
     
    ////check in-plane bonds that each bead (i,j) is attached to
    for (it=it_fnet3d->bond_a[istk].begin();it!=it_fnet3d->bond_a[istk].end();it++) {
      bead1=it->first;
      bt_range=fnet_stack[istk].bond.equal_range(bead1);
      for (bt=bt_range.first;bt!=bt_range.second;bt++) {
	if (bt->second < bead1) continue;
	bead2=bt->second;
	zt_range=it_fnet3d->bond_a[istk].equal_range(bead2);
	zt0_range=it_fnet3d->bond_a[istk].equal_range(bead1);
	for (zt=zt_range.first;zt!=zt_range.second;zt++) {
	  bead3=zt->second;
	  for (zt0=zt0_range.first;zt0!=zt0_range.second;zt0++) {
	    bead4=zt0->second;
	    kt_range=fnet_stack[kstk].bond.equal_range(bead3);
	    for (kt=kt_range.first;kt!=kt_range.second;kt++) {
	      if (kt->second==bead4) {
		quadbeads.at(0)=iQuad; 
		quadbeads.at(1)=bead1;quadbeads.at(2)=bead2;
		quadbeads.at(3)=bead3;quadbeads.at(4)=bead4;
		addOneQuad(istk,quadbeads); iQuad++;
	      }
	    }//for (kt=kt_range.first;kt!=kt_range.second;kt++) 
	  }//  for (zt0=zt0_range.first;zt0!=zt0_range.second;zt0++) 
	}// for (zt=zt_range.first;zt!=zt_range.second;zt++) 
      }//for (bt=bt_range.first;bt!=bt_range.second;bt++)
    }//  for (it=it_fnet3d->bond_a[istk].begin();it!=it_fnet3d->bond_a[istk].end();it++) 
  }//for (istk=fn3d_id_vec[0];istk<fn3d_id_vec[1];istk++) {

  nQuad=iQuad;
  cout<<"  register_quads:: "<<nQuad<<" quadrilaterals formed"<<endl;
}

/************************************************************************/
// Non-member functions
/************************************************************************/

void checkDuplicateSurfaceTag(vector<surface> * vsurface, string tag)
{ /* In vfnet3d, find iterator for the element with tag */
  vector<surface>::iterator it_surface;
  for (it_surface=(*vsurface).begin();
       it_surface!=(*vsurface).end();++it_surface) {
    if ( (*it_surface).tag==tag)
      error_dodri("Duplicate surface tag: " + tag);
  }
  return;
}

/************************************************************************/
vector<surface>::iterator findSurface(vector<surface> *vsurface,string tag)
{ /* In vsurface, find iterator for the element with tag. */ 
  vector<surface>::iterator it_surface;
  for (it_surface=(*vsurface).begin();
       it_surface!=(*vsurface).end();++it_surface) {
    if ( (*it_surface).tag==tag) break;
  }
  return it_surface;
}
/************************************************************************/
