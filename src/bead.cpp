#include "dodri.h"

/************************************************************************/
// Constructors
/************************************************************************/

bead::bead(const bead &b0,double **cm1, double *beadMass1,string tag1)
{ /* Bead copy constructor. Creates bead object identical to input bead, 
     unless cm1 is different. New positions and beadMass passed to function.
     iBead for b0.cm and cm1 must be same (not checked).  */

  int j,iBead,iCell; 
  multimap<int,int> cell_bead_tmp, bead_cell_tmp;

  tag = tag1; //tag of this bead
  img_tag = b0.img_tag; //original image tag
  maxBead = b0.maxBead; 
  nBead = b0.nBead; beadDiam = b0.beadDiam;
  zCoord = b0.zCoord;
  rows = b0.rows; columns = b0.columns;
  dRow = b0.dRow; dCol = b0.dCol;
  
  //Initialize arrays
  cm=new double*[maxBead]; beadMass=new double[maxBead];
  for(iBead=0;iBead<(int)maxBead;iBead++) cm[iBead]=new double[2];
  
  for (iBead=0;iBead<nBead;iBead++) {
    for (j=0;j<2;j++)  cm[iBead][j]=cm1[iBead][j];
    beadMass[iBead]=beadMass1[iBead]; 
  }
  
  //re-calculate cell number in case of major movement 
  for (iBead=0;iBead<nBead;iBead++) {    
    iCell=nCell*(int)(cm[iBead][0]/dRow)+(int)(cm[iBead][1]/dCol); 
    cell_bead_tmp.insert(pair<int,int>(iCell,iBead));
    bead_cell_tmp.insert(pair<int,int>(iBead,iCell));
  }

  cell_bead.clear();   bead_cell.clear();
  cell_bead=cell_bead_tmp;
  bead_cell=bead_cell_tmp;

  cout<<"  Created bead "<<tag1<<" based on bead "<< b0.tag<<endl;  
  cout<<"  maxBead= "<<maxBead<<"  nBead= "<<nBead<<endl;
}

/************************************************************************/
bead::bead(img *img0, string tag0, int iDiam, string tag1, int dir)
{  /* Area-based constructor: from img0, assign beads of size iDiam.
      tag1: if not empty, save modified grayscale image after bead assignment
      to another img tag1. If empty, img0->pxl_gray remains.   */
  
  img_tag=img0->tag; // the original image tag
  tag=tag0; // tag of this bead
  nBead=0; beadDiam=iDiam;
  rows=img0->rows; columns=img0->columns;
  dRow=img0->dRow; dCol=img0->dCol;
  zCoord=0.; 

  int i,j,irow,icol;
  int rStart=0, cStart=0, rEnd=0, cEnd=0, rDelta=0, cDelta=0; //initialize
  int iBead=0, iArea=iDiam*iDiam;
  int rbox,cbox,box_row,box_col, pxl_cnt=0;
  double cmr=0.,cmc=0.,mass=0; //newBead assignment
  unsigned char **pxl_tmp=new unsigned char*[rows];
  string sdum=(tag1=="")?tag0:tag1;
  cout<<"  area build: iDiam= "<<iDiam<<" save= "<<sdum<<
      " dir= "<<dir<<endl;

  //direction-based variables
  std::list<int> plist[4]; // similar to pcount[*], list of pixel intensities
  std::list<int>::iterator it;
  
  for (irow=0;irow< rows;irow++) {
    pxl_tmp[irow]=new unsigned char[columns];
  }
  maxBead=0;
  for (irow=0;irow<rows; irow++) { // copy pxl_gray to pxl_tmp
    for (icol=0;icol<columns;icol++) {
      pxl_tmp[irow][icol]=img0->pxl_gray[irow][icol]; 
      if (pxl_tmp[irow][icol]>0) ++maxBead;
    }
  }
  // Initialize arrays
  maxBead = maxBead/ iDiam;
  cm=new double*[maxBead]; beadMass=new double[maxBead];
  for(i=0;i<(int)maxBead;i++) cm[i]=new double[2];
  
  setScanDir(iDiam,dir,rStart,cStart,rEnd,cEnd,rDelta,cDelta);
  rbox=rDelta*iDiam; cbox=cDelta*iDiam;
  
  //area-based bead assignment check  
  for (irow=rStart;irow!=rEnd; irow+=rDelta) {
    for (icol=cStart;icol!=cEnd;icol+=cDelta) { // newBead modifies p
      if (pxl_tmp[irow][icol]>0) {
	iBead= nBead; // iBead counts from previously found number
	box_row=irow+rbox; box_col=icol+cbox;
	cmr=0.; cmc=0.; mass=0.; pxl_cnt=0; 	
	
	for (i=irow;i!=box_row;i+=rDelta) {
	  for (j=icol;j!=box_col;j+=cDelta) {
	    if (pxl_tmp[i][j]>0) ++pxl_cnt;
	  }
	}	

	// register new bead if more than half of the box pixels are occupied
	if (pxl_cnt>=(iArea/2)) {
	  for (i=irow;i!=box_row;i+=rDelta) {
	    for (j=icol;j!=box_col;j+=cDelta) {
	      mass+=(double)pxl_tmp[i][j];
	      cmr+=(double)pxl_tmp[i][j] * (double)i;
	      cmc+=(double)pxl_tmp[i][j] * (double)j;
	      pxl_tmp[i][j]=0; // exclude this pixel from further consideration 
	    }
	  }
	  cmr/=mass; cmc/=mass;
	  mass/=(double)iArea;
	  createBead(mass,cmr,cmc);
	  if (nBead>=(int)maxBead) {
	    cout<< "bead::newBead: iBead "<<iBead<< " too large."<<
	      " Increase maxBead (current= "<<maxBead<<" )."<<endl;
	    exit(-1);
	  }
	} //	if (pxl_cnt>=(iArea/2)) {
      } //     if (pxl_tmp[irow][icol]>0) {	
    } //     for (icol=cStart;icol!=cEnd;icol+=cDelta) { // newBead modifies p
  } //  for (irow=rStart;irow!=rEnd; irow+=rDelta) {
  
  if (tag1!="") { // copy the result to a new image with tag1
    cout<<"  Copying result to new image tag"<<endl;
    vimg.push_back(img(*img0, pxl_tmp, tag1));
  }
  
  if (tag1=="") { // do not delete if creating another img.
    for (irow=0;irow<rows; irow++)
      delete [] pxl_tmp[irow]; // free up memory
    delete [] pxl_tmp;
  }
  cout<<"  maxBead= "<<maxBead<<"  nBead= "<<nBead<<endl;  
}

/************************************************************************/
bead::bead(img *img0, string tag0, int iDiam,string tag1,
	   double threshold, double loCut, int flag,  //for center
	   int hiCut, int delCut, double wFrac, //for boundary
	   int verbosity,string mode)
{
  /* Multi-dir constructor for all bead build modes. 
     
     Generate temporary beads based on img0 for scanning 4 directions,
     combine them into a new bead, then erase temporary beads 
  */
  
  time_t timer;
  stringstream ss;
  int i; string sdum[4];
  vector<bead>::iterator it_bead;
  vector<bead> vbeadtemp;

  cout<<"  Finding beads in 4 scanning directions (temporary beads created)."
      <<endl;
  time(&timer); // get current time for temporary bead name
  for (i=0;i<4;i++) {
    ss.str(""); ss.clear();
    ss<<"tmp"<< timer<<"_"<<i;
    sdum[i]=ss.str();
    if (mode=="area")
      vbeadtemp.push_back(bead(img0,sdum[i],iDiam,"",i));
    //if (mode=="center") 
    //  vbeadtemp.push_back(bead(img0,sdum[i],iDiam,"",threshold,loCut,flag,i));
    //if (mode=="boundary")
    //  vbeadtemp.push_back(bead(img0,sdum[i],iDiam,"",hiCut,delCut,wFrac,i));
    else error_dodri("Invalid bead build mode."); 
  }
  
  it_bead=findBead(&vbeadtemp,sdum[0]);
  for (i=1;i<4;i++) 
    it_bead->merge(&(*findBead(&vbeadtemp,sdum[i])),0.8,verbosity);
  
  //copy vbeadtemp into global vbead
  copyBead( (*it_bead),it_bead->cm,it_bead->beadMass,tag0);
  vbeadtemp.erase(vbeadtemp.begin(),vbeadtemp.end()); 
  
  //print built bead info
  it_bead=findBead(&vbead,tag0);
  if (it_bead==vbead.end())
    error_dodri(" Multi-directional bead build failed.");
}

/************************************************************************/
bead::bead(string tag, string ifname, string image_tag0, ifstream *fin)
{ /* Constructs a bead from the data in a .dat file with the name ifname */
  this->tag = tag;
  this->readData(ifname, image_tag0, fin);
}

/************************************************************************/
// Member functions
/************************************************************************/

void bead::createBead(double mass, double cmr, double cmc)
{
  cm[nBead][0]=cmr; cm[nBead][1]=cmc;
  beadMass[nBead]=mass;
  int iCell=(int)(nCell*(int)(cmr/dRow)+(int)(cmc/dCol));
  cell_bead.insert( pair<int,int>(iCell,nBead) );
  bead_cell.insert( pair<int,int>(nBead,iCell) );
  nBead++; 
}

/************************************************************************/
int bead::decimate(double rCut,int verbosity)
{
  /*Combines beads that are within rCut*beadDiam as a weighted average of
    the bead. It uses withinRadius to find nearby beads, sorted in order of
    distance to the current bead (nearest to farthest). The nearest bead is
    merged with the current bead, and the merged beads are placed back in
    the set of beads that can be matched by other beads. 
    
    Returns: number of beads removed.
    Verbosity flag: amount of information printed - 1 - max, 0 - min
  */

  int i, j, n = 0, closest, idum, numRemoved = 1,
    iCell, startBeads = nBead, curN;
  bool flag;
  double r = rCut*beadDiam, mass_tot, *beadMass_tmp, **cm_tmp;
  vector<pair<double, int> > closeBeads;
  multimap<int, int> usedPairs, oldToNew, bead_cell_tmp, cell_bead_tmp;
  pair<multimap<int, int>::iterator, multimap<int, int>::iterator> range;
  multimap<int, int>::iterator it;
  pair<int,int> bdPair;

  while (numRemoved > 0) { // Repeat decimation until no beads are removed
    n = 0;  // initialize
    bead_cell_tmp.clear(); cell_bead_tmp.clear();
    usedPairs.clear();     oldToNew.clear();
    cm_tmp = new double*[maxBead]; beadMass_tmp = new double[maxBead];
    
    for (i = 0; i < (int)maxBead; i++) cm_tmp[i] = new double[2];
    
    for (i = 0; i < nBead; i++) {
      // get beads within radius given by rCut sorted in order of
      // distance (nearest to farthest). Returns number of matches
      closeBeads.clear();
      idum = withinRadius(i, r, closeBeads);

      // direct copy of bead if no matches within radius
      if (!idum) { 
	for (j = 0; j < 2; j++) cm_tmp[n][j] = cm[i][j];
	beadMass_tmp[n] = beadMass[i];
	iCell=bead_cell.find(i)->second;
	cell_bead_tmp.insert(pair<int, int>(iCell, n));
	bead_cell_tmp.insert(pair<int, int>(n,iCell));
	++n;
	continue;
      }

      // Decimation with replacement
      closest = closeBeads.at(0).second;
      bdPair.first = closest; bdPair.second = i;

      // Check to see if pair has already been merged
      // This may happen due to i<->closest switching
      flag = false; 
      range = usedPairs.equal_range(i);
      for (it = range.first; it != range.second; it++) {
	if (it->second == closest) { 
	  flag = true;  break;
	}
      }
      if (flag) continue;
      usedPairs.insert(bdPair);
      
      it = oldToNew.find(closest);
      if (it!=oldToNew.end()) { // case when closest is previously used
	curN = --n;             // after removing previous closest, this may be 
	n = it->second;         // the next closest: reduce n and search again
      }
      else { 
	curN = n;
	bdPair.first = i; bdPair.second = n;
	oldToNew.insert(bdPair);
      }

      // Weighted average of the pair of beads
      mass_tot = beadMass[i] + beadMass[closest];
      for (j = 0; j < 2; j++) {
	cm_tmp[n][j] = (beadMass[i]*cm[i][j] +
			beadMass[closest]*cm[closest][j]) / mass_tot;
	cm[i][j] = cm[closest][j] = cm_tmp[n][j];
      }
      beadMass_tmp[n] = (beadMass[i]*beadMass[i] +
			 beadMass[closest]*beadMass[closest]) / mass_tot;

      beadMass[i] = beadMass[closest] = beadMass[n];

      // Calculating new cell for bead and recording that information
      iCell=(int)(nCell*(int)(cm_tmp[n][0]/dRow)+(int)(cm_tmp[n][1]/dCol));

      // Removes bead-cell data from bead_cell_tmp and cell_bead_tmp
      // in the case that a bead that has already been merged is
      // merged again. Prevents duplicates of bead info and/or
      // conflicting info for a merged bead.
      it = bead_cell_tmp.find(n);
      if (it != bead_cell_tmp.end()) {
	idum = it->second;
	bead_cell_tmp.erase(it);
	range = cell_bead_tmp.equal_range(idum);
	for (it = range.first; it != range.second; it++) {
	  if (it->second == n) {
	    cell_bead_tmp.erase(it);
	    break;
	  }
	}
      }
      cell_bead_tmp.insert(pair<int, int>(iCell, n));
      bead_cell_tmp.insert(pair<int, int>(n, iCell));
      n = curN;
      ++n;
    } // for (i = 0; i < nBead; i++) {
    
    // Deleting and replacing old bead data with new bead data
    for (i = 0; i < (int)maxBead; i++)  delete[] cm[i];
    delete[] cm; delete[] beadMass;
    
    cm = cm_tmp;    
    beadMass = beadMass_tmp;
    bead_cell.clear();
    bead_cell = bead_cell_tmp;
    cell_bead.clear();
    cell_bead = cell_bead_tmp;
    numRemoved = nBead - n;
    nBead = n;
  } // while (numRemoved > 0) {
  
  numRemoved = startBeads - n;
  
  if (verbosity)
    cout << "  decimate: rCut= "<<rCut<<"  Number of beads removed= "
	 <<numRemoved<<"  new nBead= "<<nBead<<endl;
  
  return numRemoved;
}

//************************************************************************/
double bead::getAngleOffset()
{ /* Get angle offset from x-axis of bead data.   */

  int i,j;
  double bm[2],sig, pos[2][nBead];
  double Rx,Ry,Rxy,theta,rad[2]; 
  
  // save cm in [dim][nBead] format (beadMass is same for x and y
  for (i=0;i<nBead;i++) {
    for (j=0;j<2;j++) pos[j][i]=cm[i][j];
  }
  
  for (j=0;j<2;j++) //calc weighted average (centroid) of body
    getavgw(&bm[j],&sig,pos[j],beadMass,nBead);
  
  Rx=0;Ry=0;Rxy=0;
  for (i=0;i<nBead;i++){
    for (j=0;j<2;j++)  rad[j]=cm[i][j]-bm[j];
    rad[j]/=sqrt(rad[0]*rad[0]+rad[1]*rad[1]);
    Rx+=(rad[0]*rad[0]);
    Ry+=(rad[1]*rad[1]);
    Rxy+=(rad[0]*rad[1]);
  }
  
  theta=atan((2.0*Rxy)/(Ry-Rx))/2.0;
  return theta; 
}

/************************************************************************/
int bead::getBeadDir(double rCut, multimap<int,array<double,2>> &cm_dir,
		     multimap<int,array<double,2>> &cm_buf)
{ /* For each bead, look for pixels whose distance are within
     rCut*beadDiam, get pixel intensity-weighted center of mass, direction
     (via lin_fit), rmsd from the linear fit line. -> save the center of mass
     coord to cm_buf, slope (in radians) and chi2 to cm_dir.
     
     returns: number of beads whose directions are calculated
     
     cm_buf: fit coord
     cm_dir[*]: [0]: stores direction angle (in rad), [1]: chi2  */
  
  // maskMax: max number of pixels in the mask
  // maskSize: number of pixels to consider in the mask
  int ibead,i,ix,iy, maskSize, maskMax, nBeadDir=0;
  int pxl_range[2][2]; // pxl_range[x/y][min/max]
  unsigned char cdum;
  double maskRad=rCut*(double)beadDiam; // maskRad: mask radius
  double maskRadSq=maskRad*maskRad, r0[2],rdum, vdum[2], **mask;
  double x0,y0,a0,b0,siga,sigb,chi2,sx,sy;
  double *rvec=new double[nBead];
  array<double,2> cm_dir_arr, cm_buf_arr; 
  
  vector<img>::iterator it_img;
  
  it_img=findImg(img_tag);
  if (it_img==vimg.end())
    error_dodri("getBeadDir: cannot find image "+ img_tag);
  cout<<"  getBeadDir: Mask radius= "<< maskRad <<endl;
  
  i=(int)ceil(maskRad); maskMax=4*i*i; 
  mask=new double*[3];
  for (i=0;i<3;i++) mask[i]=new double[maskMax];
  
  cm_dir.clear(); cm_buf.clear(); //clear previous cm_dir, cm_buf multimap 
  for (ibead=0;ibead<nBead;ibead++) {
    // Default values for cm_buf, cm_dir
    for (i=0;i<2;i++) cm_buf_arr[i]=cm[ibead][i];
    cm_dir_arr[0]=0.; cm_dir_arr[1]=-1; 
    cm_dir.insert(pair<int,array<double,2>>(ibead,cm_dir_arr));
    cm_buf.insert(pair<int,array<double,2>>(ibead,cm_buf_arr));
    
    // setup mask
    for (i=0;i<2;i++) {
      r0[i]=cm[ibead][i];
      pxl_range[i][0]=(int)floor(r0[i]-maskRad);
      if (pxl_range[i][0]<0) pxl_range[i][0]=0;
      pxl_range[i][1]=(int)ceil(r0[i]+maskRad);
    }
    if (pxl_range[0][1]>rows) pxl_range[0][1]=rows;
    if (pxl_range[1][1]>columns) pxl_range[1][1]=columns;

    maskSize=0;
    for (ix=pxl_range[0][0];ix<pxl_range[0][1];ix++) {
      for (iy=pxl_range[1][0];iy<pxl_range[1][1];iy++) {
	vdum[0]=r0[0]-(double)ix; vdum[1]=r0[1]-(double)iy;
	rdum=dotprod(vdum,vdum,2);
	if (rdum<maskRadSq) { // within mask
	  cdum=it_img->pxl_gray[ix][iy];
	  if (cdum>0) { // only count nonzero pixel intensity
	    mask[0][maskSize]=(double)ix;
	    mask[1][maskSize]=(double)iy;
	    mask[2][maskSize]=(double)cdum;
	    ++maskSize;
	  }
	}
      }
    } // for (ix=pxl_range[0][0];ix<pxl_range[0][1];ix++) {
    
    // perform linear fit
    if (maskSize<2) continue; // too small number of pixels
    fit1(mask[0],mask[1],maskSize,mask[2],1,&a0,&b0,&siga,&sigb,&chi2);
    if (fabs(b0)>1) { // larger slope. fit to x=a+by
      fit1(mask[1],mask[0],maskSize,mask[2],1,&a0,&b0,&siga,&sigb,&chi2);
      b0=1.0/b0;
    }
    getavgw(&x0,&sx,mask[0],mask[2],maskSize);
    getavgw(&y0,&sy,mask[1],mask[2],maskSize);
    cm_buf.find(ibead)->second.at(0)=x0;
    cm_buf.find(ibead)->second.at(1)=y0;
    cm_dir.find(ibead)->second.at(0)=atan(b0);
    cm_dir.find(ibead)->second.at(1)=chi2;
    nBeadDir++;
  } //   for (ibead=0;ibead<nBead;ibead++) {
  cout<<"    Number of bead directions calculated: "<<nBeadDir<<endl;

  // find min/max of chi2
  i=0;
  for (ibead=0;ibead<nBead;ibead++) {
    if (cm_dir.find(ibead)->second.at(1)>0.)
      rvec[i++]=cm_dir.find(ibead)->second.at(1);
  }
  sx=0;sy=0; //re-set 
  MinMax(i,rvec,&ix,&iy,&sx,&sy);
  cout<<"    Min/Max of chi^2: "<<sx<<"  "<<sy<<endl;

  delete rvec;
  for (i=0;i<3;i++) delete [] mask[i];
  delete [] mask;
  
  return nBeadDir;  
}

/************************************************************************/
void bead::gradient(double rCut,int iRad,
		    multimap<int,array<double,2>> &cm_dir)
{ /* Calculates gradient at the location of each bead, based on img_tag.
     If img_tag0 is given, use that image instead of this->img_tag.     

     Call this method "Half-moon gradient": 

     For a given bead, consider pixels within a circle of radius iRad.
     Divide the circle in half for a given direction, and calculate average
     pixel intensity in each half circle. Set the direction that yields the
     greatest pixel intensity difference between half circles, to be
     perpendicular to the gradient direction.

     Due to the discreetenss of pixels, for a given iRad, only a finite number
     of directions are considered. When calculating the avg pixel intensity,
     discard pixels that lie on the directional line dividing the circle.
     
     cut_dir: stores pixels at the periphery of the circle in 1st quadrant.
     for negative slopes, switch sign of its y-value, scanning the range
     [-pi/2,pi/2]. These points are used for the line dividing the circle.

     mask0/1: pixel masks for two half-circles for a given direction.

     60 in mask: max number of directions to consider
     When all 60 dirs are used, angular resolution of gradient is 180/60=3
     degrees.  The actual number of elements in mas0/1 to use is twice 
     cut_dir.size().

     storing points in mask0/1: 
     if slope from cut_dir is greater than or equal to 1 in magnitude 
     (pi/4; steep line), store pixels on the left/right into mark0/1
     if the slope is less than 1 in magnitude, store pixels below/above 
     the line to mark0/1.
*/

  int i,j,ix,iy,irow,icol,iBead,dcount,idir,idum;
  int cx,cy,imax,pmax,pmax0,iBead0, **contrast;
  int static iRad2=iRad*iRad;
  double s1; // slope of the dividing line
  double rdum, cm0[2][2];
  string sdum;
  vector< pair<int,int> > cut_dir, mask0[60],mask1[60];
  vector< pair<int,int> >::iterator it;
  short **mask_circle=new short*[2*iRad+1]; // =1 if inside circle, 0: if not.

  multimap<int,array<double,2>> cm_buf;
  getBeadDir(rCut,cm_dir,cm_buf);   //find bead direction

  // First determine directions of cut
  if (iRad<=0)
    error_dodri("Window radius for gradient must be a positive integer");
  
  vector<img>::iterator it_img;
  it_img=findImg(img_tag);
  if (it_img==vimg.end())
    error_dodri("getBeadDir: cannot find image "+ img_tag);
  if (nBead==0) error_dodri("No beads found in "+tag);
  cout<<"  gradient: radius= "<<iRad<<endl; 

  // Find cut direction in the first quadrant relative to (irow,icol)
  cut_dir.push_back(pair<int,int>(iRad,0)); // horizontal 
  i=1;
  for (ix=iRad-1;ix>=0;--ix) { // vary x-coord
    for (iy=i;iy<=iRad;++iy) { // y-coord
      j=ix*ix+iy*iy;
      if (j<=iRad2) cut_dir.push_back(pair<int,int>(ix,iy));
    }
    ++i; // look for the next level
  }
  
  // populate mask_circle:
  for (i=-iRad;i<=iRad;i++) {
    mask_circle[i+iRad]=new short[2*iRad+1];
    for (j=-iRad;j<=iRad;j++) {
      mask_circle[i+iRad][j+iRad]=((i*i+j*j)>iRad2)?0:1;
    }
  }

  /************ CREATE MASK ************/
  // Below, consider positive and negative cutting directions separately
  dcount=0; // index for cut_dir
  for (it=cut_dir.begin();it!=cut_dir.end();it++) { // positive direction
    cx=it->first; cy=it->second; // positive slope
    // idir=0/1: dividing line less/greather than pi/4 
    if (cx>cy) { s1=(double)cy/(double)cx; idir=0; }
    else { s1=(double)cx/(double)cy; idir=1;}
    for (i=-iRad;i<=iRad;i++) {
      for (j=-iRad;j<=iRad;j++) {
	if (i==0 && j==0) continue; // skip origin
	if (mask_circle[i+iRad][j+iRad]==0) continue; // outside the circle
	if (abs(i)>abs(j)) { // angle less than pi/4
	  if (idir==1) { // dir>=pi/4, (i,j)<pi/4
	    if (i>0) mask1[dcount].push_back(pair<int,int>(i,j)); // right
	    else mask0[dcount].push_back(pair<int,int>(i,j)); // left
	  }
	  else { // both dir & (i,j) < pi/4: consider 4 quadrants
	    if ((j==0) && (s1<TOL)) continue; // horizontal
	    if (j>=0) {
	      if (i<0) mask1[dcount].push_back(pair<int,int>(i,j)); // Q2: up
	      else { // Q1
		rdum=(double)j/(double)i; // i!=0 here
		if (rdum<(s1-TOL)) 
		  mask0[dcount].push_back(pair<int,int>(i,j)); // down
		else if (rdum>(s1+TOL))
		  mask1[dcount].push_back(pair<int,int>(i,j)); // up
		else continue; // aligns w/ cut_dir
	      }
	    }
	    else { // j<0
	      if (i>0) mask0[dcount].push_back(pair<int,int>(i,j)); // Q4: down
	      else { // Q3
		rdum=(double)j/(double)i; // i!=0 here
		if (rdum<(s1-TOL)) 
		  mask1[dcount].push_back(pair<int,int>(i,j)); // up
		else if (rdum>(s1+TOL))
		  mask0[dcount].push_back(pair<int,int>(i,j)); // down
		else continue;
	      }
	    } // else { // j<0
	  } // else { // both dir & (i,j) < pi/4
	} // if (abs(i)>abs(j)) {
	else { // abs(i)<=abs(j): angle > pi/4
	  if (idir==0) { // dir<pi/4
	      if (j>0) mask1[dcount].push_back(pair<int,int>(i,j)); // up
	      else mask0[dcount].push_back(pair<int,int>(i,j)); // down
	  }
	  else { // both dir & (i,j) > pi/4: consider 4 quadrants
	    if (i==0 && s1<TOL) continue; // vertical
	    if (j>=0) {
	      if (i<0) mask0[dcount].push_back(pair<int,int>(i,j)); // Q2: left
	      else { // Q1
		rdum=(double)i/(double)j;
		if (rdum<(s1-TOL)) 
		  mask0[dcount].push_back(pair<int,int>(i,j)); // left
		else if (rdum>(s1+TOL))
		  mask1[dcount].push_back(pair<int,int>(i,j)); // right
		else continue;
	      }
	    }
	    else { // j<0
	      if (i>0) mask1[dcount].push_back(pair<int,int>(i,j)); // Q4: right
	      else { // Q3
		rdum=(double)i/(double)j;
		if (rdum<(s1-TOL)) 
		  mask1[dcount].push_back(pair<int,int>(i,j)); // right
		else if (rdum>(s1+TOL))
		  mask0[dcount].push_back(pair<int,int>(i,j)); // left
		else continue;
	      }
	    } // else { // j<0
	  } // else { // both dir & (i,j) > pi/4: consider 4 quadrants
	} //	else { // abs(i)<=abs(j): angle > pi/4
      } // for (j=-iRad;j<=iRad;j++) {
    } // for (i=-iRad;i<=iRad;i++) {
    ++dcount;
  } // for (it=cut_dir.begin();it!=cut_dir.end();it++) { // positive direction
  idum=dcount; // number of positive slope elements (used for debug only)

  // dcount continues to increase
  for (it=cut_dir.begin();it!=cut_dir.end();it++) { // negative direction
    cx=it->first; cy=it->second; // cx,cy >0
    if (cx==0 || cy==0) continue; // skip horizontal/vertical
    // idir=0/1: dividing line less/greather than pi/4 
    if (cx>cy) { s1=-(double)cy/(double)cx; idir=0; }
    else { s1=-(double)cx/(double)cy; idir=1;}
    for (i=-iRad;i<=iRad;i++) {
      for (j=-iRad;j<=iRad;j++) {
	if (i==0 && j==0) continue; // skip origin
	if (mask_circle[i+iRad][j+iRad]==0) continue; // outside the circle
	if (abs(i)>abs(j)) { // angle less than pi/4
	  if (idir==1) { // dir>=pi/4, (i,j)<pi/4
	      if (i>0) mask1[dcount].push_back(pair<int,int>(i,j)); // right
	      else mask0[dcount].push_back(pair<int,int>(i,j)); // left
	  }
	  else { // both dir & (i,j) < pi/4: consider 4 quadrants
	    if (j==0 && fabs(s1)<TOL) continue; // horizontal
	    if (j>=0) {
	      if (i>0) mask1[dcount].push_back(pair<int,int>(i,j)); // Q1: up
	      else { // Q2
		rdum=(double)j/(double)i; // <0
		if (rdum>(s1+TOL)) 
		  mask0[dcount].push_back(pair<int,int>(i,j)); // down
		else if (rdum<(s1-TOL))
		  mask1[dcount].push_back(pair<int,int>(i,j)); // up
		else continue;
	      }
	    }
	    else { // j<0
	      if (i<0) mask0[dcount].push_back(pair<int,int>(i,j)); // Q3: down
	      else { // Q4
		rdum=(double)j/(double)i; // i!=0 here
		if (rdum>(s1+TOL)) 
		  mask1[dcount].push_back(pair<int,int>(i,j)); // up
		else if (rdum<(s1-TOL))
		  mask0[dcount].push_back(pair<int,int>(i,j)); // down
		else continue;
	      }
	    } // else { // j<0
	  } // else { // both dir & (i,j) < pi/4
	} // if (abs(i)>abs(j)) {
	else { // abs(i)<=abs(j): angle > pi/4
	  if (idir==0) { // dir<pi/4
	      if (j>0) mask1[dcount].push_back(pair<int,int>(i,j)); // up
	      else mask0[dcount].push_back(pair<int,int>(i,j)); // down
	  }
	  else { // both dir & (i,j) > pi/4: consider 4 quadrants
	    if (i==0 && fabs(s1)<TOL) continue; // vertical
	    if (j>=0) {
	      if (i>0) mask1[dcount].push_back(pair<int,int>(i,j)); // Q1: right
	      else { // Q2
		rdum=(double)i/(double)j;
		if (rdum<(s1-TOL)) 
		  mask0[dcount].push_back(pair<int,int>(i,j)); // left
		else if (rdum>(s1+TOL))
		  mask1[dcount].push_back(pair<int,int>(i,j)); // right
		else continue;
	      }
	    }
	    else { // j<0
	      if (i<0) mask0[dcount].push_back(pair<int,int>(i,j)); // Q3: left
	      else { // Q4
		rdum=(double)i/(double)j;
		if (rdum<(s1-TOL)) 
		  mask1[dcount].push_back(pair<int,int>(i,j)); // right
		else if (rdum>(s1+TOL)) 
		  mask0[dcount].push_back(pair<int,int>(i,j)); // left
		else continue;
	      }
	    } // else { // j<0
	  } // else { // both dir & (i,j) > pi/4: consider 4 quadrants
	} //	else { // abs(i)<=abs(j): angle > pi/4
      } // for (j=-iRad;j<=iRad;j++) {
    } // for (i=-iRad;i<=iRad;i++) {
    ++dcount;
  } // for (it=cut_dir.begin();it!=cut_dir.end();it++) { // negative direction
  //dneg=dcount-idum; // number of negative slope elements

  /************ DONE WITH CREATING MASK ************/
  contrast=new int*[dcount];
  for (i=0;i<dcount;i++) contrast[i]=new int[2];
  iBead0=-1; pmax0=0;
  
  for (iBead=0;iBead<nBead;iBead++) {
    irow=round(cm[iBead][0]); icol=round(cm[iBead][1]);
    assert(irow<=(int)rows); assert(icol<=(int)columns);
    // adjust if bead is too close to the boundary
    if (irow<iRad || icol<iRad ) { 
      cm_dir.find(iBead)->second.at(0)=-1.;             
      cm_dir.find(iBead)->second.at(1)=0;       
      continue;
    }
    i=irow+iRad; j=icol+iRad;
    if (i>=rows || j>=columns ) {
      cm_dir.find(iBead)->second.at(0)=1.;             
      cm_dir.find(iBead)->second.at(1)=0;       
      continue; 
    }
    
    for (i=0;i<dcount;i++) for (j=0;j<2;j++) contrast[i][j]=0;

    for (i=0;i<dcount;i++) { // sum of pixel intensities
      for (it=mask0[i].begin();it!=mask0[i].end();it++) {
	ix= it->first; iy= it->second;
	contrast[i][0]+=(int)it_img->pxl_gray[irow+ix][icol+iy];
      } 
      for (it=mask1[i].begin();it!=mask1[i].end();it++) {
	ix= it->first; iy= it->second;
	contrast[i][1]+=(int)it_img->pxl_gray[irow+ix][icol+iy];
      }
    } // for (i=0;i<dcount;i++) { // sum of pixel intensities
    imax=-1; pmax=0;
    for (i=0;i<dcount;i++) {
      j=abs(contrast[i][0]-contrast[i][1]);
      if (j>pmax) { pmax=j; imax=i; } // direction with max contrast
    }
    if (imax==-1) { // no contrast
      cout<<"  Zero intensity contrast for bead "<<iBead<<endl;
      cm_dir.find(iBead)->second.at(0)=cm_dir.find(iBead)->second.at(1)=0.;
      continue; 
    }
    
    if (pmax>pmax0) {pmax0=pmax; iBead0=iBead;} // global max
    
    // get center of masses
    for (i=0;i<2;i++) for (j=0;j<2;j++) cm0[i][j]=0.;
    idum=0;
    for (it=mask0[imax].begin();it!=mask0[imax].end();it++) {
	ix= irow+ it->first; iy= icol+ it->second;
	if (contrast[imax][0]>0) { // has nonzero pixel
	  cm0[0][0]+=(double)ix*(double)it_img->pxl_gray[ix][iy];
	  cm0[0][1]+=(double)iy*(double)it_img->pxl_gray[ix][iy];
	}
	else { // half moon is all zero; find centroid
	  cm0[0][0]+=(double)ix;  cm0[0][1]+=(double)iy;
	  ++idum;
	}
    }
    for (it=mask1[imax].begin();it!=mask1[imax].end();it++) {
	ix= irow+ it->first; iy= icol+ it->second;
	if (contrast[imax][1]>0) { // has nonzero pixel
	  cm0[1][0]+=(double)ix*(double)it_img->pxl_gray[ix][iy];
	  cm0[1][1]+=(double)iy*(double)it_img->pxl_gray[ix][iy];
	}
	else { // half moon is all zero; find centroid
	  cm0[1][0]+=(double)ix;  cm0[1][1]+=(double)iy;
	  ++idum;
	}
    }
    for (i=0;i<2;i++) { // mask0/1
      for (j=0;j<2;j++) { // x,y components
	if (contrast[imax][i]>0) cm0[i][j]/=(double)contrast[imax][i];
	else  cm0[i][j]/=(double)idum; // centroid
      }
    }
    rdum=getDistPos(cm0[0],cm0[1],0); // distance between cm
    for (i=0;i<2;i++) {
      cm0[1][i]-=cm0[0][i]; // change to displacement
      cm0[1][i]/=rdum; // normalize
    }
    // intensity difference opposite in dir to cm0[1]
    s1=(double)(contrast[imax][0]-contrast[imax][1]); 
    s1/=(double)iRad2; // divide by area
    for (i=0;i<2;i++) 
      cm_dir.find(iBead)->second.at(i)=s1*cm0[1][i]/rdum; // assign gradient
  } //   for (iBead=0;iBead<nBead;iBead++) {
  
  cout<<"    Max gradient for bead= "<<iBead0<<"  intensity contrast= "
      <<pmax0<<endl;
  cout<< "    Gradient= "<<cm_dir.find(iBead0)->second.at(1)<<" "
      <<-cm_dir.find(iBead0)->second.at(0)<< endl;
  
  
  for (i=0;i<(2*iRad+1);i++) delete [] mask_circle[i];
  for (i=0;i<dcount;i++) delete [] contrast[i];
  delete [] mask_circle; delete [] contrast;
  return;
}

/************************************************************************/
void bead::merge(bead *bd, double rCut,int verbosity)
{
  /* Combines the beads of bead b with the beads of this bead.  Only this
     bead is modified. If rCut>0, the decimate(rCut) is called to remove
     beads that are within rCut*beadDiam of one another and replace them
     with a bead or beads that are weighted averages of the ones removed.
  */
  
  int i, j, iCell, oldMaxBead = maxBead;
  double *beadMass_tmp, **cm_tmp;
  string sdum=tag;

  if (verbosity)
    cout<<"  merge_bead:  rCut= "<<rCut<<" merge= "<<bd->tag<<
      "  with= "<<tag<<"  save="<<sdum<<endl;
  
  maxBead = (maxBead>bd->maxBead)?maxBead:bd->maxBead;
  maxBead = (maxBead>(int)(nBead+bd->nBead))?maxBead:maxBead+bd->maxBead;
  
  cm_tmp = new double*[maxBead]; beadMass_tmp = new double[maxBead];  
  for (i=0;i<(int)maxBead;i++)     cm_tmp[i]=new double[2];

  // Copies data from both beads into new arrays
  for (i = 0; i < nBead + bd->nBead; i++) {
    if (i < nBead) { // Data from this bead
      for (j = 0; j < 2; j++) 	cm_tmp[i][j] = cm[i][j];
      beadMass_tmp[i] = beadMass[i];
    }
    else { // Data from bead b
      for (j = 0; j < 2; j++) 	cm_tmp[i][j] = bd->cm[i-nBead][j];
      beadMass_tmp[i] = 1.; //bd->beadMass[i-nBead];
      
      iCell=(int)(nCell*(int)(cm_tmp[i][0]/dRow)+(int)(cm_tmp[i][1]/dCol));
      cell_bead.insert( pair<int,int>(iCell,i) );
      bead_cell.insert( pair<int,int>(i,iCell) );
    }
  }

  for (i = 0; i < oldMaxBead; i++) delete[] cm[i];  
  delete[] cm;  
  delete[] beadMass;

  cm = cm_tmp;
  beadMass = beadMass_tmp;
  nBead += bd->nBead;
  if (rCut > 0) decimate(rCut,verbosity); 
  return;
}

/************************************************************************/
int bead::moveBeadPos(int p, double *newPos)
{ /* Moves position of bead p to newPos. Returns: cell number after move */
  
  int iCell, jCell;
  double x0,y0;
  multimap<int,int>::iterator it;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> it_range;

  cm[p][0]=x0=newPos[0]; cm[p][1]=y0=newPos[1];
  iCell=nCell*(int)(x0/dRow)+(int)(y0/dCol);
  
  it=bead_cell.find(p); jCell=it->second;
  if (jCell == iCell) return iCell; 
  else { // cell number changed
    bead_cell.erase(it);
    bead_cell.insert( pair<int,int>(p,iCell) );
    it_range=cell_bead.equal_range(jCell);
    for (it=it_range.first; it!=it_range.second; it++) {
      if (it->second ==p) { cell_bead.erase(it); break; }
    }
    cell_bead.insert( pair<int,int>(iCell,p) );
  }
  return iCell;
}

/************************************************************************/
void bead::readData(string fname, string image_tag0, ifstream *fin)
{ /* Reads in a bead data file to fill in the data fields 
     of the new bead object.   */
  
  stringstream ss;
  string image_tag, line;
  int i,ivec[1], bead, nBead_temp=0;
  double avec[1], mass, cmr, cmc;
  bool newFile = (fin == NULL);
  bool i_tag, mBead, diam, r, c, dr, dc, flag;
  vector<img>::iterator it_img;
  flag = i_tag = mBead = diam = r = c = dr = dc = false;
  
  if (img_tag!="") {
    cout << "WARNING: bead " << tag << " is being overwritten.\n";
    delete[] beadMass;
    for (i = 0; i < (int)maxBead; i++)  delete[] cm[i];
    delete[] cm;    
    bead_cell.clear();    cell_bead.clear();
  }
  
  nBead = 0;
  if (newFile) {
    fin = new ifstream;
    fin->open(fname.c_str(), ios::in);
  }
  
  if (fin->bad()) error_dodri("dodri: Invalid bead file");
  
  while(!flag && fin->good()) {
    
    flag = i_tag && mBead && diam && r && c && dr && dc;
    if (flag) break;
    
    if (line.find("Bead ")!=std::string::npos) error_dodri("READ:: Incomplete bead header");
    
    getline(*fin, line);
    if (line == "") continue;
    ss.str(line);
    
    if (line.find("nBeads=")!=string::npos) {
      getCommOptInt(ss,"nBeads=",ivec,1,-1);
      nBead_temp=ivec[0];
      if (nBead_temp==0) {
	cout<<"  "<<std::right<<"WARNING: file contains 0 beads."<<endl;
	return;
      }
    }      
    
    if (!i_tag && line.find("img_tag=")!=string::npos) {
      if (image_tag0=="")
	image_tag=getCommOptString(ss,"img_tag=","");
      else image_tag=image_tag0; 
      if (image_tag!="") {
	img_tag = image_tag;
	i_tag = true;
	continue;
      }
      else error_dodri("READ:: missing image tag");
    }
    it_img = findImg(img_tag);  //check if img_tag exists in vimg
    if (it_img==vimg.end()) error_dodri("READ:: image tag not found");
    
    if (!mBead && line.find("maxBead=")!=string::npos) {
      getCommOptInt(ss, "maxBead=", ivec, 1, -1);
      if (ivec[0] != -1) {
	mBead = true;
	maxBead = ivec[0];
	cm=new double*[maxBead]; beadMass=new double[maxBead];
	for (i=0;i<(int)maxBead;i++) 	  cm[i]=new double[2];
	continue;
      }
      else error_dodri("READ:: invalid or missing maxBead");
    }
    
    if (!diam && line.find("beadDiam=")!=string::npos) {
      getCommOptInt(ss, "beadDiam=", ivec, 1, -1);
      if (ivec[0] != -1) {
	beadDiam = ivec[0];
	diam = true;
	continue;
      }
      else error_dodri("READ:: invalid or missing beadDiam");
    }
    
    if (!r && line.find("rows=")!=string::npos) {
      getCommOptInt(ss, "rows=", ivec, 1, -1);
      if (ivec[0] != -1) {
	rows = ivec[0];
	r = true;
	continue;
      }
      else error_dodri("READ:: invalid or missing rows");
    }
    
    if (!c && line.find("columns=")!=string::npos) {
      getCommOptInt(ss, "columns=", ivec, 1, -1);
      if (ivec[0] != -1) {
	columns = ivec[0];
	c = true;
	continue;
      }
      else error_dodri("READ:: invalid or missing columns");
    }

    if (!dr && line.find("dRow=")!=string::npos) {
      getCommOptDouble(ss, "dRow=", avec, 1, -1);
      if (avec[0] != -1) {
	dRow = avec[0];
	dr = true;
	continue;
      }
      else error_dodri("READ:: invalid or missing dRow");
    }
    
    if (!dc && line.find("dCol=")!=string::npos) {
      getCommOptDouble(ss, "dCol=", avec, 1, -1);
      if (avec[0] != -1) {
	dCol = avec[0];
	dc = true;
	continue;
      }
      else error_dodri("READ:: invalid or missing dCol");
    }
  }
  
  while(!ss.fail() && getline(*fin, line)) {
    if (line == "" || line.find("Bead ")!=std::string::npos) continue;
    ss.str(line);
    ss >> bead >> mass >> cmr >> cmc >> zCoord;
    if (ss.fail()) error_dodri("READ:: error in reading beads");
    createBead(mass, cmr, cmc);
    if (!isdigit(fin->peek())) break;
  }
  
  cout<<"  read: maxBead= "<<maxBead<<"  nBead= "<<nBead<<endl;
  if (newFile) fin->close();
}

/************************************************************************/
void bead::rotate(double theta,double *pt)
{ /* Rotates this bead set by theta (in deg) about either centroid or pt[].
     Updates bead positions.  */
  
  //cout<<"  rotate: angle of rotation= "<<theta<<endl;
  //cout<<"   MESSAGE: bead " << tag << " will be modified.\n";
  
  int i,j;
  double pos[2][nBead],bm[2],sig; 
  double phi=theta*degree2rad,cm_temp[nBead][2]; 
  
  // save cm in [dim][nBead] for getavgw() 
  for (i=0;i<nBead;i++) {
    for (j=0;j<2;j++) pos[j][i]=cm[i][j];
  }
  
  // set point about which to rotate
  if (pt==NULL) {
    for (j=0;j<2;j++) //calc weighted average (centroid) of body
      getavgw(&bm[j],&sig,pos[j],beadMass,nBead);
  }
  else {
    for (j=0;j<2;j++) bm[j]=pt[j];
  }
  
  // translate bead pos s.t. centroid is at (0,0)
  for (i=0;i<nBead;i++) {
    for (j=0;j<2;j++)  pos[j][i]=pos[j][i]-bm[j]; 
  }
  
  // rotate pos by angle phi about (0,0)
  rotate2dim(pos[0],pos[1],nBead,phi); 
  
  // translate beadpos back to original centroid and save
  for (i=0;i<nBead;i++) {
    for (j=0;j<2;j++) cm_temp[i][j]=pos[j][i]+bm[j];
    moveBeadPos(i,cm_temp[i]);
  }
  return; 
}

/************************************************************************/
void bead::setScanDir(int iDiam,int dir,int& rStart,int& cStart,int& rEnd,
		  int& cEnd,int& rDelta, int& cDelta)
{ /* Set scanning direction. */ 
  
  if (dir == 0) { // left to right, top to bottom
    rStart = 0;
    rEnd = rows - iDiam;
    cStart = 0;
    cEnd = columns - iDiam;
    rDelta = 1;
    cDelta = 1;
  }
  else if (dir == 1) { // right to left, top to bottom
    rStart = 0;
    rEnd = rows - iDiam;
    cStart = columns-1;
    cEnd = iDiam - 1;
    rDelta = 1;
    cDelta = -1;
  }
  else if (dir == 2) { // left to right, bottom to top
    rStart = rows - 1;
    rEnd = iDiam - 1;
    cStart = 0;
    cEnd = columns - iDiam;
    rDelta = -1;
    cDelta = 1;
  }
  else if (dir == 3) { // right to left, bottom to top
    rStart = rows - 1;
    rEnd = iDiam - 1;
    cStart = columns - 1;
    cEnd = iDiam - 1;
    rDelta = -1;
    cDelta = -1;
  }
}

/************************************************************************/
int bead::withinRadius(int beadNum, double radius,
		       vector<pair<double,int> >& closeBeads)
{
  /* This function checks for beads that are within a radius of the
     bead with the id beadNum. This function returns a vector of the
     pairs of the distance to the bead from bead beadNum followed by
     the bead number for all of the beads within the radius. The
     vector is sorted so closest bead is first in the vector.
   */
  
  multimap<int, int>::iterator it;
  double r_sq = radius*radius;
  double dist_sq, dx, dy, cmr = cm[beadNum][0], cmc = cm[beadNum][1];
  int i, j, cell, iCell=(int)(nCell*(int)(cmr/dRow)+(int)(cmc/dCol));
  pair<multimap<int, int>::iterator, multimap<int, int>::iterator> range;
  
  for (i = -nCell; i <= nCell; i+=nCell) {
    for (j = -1; j <= 1; j++) {
      cell = iCell + i + j;
      range = cell_bead.equal_range(cell);
      for (it = range.first; it != range.second; it++) {
	if (it->second == beadNum) continue;
	dx = cm[it->second][0] - cmr;
	dy = cm[it->second][1] - cmc;
	dist_sq = dx*dx + dy*dy;
	if (dist_sq < r_sq) {
	  closeBeads.push_back(pair<double, int>(dist_sq, it->second));
	}
      }
    }
  }
  if (closeBeads.size() > 1) {
    sort(closeBeads.begin(), closeBeads.end());
  }
  return closeBeads.size();
}

/************************************************************************/
void bead::writeCor(string fname, string outMode)
{ /* Write coor file for beads. */ 
  
  FILE *ff;
  string corname,segname;
  
  double x1, y1, z1, massMin, massMax, occupancy;
  int imin, imax;
  multimap<int,int>::iterator it;
  
  // determine filename extension
  corname = addFileNameExt(fname,".cor");
  cout<<"  Write bead coor to "<<corname<<" outMode= " <<outMode<<endl;
  ff=fopen(corname.c_str(),"w");

  // find min/max mass: normalized mass entered as occupancy (col 55-60)
  MinMax(nBead, beadMass, &imin, &imax, &massMin, &massMax);

  int i,k;
  // use extended format for more than 99999 atoms
  if (nBead>99999) fprintf(ff,"%10d%s\n",nBead,"  EXT");
  else fprintf(ff,"%d\n",nBead);

  segname=tag.substr(0,4);  // 4-letter segname
  for (i=0;i<nBead;i++) {
    k=i+1;
    if (outMode=="verbose") {
      x1=cm[i][0]; y1=cm[i][1]; z1=zCoord;}
    else { // outMode=="image"
      y1=((double)rows-cm[i][0]); x1=cm[i][1]; z1=zCoord;} 
    occupancy=beadMass[i]/massMax;

    /* From coorio.src:
       fm2='(2I5,1X,A4,1X,A4,3F10.5,1X,A4,1X,A4,F10.5)'
       I,IRES,RES(IRES),ATYPE(I),Xyz(1,I),xYz(2,I),xyZ(3,I),SID,RID,WMAIN(I)
    */
    if (nBead>99999) { // use extended format
      fprintf (ff,
	       "%10d%10d  %-8s  %-8s%20.10f%20.10f%20.10f  %-8s  %-8d%20.10f\n",
	       i+1,k,"BEAD","O",x1,y1,z1,segname.c_str(),k,occupancy);
    }
    else {
    //    fprintf (ff,"%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %4s %-4d%10.5f\n",
      fprintf (ff,"%5d%5d %-4s %-4s%9.4f %9.4f %9.4f  %4s %-4d%10.5f\n",
	       i+1,k,"BEAD","O",x1,y1,z1,segname.c_str(),k,occupancy);
    }
  }
  fclose(ff);
}

/************************************************************************/
void bead::writeData(string fname, FILE *dat)
{  /* Writes bead data to a dat file with the name fname. The dat file 
      contains all of the information necessary to construct an identical 
      bead object.   */

  int iBead;
  bool newFile = (dat == NULL);
  string datname;

  datname = addFileNameExt(fname,".dat");

  if (newFile) {
    dat = fopen(datname.c_str(), "w");
    cout<<"  Write bead data to "<<datname<<endl;
  }

  fprintf(dat, "img_tag= %s\ntag= %s\nmaxBead= %d\nnBeads= %d\nbeadDiam= %d\nrows= %d\ncolumns= %d\ndRow= %.2f\ndCol= %.2f\n",
	  img_tag.c_str(), tag.c_str(), maxBead, nBead, beadDiam, rows, columns, dRow, dCol);

  fprintf(dat, "Bead     Mass               X                  Y                  Z\n");

  for (iBead = 0; iBead < nBead; iBead++) {
    fprintf(dat, "%-8d %-15.12e %-15.12e %-15.12e %-15.12e \n",
	    iBead, beadMass[iBead], /*cm[iBead][1], (double)rows-cm[iBead][0],*/
	    cm[iBead][0], cm[iBead][1], zCoord);
  }

  if (newFile) fclose(dat);
  return;
}

/************************************************************************/
void bead::writeImage(int beadL, int nColor, double *colorBg,
		      double *colorBead, string fname,string fmt)
{ /* Write bead coordinates to an image file.
     beadL: Sizes of beads. If <=0, do not draw. 
     nColor=1: grayscale, nColor=3: rgb
     colorB{bg,bead}: corresponding colors.
     fname,fmt: output file name.  */
  
  int i,flag;
  string ofname;

  Image img0(Geometry(0,0),"black"); 
  list<Magick::Drawable> drawBead; //construct drawing list
  
  ///////////////////////////////////////////////////////
  // Preparation
  ///////////////////////////////////////////////////////

  if (fmt=="") { // if no format given, check if fname already has an extension
    flag=checkImageNameExt(fname);
    ofname=(flag==0)?(addFileNameExt(fname,"tif")):fname;
  }
  else ofname=addFileNameExt(fname,fmt); //set image file name

  // check sizes
  if (nBead==0) {
    cout<<"  "<<"WARNING: nBead= 0. No image written."<<endl;
    return;
  }
  if (beadL<=0) {
    cout<<"  "<<"WARNING: Zero bead size. No image written"<<endl;
    return;
  }

  ///////////////////////////////////////////////////////
  // Write parameters and render
  ///////////////////////////////////////////////////////

  cout<<"  Write image to "<<ofname<<endl;
  cout<<"  Bead size: "<<setw(5)<<beadL<<endl;
  if (nColor==3) { // rgb
    cout<<endl<<"    RGB colors:"<<endl;
    cout<<"    Background: ";
    for (i=0;i<3;i++) cout<<setw(5)<<setprecision(2)<<colorBg[i]<<" ";
    cout<<endl;
    cout<<"    Bead: ";
    for (i=0;i<3;i++) cout<<setw(5)<<setprecision(2)<<colorBead[i]<<" ";
    cout<<endl;
    img0.extent(Geometry(columns,rows),ColorRGB(colorBg[0],colorBg[1],colorBg[2]));
  } //   if (nColor==3) { // rgb
  else { // gray: only use the first color
    cout<<endl<<"    Grayscale colors:"<<endl;
    cout<<"    Background: "<<setw(5)<<setprecision(2)<< *colorBg <<endl;
    cout<<"    Bead: "<<setw(5)<<setprecision(2)<< *colorBead <<endl;

    img0.extent(Geometry(columns,rows), ColorGray(*colorBg)); //base image
  } //   else { // gray: only use the first color

  img0.modifyImage(); //ensure there is only 1 reference to underlying image
  writeImageBead(drawBead,cm,nBead,beadL,nColor,colorBead);
  img0.draw(drawBead); //draw all objects on list
  img0.depth(8); //write to 8-bit image
  img0.write(ofname);  

  return; 
}

/************************************************************************/
void bead::writePsf(string fname)
{ /* Write psf file for beads. */ 
  int j,resid=1;
  FILE *psf;
  string segname,psfname;
  multimap<int,int>::iterator it;

  psfname = addFileNameExt(fname,".psf");  // determine filename extension
  cout<<"  Write bead psf to "<<psfname<<endl;
  
  psf=fopen(psfname.c_str(),"w");
  
  fprintf(psf,"PSF \n\n !NTITLE\n");
  fprintf(psf,"* PSF for %s \n* \n\n",fname.c_str());
  
  fprintf(psf,"%6d !NATOM\n",nBead);
  
  segname=tag.substr(0,4);  // 4-letter segname
  for (j=0;j<nBead;j++) { 
    fprintf(psf,"%8d %4s %04d BEAD O       0 %14.6f %14.6f %d\n",
	    j+1,segname.c_str(),resid,0.0,beadMass[j],0);
    /* from psfres.src: fmt01=(I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8)
       I(8),SEGID(ISEG)(4),RESID(IRES)(4),RES(IRES)(4),
       ATYPE(I)(4),IAC(I)(4; atomID),CG(I)(14.6;charge),
       AMASS(I)(14.6;mass),IMOVE(I)(8)
    */
  }
  fprintf(psf,"\n%8d !NBOND\n",0);
  fclose(psf);
  return;
}

/************************************************************************/
// Non-member functions
/************************************************************************/

void checkDuplicateBeadTag(vector<bead> * vbead,string tag)
{ // In vbead, find iterator for the element with tag
  vector<bead>::iterator it_bead;
  for (it_bead=(*vbead).begin();it_bead!=(*vbead).end();++it_bead) {
    if ( (*it_bead).tag==tag) error_dodri("Duplicate bead tag: " + tag);
  }
  return;
}

/************************************************************************/
void copyBead(const bead &b0, double **cm1, double *beadMass1,string tag1)
{ /* Copy b0 to new bead with tag1, with center of mass values cm. 
     If tag1 exists, that image is overwritten. If bead with tag1 exists,
     size of b0.cm[nBead][2] and size of cm1[nBead][2] and 
     beadMass1[nBead] must match (no check made). 

     Note: Cannot call copyBead more than once within a bead member 
     function since vbead changes. 
  */

  int j,iBead,iCell;
  multimap<int,int> cell_bead_tmp, bead_cell_tmp;
  vector<bead>::iterator it_bead=findBead(&vbead,tag1);
  cout<< "  Bead "<<b0.tag<<" being copied to "<<tag1<<endl;

  if (it_bead !=vbead.end()){ 
    cout<<"  WARNING: Bead " + tag1 + " is overwritten."<<endl;
    // The following operations are the same as the constructor
    // bead::bead(bead *b0, double **cm1, string tag1)
    it_bead->img_tag=b0.img_tag;
    if (it_bead->nBead!=b0.nBead)
      error_dodri("copyBead: Bead count of the two beads does not match.");
    it_bead->maxBead=b0.maxBead;
    it_bead->nBead=b0.nBead;
    it_bead->beadDiam=b0.beadDiam;
    it_bead->zCoord=b0.zCoord;
    it_bead->rows=b0.rows; it_bead->columns=b0.columns;
    it_bead->dRow=b0.dRow; it_bead->dCol=b0.dCol;
    for (iBead=0;iBead<b0.nBead;iBead++) {
      for (j=0;j<2;j++) it_bead->cm[iBead][j]=cm1[iBead][j];
      it_bead->beadMass[iBead]=beadMass1[iBead];
    } //     for (iBead=0;iBead<nBead;iBead++) {
    //re-calculate cell number in case of major movement
    for (iBead=0;iBead<b0.nBead;iBead++) {
      iCell=nCell*(int)(cm1[iBead][0]/(b0.dRow))+
	(int)(cm1[iBead][1]/(b0.dCol));
      cell_bead_tmp.insert(pair<int,int>(iCell,iBead));
      bead_cell_tmp.insert(pair<int,int>(iBead,iCell));
      it_bead->cell_bead=cell_bead_tmp;
      it_bead->bead_cell=bead_cell_tmp;
    } //     for (iBead=0;iBead<nBead;iBead++) {
  } //   if (it_bead !=vbead.end()) {
  else {  //new bead
    vbead.push_back(bead(b0,cm1,beadMass1,tag1));
  }
  return; 
}

/************************************************************************/
vector<bead>::iterator findBead(vector<bead> * vbead,string tag)
{ // In vbead, find iterator for the element with tag
  vector<bead>::iterator it_bead;
  for (it_bead=(*vbead).begin();it_bead!=(*vbead).end();++it_bead) {
    if ( (*it_bead).tag==tag) break;
  }
  return it_bead;
}
