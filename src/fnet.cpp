#include "dodri.h"

/************************************************************************/
//Constructors
/************************************************************************/

fnet::fnet(bead *bd, string tag0) 
{
  int ibead;
  multimap<int,int>::iterator it;

  // variables inherited from bead
  bead_tag=bd->tag; 
  img_tag=bd->img_tag;
  tag=tag0; // tag of this fnet
  nBead = bd->nBead ; beadDiam= bd->beadDiam;
  rows=bd->rows; columns=bd->columns;
  dRow=bd->dRow; dCol=bd->dCol;
  zCoord=bd->zCoord; 
  maxBead=bd->maxBead;
  maxFilLength=maxBead;

  cell_bead.clear(); bead_cell.clear(); 
  for (it=(bd->bead_cell).begin(); it!=(bd->bead_cell).end(); it++)  
    bead_cell.insert(pair<int,int>(it->first,it->second));

  for (it=(bd->cell_bead).begin(); it!=(bd->cell_bead).end(); it++) 
    (this->cell_bead).insert(pair<int,int>(it->first,it->second));

  cm=new double*[maxBead];  
  beadMass=new double[maxBead];
  for(ibead=0;ibead<(int)maxBead;ibead++) cm[ibead]=new double[2];
  for (ibead=0;ibead<nBead;ibead++) {
    cm[ibead][0]=bd->cm[ibead][0]; cm[ibead][1]=bd->cm[ibead][1];
    beadMass[ibead]=bd->beadMass[ibead];
  }

  // fnet variables
  bond.clear();  bead_fil.clear(); fil_bead.clear();
  nBond=nFilament=0;

  cout<<"  fnet_tag= "<< tag << " initialized."<<endl;
}

/************************************************************************/
fnet::fnet(img *img0, int iDiam,unsigned char pxlCut,string tag0)
{ /* Construct fnet directly from img0.
 
     For each cluster, locate a boundary pixel. Travel CCW, marking 
     boundary pixels. Track distance between marked boundary pixels, 
     add a bead and bond (assigned simultaneously) between those two 
     pixels if distance cutoff is met. If search does not return to 
     starting position, make a CW search from starting pos. 
     NOTE: Beads created here do not have the corresponding bead object.
     Intensity checks currently assume input image is binary. 

     iDiam: feature scale
  */
  
  int ncluster=0,nmax_clstr,idum,idum1,flag,ipBd,ipBd0;
  int i,j,k,ix=-1,iy=-1,ix1=-1,iy1=-1,dir0=0,dir;
  double x0,y0,x1=-1,y1=-1,xch,ych,rdum[2],rdum1[2],bondL;
  pair<double,double> pdum; 
  vector<pair<double,double>> *edge_pos,corner; 
  // edge_pos: list of pixel CORNERS shared with different intensity pixel (ix,iy)
  
  vector<vector<pair<int,int>>> cluster;
  vector<pair<double,double>>::iterator it,jt;
  vector<pair<int,int>>:: iterator iit; 
  
  // set tag and img global vars
  bead_tag=""; // no corresponding bead
  img_tag=img0->tag;
  tag=tag0; // tag of this fnet
  nBead = 0;
  beadDiam=iDiam; 
  rows=img0->rows; columns=img0->columns;
  dRow=img0->dRow; dCol=img0->dCol;
  zCoord=0.; // set to zero for now.
  
  ncluster=img0->findCluster(pxlCut, cluster);
  edge_pos=new vector<pair<double,double>>[ncluster];   
  
  nmax_clstr=0; // max cluster size  
  for (i=0;i<ncluster;i++) {
    edge_pos[i].clear();     
    k=cluster[i].size(); 
    nmax_clstr=(k>nmax_clstr)?k:nmax_clstr;
    for (j=0;j<k;j++) {
      ix=cluster[i][j].first; iy=cluster[i][j].second; // pixel to be checked
      assert((int)img0->pxl_gray[ix][iy]!=0); 
      getPxlEdge(ix,iy,pxlCut,edge_pos[i],img0);
    }
  }
  
  // remove duplicate entries within edge_pos[i]
  for (i=0;i<ncluster;i++) {
    if (edge_pos[i].empty())
      error_dodri(" No edges identified in cluster."); 
    corner.clear(); 
    corner=edge_pos[i]; // temporary copy
    sort(corner.begin(),corner.end());
    corner.erase( unique( corner.begin(), corner.end() ), corner.end() );
    edge_pos[i].clear();
    edge_pos[i]=corner; 
  }

  // initialize bead, fnet objects
  maxBead=nmax_clstr*ncluster; // maxBead = (max cluster size)\times(ncluster)
  maxFilLength=maxBead;
  cm=new double*[maxBead];  
  beadMass=new double[maxBead];
  for(i=0;i<(int)maxBead;i++) cm[i]=new double[2];
  
  cell_bead.clear(); bead_cell.clear(); 
  bond.clear();  bead_fil.clear(); fil_bead.clear();
  nBond=nFilament=0;
  bondL=(double)beadDiam;  // scale bond length by beadDiam
  
  /*** Begin search ***/
  for (i=0;i<ncluster;i++) {
    if ((int)edge_pos[i].size()<=4*beadDiam) continue;
    // at min, if cluster size=1, and beadDiam==1, then 4 edges are needed  
  newSeed: 
    ++nFilament; 
    x0=edge_pos[i][0].first; y0=edge_pos[i][0].second; // start search, store seed
    edge_pos[i].erase(edge_pos[i].begin());

    // add a new bead
    bead::createBead(1.,x0,y0);
    beadMass[nBead-1]=1.0;     
    putBeadToFilament(nBead-1,nFilament);
    ipBd0=ipBd=nBead-1; // ipBd0= save nBead of seed 
    
    xch=x0; ych=y0;  // set search start 
    dir=dir0;    flag=0;

    while (!flag) {
      // check presumptive edges
      corner.clear(); 
      listCorners(xch,ych,dir,1.,corner);  // dir=0, dr=1; list corners dr from (xch,ych)
      for (it=corner.begin();it!=corner.end();it++) {
	x1=(*it).first; y1=(*it).second;
	pdum=make_pair(x1,y1); 
	if (find(edge_pos[i].begin(),edge_pos[i].end(),pdum)!=edge_pos[i].end()) { // cluster edge
	  // check if this potential edge b/en two pixels is a TRUE boundary edge
	  // need to locate ix,iy of two pixels to compare
	  if (x1==xch) { // same x, comp pixels top/bottom
	    ix=(int)(x1-0.5); iy=(int)((y1+ych)/2); // (y= midpoint)
	    ix1=(int)(x1+0.5); iy1=(int)((y1+ych)/2); // position of the two pixels
	  }
	  else if (y1==ych) { // same y, comp pixels left/right
	    ix=(int)((x1+xch)/2); iy=(int)(y1-0.5);
	    ix1=(int)((x1+xch)/2); iy1=(int)(y1+0.5); 	    
	  }
	  // ix,iy and ix1,iy1 may be SLIGHTLY out of boundary (0.5 units). Adjust.
	  // since a "stop" if search goes out of boundary is done below, this is fine.
	  ix=(ix<0)?(0):(ix);	  ix=(ix>=rows)?(rows-1):(ix);
	  ix1=(ix1<0)?(0):(ix1);  ix1=(ix1>=rows)?(rows-1):(ix1);	  
	  iy=(iy<0)?(0):(iy);	  iy=(iy>=rows)?(rows-1):(iy);
	  iy1=(iy1<0)?(0):(iy1);  iy1=(iy1>=rows)?(rows-1):(iy1);	  
	  
	  idum=(int)(img0->pxl_gray[ix][iy]);
	  idum1=(int)(img0->pxl_gray[ix1][iy1]);
	  if (idum!=idum1) { // true edge found
	    // remove this entry from edge_pos, update tracking position
	    jt=find(edge_pos[i].begin(),edge_pos[i].end(),pdum);
	    edge_pos[i].erase(jt);	    
	    xch=x1; ych=y1;  // position for next search
	    
	    // check if distance cutoff met to add a bead
	    rdum[0]=cm[ipBd][0]; rdum[1]=cm[ipBd][1]; // pos of previously added bead
	    rdum1[0]=x1; rdum1[1]=y1;  
	    if (getDistPos(rdum,rdum1,0)>bondL) {
	      bead::createBead(1.,x1,y1); // nBead-1=index of this bead
	      beadMass[nBead-1]=1.0;     	      
	      putBeadToFilament(nBead-1,nFilament);
	      registerOneBond(ipBd,(nBead-1),-1,-1,img0);
	      ipBd=nBead-1; // save index of last added bead
	    }
	    break; 
	  }
	} 
      } //       for (it=corner.begin();it!=corner.end();it++) {
      
      /* Check ending conditions for while loop  */
      if ((int)edge_pos[i].size()<=beadDiam) flag=1; 
      if (!(xch==x1 && ych==y1)) {       // a new edge not found      
	if (dir==dir0) { // try next search dir
	  dir=(dir==0)?(1):(0); 
	  xch=x0; ych=y0; // use original seed   // goes to start of while
	  ipBd=ipBd0;
	}
	else {
	  if (!edge_pos[i].empty()) goto newSeed; 
	  else flag=1; // end while loop
	}
      }
      if (xch<0 || xch>=rows || ych<0 || ych>=columns) {  // ends of image reached
	if (dir==dir0) { // try next search dir
	  dir=(dir==0)?(1):(0); 
	  xch=x0; ych=y0; // use original seed    // goes to start of while
	  ipBd=ipBd0; // switch direction 
	}
	else goto newSeed; // look for new seed
      }
    } //     while (!flag) {
  } //   for (i=0;i<ncluster;i++) {
  
  // info on generated BOC
  cout<<"         fnet_tag= "<<tag<<" created: "<<endl;
  cout<<"           "<<nBead<<" beads, "<< nBond<<" bonds, "<<
    nFilament<<" filaments assigned."<<endl;    
  
  return;
}

/************************************************************************/
fnet::fnet(string tag,string ifname, string bead_tag0,ifstream *fin)
{ /* Construct fnet from the data in a .dat file with name ifname. 
     Does not build bead. Bead should be built/loaded separately.   */
  this->tag = tag;
  if (findFnet(&vfnet,tag)==vfnet.end())
    error_dodri("Fnet tag not initialized."); 
  this->readData(ifname,bead_tag0,fin);
}

/************************************************************************/
// Member functions
/************************************************************************/

void fnet::addBead2Fn(double mass,double cmr,double cmc)
{ /* Copy of bead::createBead() for fnet edits of bead 
     when bead is not built. */
  
  cm[nBead][0]=cmr; cm[nBead][1]=cmc;
  beadMass[nBead]=mass;
  int iCell=(int)(nCell*(int)(cmr/dRow)+(int)(cmc/dCol));
  cell_bead.insert( pair<int,int>(iCell,nBead) );
  bead_cell.insert( pair<int,int>(nBead,iCell) );
  nBead++;   
}

/************************************************************************/
void fnet::assignFilSeqPos(int filSeq[], double x[], double y[], int N,
			   signed char dir)
{ /* For bead indexes stored in filSeq, get x,y positions for 
     n beads, then reverse if dir<0 */
  
  int j,k;
  for (j=0;j<N;j++) {
    if (dir>=0) k=j; else k=N-1-j;
    x[j]=cm[filSeq[k]][0]; y[j]=cm[filSeq[k]][1];
  }
}

/************************************************************************/
double fnet::contourLength(int *filSeq, int ipos, int jpos) 
{ /* Get contour length from ipos to jpos (inclusive) in filSeq. */
  int i;  double s=0., rdum;
  for (i=ipos;i<jpos;i++) {
    rdum=getDist(filSeq[i],filSeq[i+1],1); s+=rdum;
  }
  return s;
}

/************************************************************************/
void fnet::deleteBead(int iBead, char sort)
{ /* Delete iBead. Procedure:
     1) Delete the cell association
     2) Delete it from the filament. If the bead is in the middle, remove bonds
     to two neighboring beads and connect the neighbors.
     3) Delete any remaining bond.
     4) if sort==0; stop at 3. (nBead doesn't decrease)
     if sort!=0: Change the last bead (index nBead) to index iBead;  
     Update its cell association, filament, other bonds, and **cm array and *beadMass
  */
  
  multimap<int,int>::iterator it,jt,jt1;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> it_range,
    jt_range;

  int i,j,k,filLength, filLength0, iCell,jCell, iFil,jFil, jBead,kBead;
  int *filSeq=new int[maxFilLength]; // sequential bead index
  int  *bondVecList=new int[maxVec],nVec; 

  // delete cell association
  it=bead_cell.find(iBead); iCell=it->second;
  it_range=cell_bead.equal_range(iCell);
  for (it=it_range.first;it!=it_range.second;it++)  {
    if (it->second == iBead) { cell_bead.erase(it); break; }
  }
  bead_cell.erase(iBead);

  // delete from filament 
  it=bead_fil.find(iBead);
  if (it!=bead_fil.end()) { // iBead is a part of a filament
    iFil=it->second;
    filLength=findFilamentSeq(iFil,filSeq);
    filLength0=filLength; // track if cyclic 
    filLength=abs(filLength); 
    // i becomes index of iBead along iFil
    for (i=0;i<filLength;i++) if (filSeq[i]==iBead) break;
    jt_range=fil_bead.equal_range(iFil);
    for (jt=jt_range.first;jt!=jt_range.second;jt++) 
      if (jt->second == iBead) { fil_bead.erase(jt); break; }
    bead_fil.erase(it);
    // If iBead is in the middle of the filament, connect neighbor beads
    //  if (filLength0>0) { // non-cyclic fil
    if ((i!=0)&&(i!=(filLength-1))) {
      jBead=filSeq[i-1]; kBead=filSeq[i+1];
      jt_range=bond.equal_range(jBead);
      for (jt=jt_range.first;jt!=jt_range.second;jt++) 
	assert(jt->second != kBead); // ensure no bond bet jBead-kBead
      registerOneBond(jBead,kBead,-1,-1,&(*findImg(this->img_tag)));
    }
    else if (filLength0<-5) { // cyclic filament w/ at min 3 beads
      if (i==0) { // index is 0 ("first" bead in fil)
    	j=filLength-1; k=i+1; }
      else if (i==filLength-1) { // index is filLength-1 ("last" bead in fil)
    	j=i-1; k=0;
      }
      else goto nomatch; 
      jBead=filSeq[j]; kBead=filSeq[k];
      jt_range=bond.equal_range(jBead);
      for (jt=jt_range.first;jt!=jt_range.second;jt++) 
	assert(jt->second != kBead); // ensure no bond bet jBead-kBead
      registerOneBond(jBead,kBead,-1,-1,&(*findImg(this->img_tag)));
    }
  nomatch: 
    // for a single-bead filament, move nFilament to iFil after deletion
    if ((fil_bead.count(iFil)==0)&&(iFil<nFilament)) {
      jt_range=fil_bead.equal_range(nFilament);
      for (jt=jt_range.first;jt!=jt_range.second;jt++) {
	jBead=jt->second;
	bead_fil.erase(jBead); // erase association to nFilament
	putBeadToFilament(jBead,iFil);
      }
      fil_bead.erase(nFilament);
      --nFilament;
    } // if ((fil_bead.count(iFil)==0)&&(iFil<nFilament)) {
  }// if (it!=bead_fil.end()) { // iBead is a part of a filament
  // delete bond
  nVec=getBondVec(iBead,-1,bondVecList);
  for (i=0;i<nVec;i++) {
    kBead=bondVecList[i];
    removeOneBond(iBead,kBead);
  }

  /* Re-index nBead-1 to iBead */
  if (sort==0) return; // nBead doesn't change
  jBead=nBead-1;
  if (jBead==iBead) {goto dB0;} // Skip if iBead is the last bead

  // Cell index
  it=bead_cell.find(jBead);  jCell=it->second;
  it_range=cell_bead.equal_range(jCell);
  for (it=it_range.first;it!=it_range.second;it++) 
    if (it->second == jBead) { cell_bead.erase(it); break; }
  bead_cell.erase(jBead);
  bead_cell.insert(pair<int,int>(iBead,jCell));
  cell_bead.insert(pair<int,int>(jCell,iBead));
  // filament index
  it=bead_fil.find(jBead);
  if (it!=bead_fil.end()) { // jBead is a part of a filament
    jFil=it->second;
    jt_range=fil_bead.equal_range(jFil);
    for (jt=jt_range.first;jt!=jt_range.second;jt++) 
      if (jt->second == jBead) {fil_bead.erase(jt); break;}
    bead_fil.erase(jBead);
    bead_fil.insert(pair<int,int>(iBead,jFil));
    fil_bead.insert(pair<int,int>(jFil,iBead));
  }
  else  { // save unbound bead to filament=0
    fil_bead.insert(pair<int,int>(0,iBead));
    bead_fil.insert(pair<int,int>(iBead,0));
  }
    
  // bond association
  nVec=getBondVec(jBead,-1,bondVecList);
  for (i=0;i<nVec;i++) {
    kBead=bondVecList[i];
    removeOneBond(jBead,kBead);
    registerOneBond(iBead,kBead,-1,-1,&(*findImg(this->img_tag)));
  }
  
  // Change **cm array, *beadMass
  cm[iBead][0]=cm[jBead][0];cm[iBead][1]=cm[jBead][1];
  beadMass[iBead]=beadMass[jBead];
  
 dB0: --nBead; 
  // free up memory
  delete [] filSeq;  delete [] bondVecList;
  return;
}

/************************************************************************/
char fnet::deleteMark(int p, int q)
{ /* Mark bonds between p and q for deletion.
     returns 0: if already marked. 1: for a new mark */ 
  multimap<int,int>::iterator it;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> it_range;
  it_range=bond_del.equal_range(p);
  for (it=it_range.first;it!=it_range.second;it++) {
    if (it->second==q) return 0;
  }
  bond_del.insert( pair<int,int>(p,q) );
  return 1;
}

/************************************************************************/
void fnet::equalizeFilamentBond(double rCut)
{ /* Make beads equally spaced in a filament. Bond length becomes approximately
     rCut*beadDiam.

  1) First get rid of narrowly spaced beads: Go over two consecutive beads
    and if the sum of their bond lengths is less than rCut1=rCut*beadDiam,
    delete the middle bead.
  2) Measure avg bond length bLeng
  3) Go over each bond: if the length is longer than 1.5*bLeng, add additional
    beads.
  4) Go over two consecutive bonds,  apply fitQuad0 to move the midle
   bead in the mid-position between 1st and 3rd bead.

  If necessary, this command can be called multiple times to enhance equality
  of bond lengths. */
  
  int iFil, filLength, *filSeq=new int[maxFilLength];
  int b1,b2,b3,i,p;
  double contLeng,mass,bLeng;
  double *x=new double[maxFilLength], *y=new double[maxFilLength];
  double *uVec=new double[2],*rVec1=new double[2],*rVec2=new double[2];
  double rCut1=rCut*beadDiam,rdum;
  double *vel=new double[2], *acc=new double[2];
  double *pVec1=new double[2], *pVec2=new double[2], *pVec3=new double[2];

  multimap<int,int>::iterator it;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> it_range;

  cout<<"  equalize_bond: rCut= "<<rCut<<endl;

  for (iFil=1;iFil<=nFilament;iFil++) {
    filLength=findFilamentSeq(iFil,filSeq);
    filLength=abs(filLength); 
    if (filLength>2) { // step 1)
      for (p=0;p<(filLength-2);p++){
	contLeng=contourLength(filSeq,p,p+2);
	if (contLeng<rCut1) { // delete middle bead
	  deleteBead(filSeq[p+1],1);
	  if (p<(filLength-3)) {
	    filLength=findFilamentSeq(iFil,filSeq);
	    filLength=abs(filLength); 	    
	    --p; // check p again
	  }
	}
      } // for (p=0;p<(filLength-2);p++){
    } // if (filLength>2) { // step 1)

    filLength=findFilamentSeq(iFil,filSeq); // step 2)
    filLength=abs(filLength);
    
    contLeng=contourLength(filSeq,0,(filLength-1));
    bLeng=contLeng/(double)(filLength-1);
    for (p=0;p<(filLength-1);p++) { // step 3)
      contLeng=contourLength(filSeq,p,p+1);
      if (contLeng>(1.5*bLeng)) {
	b1=filSeq[p]; b2=filSeq[p+1]; 
	for (i=0;i<2;i++) {
	  uVec[i]=cm[b2][i]-cm[b1][i];
	  uVec[i]*=(bLeng/contLeng);
	}
	removeOneBond(b1,b2); // remove old bond
	// give avg mass to new beads
	mass=0.5*(beadMass[b2]+beadMass[b1]);
	for (i=0;i<2;i++) pVec1[i]=cm[b1][i]+uVec[i];
	bead::createBead(mass,pVec1[0],pVec1[1]);
	putBeadToFilament(nBead-1,iFil);
	registerOneBond(b1,(nBead-1),-1,-1,&(*findImg(this->img_tag)));
	registerOneBond(b2,(nBead-1),-1,-1,&(*findImg(this->img_tag)));
	filLength=findFilamentSeq(iFil,filSeq); // refresh
	filLength=abs(filLength); 
      } // if (contLeng>(1.5*bLeng)) {
    } // for (p=0;p<(filLength-1);p++) { // step 3)
    filLength=findFilamentSeq(iFil,filSeq);
    filLength=abs(filLength); 
    for (p=0;p<(filLength-2);p++){ // step 4)
      for (i=0;i<2;i++) {
	b1=filSeq[p];	 pVec1[i]=cm[filSeq[p]][i]; 
	b2=filSeq[p+1]; pVec2[i]=cm[b2][i]; 
	b3=filSeq[p+2]; pVec3[i]=cm[b3][i];
      }
      fitQuad0(pVec1,pVec2,pVec3,vel,acc,2); // step 2b)
      bLeng=0.5*contourLength(filSeq,p,p+2);
      for (i=0;i<2;i++) {
	rdum=0.5*acc[i]*bLeng*bLeng;
	pVec2[i]=pVec1[i]+vel[i]*bLeng+rdum;
      }
      bead::moveBeadPos(b2,pVec2); 
    } // for (p=0;p<(filLength-2);p++){ // step 4)
  } //   for (iFil=1;iFil<=nFilament;iFil++) {

  cout<<"  Bead count: "<<nBead<<"; Bond count: "<<nBond<<endl;

  delete [] filSeq; delete [] x; delete [] y;
  delete [] uVec; delete [] rVec1; delete [] rVec2;
  delete [] vel; delete [] acc;
  delete [] pVec1; delete [] pVec2; delete [] pVec3;  
  return;
}

/************************************************************************/
int fnet::findCluster(multimap<int,int> *bead_cluster, 
			   multimap<int,int> *cluster_bead)
{ /* Find clusters (groups of bonded beads that are not part of any 
     filament).
     Return value: number of clusters. 
     Cluster number ranges [1,nCluster] (exclude 0)  */
  
  int i,nCluster=0,iBead,jBead,kBead, clusterSize;
  int nIter, maxIter=nBead;
  multimap<int,int>::iterator it,jt,kt;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> it_range,
    jt_range;

  for (iBead=0;iBead<nBead;iBead++){
    // skip if iBead is part of a filament or a cluster
    if (bead_fil.find(iBead)!=bead_fil.end()) continue; 
    if ((*bead_cluster).find(iBead)!=(*bead_cluster).end()) continue;
    
    it_range=bond.equal_range(iBead);
    i=(int)bond.count(iBead);
    if (i<1) continue; // isolated bead

    clusterSize=1; ++nCluster; // new cluster
    (*bead_cluster).insert(pair<int,int>(iBead,nCluster));
    (*cluster_bead).insert(pair<int,int>(nCluster,iBead));
    nIter=0;
  fC0: it_range=(*cluster_bead).equal_range(nCluster);
    for (it=it_range.first;it!=it_range.second;it++) {
      jBead=it->second;  jt_range=bond.equal_range(jBead);
      for (jt=jt_range.first;jt!=jt_range.second;jt++) {
	kBead=jt->second;
	if (bead_fil.find(kBead)!=bead_fil.end()) continue; 
	kt=(*bead_cluster).find(kBead);
	if (kt!=(*bead_cluster).end()) {
	  assert(kt->second==nCluster); 
	  continue; 
	}
	(*bead_cluster).insert(pair<int,int>(kBead,nCluster));
	(*cluster_bead).insert(pair<int,int>(nCluster,kBead));
	++clusterSize;
      }
    } // for (it=it_range.first;it!=it_range.second;it++) {
    // Check if there are leftover beads
    ++nIter;
    if (nIter==maxIter)
      error_dodri("findCluster: max number of iterations reached");
    it_range= (*cluster_bead).equal_range(nCluster);
    for (it=it_range.first;it!=it_range.second;it++) {
      jBead=it->second;  jt_range=bond.equal_range(jBead);
      for (jt=jt_range.first;jt!=jt_range.second;jt++) {
	kBead=jt->second;
	if (bead_fil.find(kBead)!=bead_fil.end()) continue; 
	kt=(*bead_cluster).find(kBead);
	if (kt!=(*bead_cluster).end()) {
	  assert(kt->second==nCluster); 
	  continue; 
	}
	goto fC0; // something left over
      }
    } // for (it=it_range.first;it!=it_range.second;it++) {
  } // for (iBead=0;iBead<nBead;iBead++){
  return nCluster;
}

/************************************************************************/
void fnet::findFilament(int flag)
{ /* Idenfity filaments only in the unbranched region. 
     Exclude single-bead filament. Single-bead filaments are fil=0.
     If flag==1 (default): delete bonds that would make filaments cyclic.  */

  int i,j, j0, nVec;
  multimap<int,int>::iterator it, it1;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> 
    it_range, it1_range;
  
  for (i=0;i<nBead;i++) {
    if (bead_fil.find(i)!=bead_fil.end()) continue; // already registered
    nVec=(int)bond.count(i);
    if (nVec==0) { //unbound bead
      fil_bead.insert(pair<int,int>(0,i)); 
      bead_fil.insert(pair<int,int>(i,0));
    }
    if (nVec==1) { // single-bond filament
      it=bond.find(i); j=it->second; // j: bonded bead
      if ((int)bond.count(j)!=1) continue; // j connected to another bead
      else { // i-j isolated
	if (bead_fil.find(j)!=bead_fil.end()) continue; // already registered
	++nFilament; 
	putBeadToFilament(i,nFilament); putBeadToFilament(j,nFilament);
	continue;
      }  // i-j isolated
    } // if (nVec==1) { // single-bond filament
    if (nVec!=2) continue; // only consider middle of a filament
    it_range=bond.equal_range(i); 
    // it1: last elment in it_range, j0,j: beads connecting to i
    it1=it_range.second; --it1;
    j0=(int)it_range.first->second; j=(int)it1->second; 
    if (((int)bond.count(j0)>2)&&((int)bond.count(j)>2))
      continue; // connecting beads do not continue as 2-bond 
    // one end belongs to another fil
    if (bead_fil.find(j)!=bead_fil.end()) { 
      if (((int)bond.count(j0)>2)||(bead_fil.find(j0)!=bead_fil.end()))
	continue; // other end belongs to another fil or not 2-bonded
    }
    // one end belongs to another fil
    if (bead_fil.find(j0)!=bead_fil.end()) { 
      if (((int)bond.count(j)>2)||(bead_fil.find(j)!=bead_fil.end()))
	continue; // other end belongs to another fil or not 2-bonded
    }
    ++nFilament; // Test for single-bead filament passed. Set new filament.
    putBeadToFilament(i,nFilament);
    bond_del.clear();
    for (it=it_range.first;it!=it_range.second;it++) { // 1 or 2 bonds
      j0=i; // previous bead
      j=it->second; // next bead
      // Check if j belongs to a filament
      if (bead_fil.find(j)!=bead_fil.end()) { // Prevent loop formation
	deleteMark(j,j0); continue;
      }
      nVec=(int)bond.count(j);
      while (nVec==2) { // follow the chain and register into nFilament
        it1=bead_fil.find(j);
        if (it1!=bead_fil.end()) break; // j already in a filament
        putBeadToFilament(j,nFilament);
        it1_range=bond.equal_range(j); 
        for (it1=it1_range.first;it1!=it1_range.second;it1++) { 
          if (it1->second==j0) continue; // previous bead
          else {
            j0=j; // reset prev bead
            j=it1->second; //  next bead
            break; // break the loop after finding the next bead
          }
        }
        nVec=(int)bond.count(j);
      }
      if (nVec==1) putBeadToFilament(j,nFilament); // Filament end 
    } // for (it=it_range.first;it!=it_range.second;it++) { // 1 or 2 bonds
    if (flag) { // do not allow cyclic fils 
      for (it=bond_del.begin();it!=bond_del.end();it++) 
	removeOneBond(it->first,it->second);
    }
  } // for (i=0;i<nBead;i++) {
  cout<<"  findFilament:: Filaments found= "<<nFilament<<endl;
  return; 
}

/************************************************************************/
int fnet::findFilamentSeq(int i, int *filSeq)
{ /* For filament i, order it and store the sequence of bead numbers in
     array filSeq. Returns: number of beads in the filament.

    For cyclic filament, returns negative of the number of beads. */
  
  int filLength=fil_bead.count(i),nVec, p, q, j,k,p_end[2], n0;
  int cflag=0; // flag for cyclic filament
  multimap<int,int>::iterator it, it1, it2;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> 
    it_range, it1_range;

  if (filLength>maxFilLength)
    error_dodri("Filament longer than maxFilLength.");
  else if (filLength==1) { // single-bead filament
    filSeq[0]=fil_bead.find(i)->second;  return 1;
  }
  else { assert(filLength>0); } // nothing to do
  //n0=0;
  for (j=0;j<filLength;j++) filSeq[j]=-1;
  p_end[0]=p_end[1]=-1;    it_range=fil_bead.equal_range(i);
  // Find end points
  for (it=it_range.first;it!=it_range.second;it++) { 
    nVec=0; // number of bonds to the current bead w/in the filament
    p=it->second;
    it1_range=bond.equal_range(p); // bonds connected to this bead
    for (it1=it1_range.first;it1!=it1_range.second;it1++) {
      q=it1->second; it2=bead_fil.find(q);
      j=it2->second;  // skip bonds to other filaments
      if ((it2!=bead_fil.end())&&(j==i)) ++nVec; 
    }
    if (nVec==1) { // end point found
      if (p_end[0]==-1) p_end[0]=p;
      else {
        if (p_end[1]==-1) p_end[1]=p;
	else error_dodri("More than 2 end points in filament "+to_string(i));
      }
    }
    else if (nVec>2) {
      error_dodri("Filament "+to_string(i)+" has "
		 +to_string(nVec)+" branches.");
    }
    else if (nVec==0) { // single-bead filament
      filSeq[0]=p; return 1; // assert(j==0); 
    }
    else  continue; // nVec==2
  }
  if (p_end[0]!=-1) { // linear filament
    filSeq[0]=p_end[0];
  } // if (p_end[0]!=-1) { // linear filament
  else { // cyclic filament
    it=it_range.first; // first bead
    p=it->second; filSeq[0]=p;
    cflag=1;
  } //   else { // cyclic filament

  for (j=1;j<filLength;j++) { // Find sequence of bonds
    it1_range=bond.equal_range(filSeq[j-1]);
    for (it1=it1_range.first;it1!=it1_range.second;it1++) {
      q=it1->second; it2=bead_fil.find(q);
      k=it2->second;  // skip bonds to other filaments
      if ((it2!=bead_fil.end())&&(k==i)) { // pick the next, not prev bead
	if ((j==1)||((j>1)&&(q!=filSeq[j-2])))  
	  {filSeq[j]=q; break;}
      }
    }
  }

  if (cflag==0) { // linear
    assert(filSeq[filLength-1]==p_end[1]);
    n0=filLength;
  }
  else { // cyclic
    j=isBonded(filSeq[0],filSeq[filLength-1]); // check cyclic
    if (j!=2) error_dodri("findFilamentSeq: cyclic filament "+to_string(i)
			 +" has disconnected ends.");
    n0=-1*filLength; 
  }
  if (n0<0) assert(filLength>2); // sanity check for cyclic filament
  return n0;
}

/************************************************************************/
int fnet::getBondVec(int iBead, int iSkip, int *bondVecList)
{ /* Find list of beads bonded to iBead and save to bondVecList, 
     except for iSkip. Returns the number of bonded beads. */
  int n=0,jBead;
  multimap<int,int>::iterator it;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> it_range;  
  it_range=bond.equal_range(iBead);
  for (it=it_range.first;it!=it_range.second;it++) {
    jBead=it->second; // iBead-jBead form a bond
    if (jBead!=iSkip) {
      bondVecList[n]=jBead;
      ++n;
    }
  }
  return n;
}

/************************************************************************/
double fnet::getDist(int i, int j, int flag)
{ // get distance between bead i & j.
  // flag=2: distance squared; flag!=2: distance
  double r, x1,y1,x2,y2;
  x1=cm[i][0]; y1=cm[i][1]; x2=cm[j][0]; y2=cm[j][1];
  r=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
  r=(flag==2)?r:sqrt(r);
  return r;
}

/************************************************************************/
void fnet::getPxlEdge(int ix0,int iy0,unsigned char pxlCut,
		      vector<pair<double,double>> &edge_pos,img *img0)
{ /* Return a list of corner positions (ie edge) for pixel ix0,iy0 
     shared with a pixel of a different intensity. Checks all four 
     non-diagonal directions. Note in image, directions are flipped. 

     Skip edges on image boundary.
  */

  int ix,iy,idum;  
  double x0=(double)ix0,y0=(double)iy0;
  double e0x,e0y;

  //if (ix0==0 || ix0+1==(int)img0->rows || iy0==0 || iy0+1==(int)img0->columns)
  //  return; // skip pixels on boundary 
  
  ix=ix0; iy=iy0+1;   // top
  if (iy>=(int)img0->columns) goto bottom; 
  else  idum=(int)(img0->pxl_gray[ix][iy]);
  if (idum<pxlCut) {
    e0x=x0-0.5; e0y=y0+0.5; // corner 1
    edge_pos.push_back(pair<double,double>(e0x,e0y));
    e0x=x0+0.5; e0y=y0+0.5; //corner 2
    edge_pos.push_back(pair<double,double>(e0x,e0y));
  }
  
 bottom: 
  ix=ix0; iy=iy0-1;
  if (iy<=0) goto right; 
  else  idum=(int)(img0->pxl_gray[ix][iy]);
  if (idum<pxlCut) {
    e0x=x0-0.5; e0y=y0-0.5; // corner 1
    edge_pos.push_back(pair<double,double>(e0x,e0y));
    e0x=x0+0.5; e0y=y0-0.5; //corner 2
    edge_pos.push_back(pair<double,double>(e0x,e0y));
  }
  
 right: 
  ix=ix0+1; iy=iy0;
  if (ix>=(int)img0->rows) goto left; 
  else  idum=(int)(img0->pxl_gray[ix][iy]);
  if (idum<pxlCut) {
    e0x=x0+0.5; e0y=y0-0.5; // corner 1
    edge_pos.push_back(pair<double,double>(e0x,e0y));
    e0x=x0+0.5; e0y=y0+0.5; //corner 2
    edge_pos.push_back(pair<double,double>(e0x,e0y));
  }

 left: 
  ix=ix0-1; iy=iy0;
  if (ix<=0) return; 
  else   idum=(int)(img0->pxl_gray[ix][iy]);
  if (idum<pxlCut) {
    e0x=x0-0.5; e0y=y0-0.5; // corner 1
    edge_pos.push_back(pair<double,double>(e0x,e0y));
    e0x=x0-0.5; e0y=y0+0.5; //corner 2
    edge_pos.push_back(pair<double,double>(e0x,e0y));
  }
  
  return;   
}

/************************************************************************/
int fnet::isBonded(int ibead,int jbead)
{ /* Check bond between bead i and j. 
     returns: 0 (no bond). 2: bonded (i->j and j->i) */
  
  int ib=0, k,flag;
  multimap<int,int>::iterator it, it1, iter,jter;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator>
    it_range, it_range1;
  stringstream ss; 
  string static sdum="checkDuplicateBond: ";
  
  iter=bond.find(ibead);
  if (iter!=bond.end()) {// bond exists
    it_range=bond.equal_range(ibead); 
    it_range1=bond.equal_range(jbead);
    for (it=it_range.first;it!=it_range.second;it++) {
      if (it->second == ibead) {
	ss.clear(); ss << sdum<< "Self-bond found for bead "<< ibead;
	error_dodri(ss.str());
      }
      else if (it->second == jbead) {
	flag=0;
	for (it1=it_range1.first;it1!=it_range1.second;it1++) {
	  k=it1->second; if (k==ibead) flag++;
	}
	if (flag!=1) { // bonds not in a pair: this shouldn't happen
	  ss.clear(); 
	  if (flag==0) ss << sdum<<"Lone-pair bond from bead "<<ibead
			  <<" to "<<jbead;
	  else ss << sdum<<flag<<" bonds from bead " <<jbead<< " to "
		  <<ibead;
	  error_dodri(ss.str());
	}
	else ib=2; // normal
      }
      else continue; // another bond
    } //     for (it=it_range.first;it!=it_range.second;it++) {
  } // if (iter!=bond.end()) {// bond exists
  return ib;
} 

/************************************************************************/
void fnet::listCorners(double x0,double y0,int dir, double dr,
		       vector<pair<double,double>> &corner ) 
{ /* Return ordered list of pixel edges from corner pos (x0,y0). Position 
     (x0,y0) is at the intersection of 4 pixels, this function returns the
     attached corners of the 4 pixels from x0,y0 s.t. a pixel edge is from
     (x0,y0)-->corner. Does not exclude image bounds.     

     CW search (dir=0) or CCW search (dir=1).
     dr = search radius.
  */ 
  
  corner.push_back(pair<double,double>(x0+dr,y0));
  corner.push_back(pair<double,double>(x0,y0+dr));
  corner.push_back(pair<double,double>(x0-dr,y0));
  corner.push_back(pair<double,double>(x0,y0-dr));

  if (dir==1) { // reverse order
    reverse(corner.begin(),corner.end());
  }
  return; 
}

/************************************************************************/
void fnet::joinCyclicFil(double rCut)
{ /* Add bond to first/last bead of filament if it is cyclic 
     and bond is missing.   */ 

  int ibd,ibd1,nJoined,nBond0; 
  int iFil,filLength,*filSeq=new int[maxFilLength];
  double rCutSq=rCut*(double)beadDiam; //scale rCut by bead diameter
  rCutSq=rCutSq*rCutSq; 
  
  vector<img>::iterator it_img;
  it_img=findImg(img_tag);
  if (it_img==vimg.end())
    error_dodri("assignBond: cannot find image "+ img_tag);
  
  // close filaments if ends w/in bondL
  nJoined=0; // count # fils joined
  for (iFil=1;iFil<=nFilament;iFil++) {
    nBond0=nBond; 
    filLength=abs(findFilamentSeq(iFil,filSeq));
    if (filLength<=2) continue;
    ibd=filSeq[0]; ibd1=filSeq[filLength-1];
    registerOneBond(ibd,ibd1,-1,rCutSq,&(*it_img));
    if (nBond>nBond0) nJoined++; // a new bond was added
  }
  
  cout<<"  "<<nJoined<<" filaments joined."<<endl; 
  
  delete[] filSeq; 
  return;
}

/************************************************************************/
void fnet::measureWidth(double arr[],int nslice,string ofname,
			vector<array<double,2>> &axvec,
			vector<double> &axwdth,
			vector<array<double,2>> &epos0,
			vector<array<double,2>> &epos1,
			string ax0,double a0) 
{ /* Use 2 reference points to determine an line on the long
     axis of a slice. If arr is not set by user, this code 
     will automatically use the min/max positions of the BOC. Preferably, 
     auto locating arr is done if the model is oriented to either the x or y 
     axis. 

     Using line, calculate distances perpendicular to line up to where 
     extensions cross a filament boundary. Measures the total distance 
     of extension. */ 
  
  int npts=2; // use 2 points
  
  int i,j,k;
  double theta, **pt=new double*[2]; //pt[x,y][ipt]
  multimap<int,int>::iterator it;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> it_range;
  
  //finding cross points & calc distance
  int ib,flag,iBead,jBead;
  double p0[2],p1[2][2],rpos[2][2],axpt[2],cmtemp[2];  
  double step,dist,norm;
  array<double,2> widtharr;
  vector<array<double,2>> width1,width2;
  FILE *ofn;
  if (ofname=="") ofname="w_dist"; 
  ofname=addFileNameExt(ofname,".txt");
  ofn=fopen(ofname.c_str(),"w");

  // for locating arr
  int imin,imax;
  double min, max, sig;
  double *bm=new double[2],**pos;

  if (arr[0]!=-1) goto measure;
  
  /*********************/
  /* LOCATE reference points: boc centroid, max&min of x or y */
  // beads should prefereably be oriented along either x or y axis

  cout<<"  MESSAGE:: 3 reference points will be automatically located."<<endl;  
  
  pos=new double*[2];
  for (j=0;j<2;j++) pos[j]=new double[nBead];

  ib=0; 
  for (i=0;i<nBead;i++) {
    for (j=0;j<2;j++) pos[j][ib]=cm[i][j];
    ++ib;
  }
  
  for (j=0;j<2;j++) 
    getavgw(&bm[j],&sig,pos[j],beadMass,nBead);  // centroid
  if (ax0=="x") MinMax(nBead,pos[0],&imin,&imax,&min,&max); 
  else if (ax0=="y")  MinMax(nBead,pos[1],&imin,&imax,&min,&max);
  else error_dodri("fnet::measureWidth: invalid axis."); 

  if (ax0=="x") {
    arr[0]=pos[0][imin]; 
    arr[2]=pos[0][imax];
    arr[1]=arr[3]=bm[1]; // hold y pos
  }
  else if (ax0=="y") {
    arr[0]=arr[2]=bm[0]; // hold x pos
    arr[1]=pos[1][imin];
    arr[3]=pos[1][imax];
  }

  arr[4]=bm[0]; arr[5]=bm[1];   // 3rd point: centroid 
  /*********************/

 measure:

  //place arr into pt  : 1&2 will be ends, all rest are cross points
  for (j=0;j<2;j++) {pt[j]=new double[npts];}
  for (i=0;i<npts;i++) {
    for (j=0;j<2;j++) {
      pt[j][i]=arr[i*2+j];}
  }
  
  cout<<"  measure_width: npts= "<<npts<<" nSlice= "<<nslice<<
    " ofname= "<<ofname<<endl;
  cout<<"    reference points (x y z)= ";
  for (i=0;i<npts;i++) {
    cout<<"{";
    for (j=0;j<2;j++) 
      cout<<pt[j][i]<<" ";
    cout<<zCoord<<"} ";
  }
  cout<<endl;    

  // boc on x or y axis
  if (abs(pt[0][1]-pt[0][0])<(a0+.01)) { // aligned to y - vertical axis
    theta=PI/2.;
    norm=0;
    step=abs(pt[1][1]-pt[1][0])/(double)nslice;
  }
  else if (abs(pt[1][1]-pt[1][0])<(a0+.01)) { //  aligned to x - horizontal axis
    theta=(pt[1][1]-pt[1][0])/(pt[0][1]-pt[0][0]);    
    norm=PI/2.;
    step=abs(pt[0][1]-pt[0][0])/(double)nslice;
  }
  else {  
    theta=(pt[1][1]-pt[1][0])/(pt[0][1]-pt[0][0]);
    norm=-1./theta;     
    step=abs(pt[0][1]-pt[0][0])/(double)nslice;
  }

  fprintf(ofn,"#centerline \n");
  fprintf(ofn,"draw color gray \n");  
  fprintf(ofn,"draw line {%9.4f %9.4f %9.4f} {%9.4f %9.4f %9.4f} width 2\n",
	  pt[0][0],pt[1][0],zCoord,pt[0][1],pt[1][1],zCoord);
  fprintf(ofn,"draw color blue \n");    

  // start search at lower x/y end
  if (pt[1][0]<pt[1][1]) axpt[1]=pt[1][0]; 
  else axpt[1]=pt[1][1];
  if (pt[0][0]<pt[0][1]) axpt[0]=pt[0][0];
  else axpt[0]=pt[0][1]; 
  
  for (i=0;i<=nslice;i++) {
    if (norm!=0) axpt[1]=theta*(axpt[0]-pt[0][1])+pt[1][1]; // pt on defined axis

    rpos[0][0]=rpos[1][0]=axpt[0]; // start positions for search
    rpos[0][1]=rpos[1][1]=axpt[1];
    cmtemp[0]=cmtemp[1]=-1;
    
    flag=0;
    while (!flag) { // search until crossing found
      if (rpos[1][0]>rows || rpos[1][1]>columns) break; //no intersection found, stop search
      
      if (norm==PI/2.) rpos[1][1]+=step;  // horizontal axis
      else if (norm==0) rpos[1][0]+=step; // vertical axis
      else  rpos[1][1]=norm*(rpos[1][0]-axpt[0]) + axpt[1];
      
      cmtemp[0]=cmtemp[1]=-1;
      for (ib=0;ib<nBead;ib++) {
	it_range=bond.equal_range(ib);
	for (it=it_range.first;it!=it_range.second;it++) {
	  iBead=ib;jBead=it->second;
	  for (j=0;j<2;j++) {
	    p1[0][j]=cm[iBead][j];  p1[1][j]=cm[jBead][j]; 
	  }

	  flag=checkBondCrossing(rpos[0],rpos[1],p1[0],p1[1],p0);	  
	  if (flag) {
	    widtharr.at(0)=p0[0]; widtharr.at(1)=p0[1];
	    width1.push_back(widtharr);
	    cmtemp[0]=rpos[0][0]; cmtemp[1]=rpos[0][1];
	    goto endloop1;
	  }
	}
      }
      if (norm==PI/2.) rpos[1][0]=axpt[0];
      else rpos[1][0]+=step;
    endloop1:  ;
    } //while (!flag)
    
    
    //repeat loop, search in other direction 
    rpos[1][0]=axpt[0]; rpos[1][1]=axpt[1];
    
    flag=0;
    while (!flag) {     //until cross point
      if (rpos[1][0]<0 || rpos[1][1] <0) break; //no intersection
    
      if (norm==PI/2.)  rpos[1][1]-=step; // horizontal axis
      else if (norm==0) rpos[1][0]-=step; // vertical axis      
      else  rpos[1][1]=norm*(rpos[1][0]-axpt[0]) + axpt[1];
      
      for (ib=0;ib<nBead;ib++) {
    	it_range=bond.equal_range(ib);
    	for (it=it_range.first;it!=it_range.second;it++) {
    	  iBead=ib;jBead=it->second;
    	  for (j=0;j<2;j++) {
    	    p1[0][j]=cm[iBead][j];  p1[1][j]=cm[jBead][j]; 
    	  }	    

	  flag=checkBondCrossing(rpos[0],rpos[1],p1[0],p1[1],p0); 
    	  if (flag==1 && cmtemp[0]==rpos[0][0] && cmtemp[1]==rpos[0][1]) {
    	    fprintf(ofn,"draw sphere {%9.4f %9.4f %9.4f}\n",
    		    rpos[0][0],rpos[0][1],zCoord); 
    	    widtharr.at(0)=p0[0]; widtharr.at(1)=p0[1];
    	    width2.push_back(widtharr);
    	    goto endloop2;
    	  }
    	}
      } //   for (ib=0;ib<nBead;ib++) {
      if (norm==PI/2.) rpos[1][0]=axpt[0];
      else  rpos[1][0]-=step;
    endloop2: ;
    } //while (!flag)
    if (flag==0 && cmtemp[0]==rpos[0][0] && cmtemp[1]==rpos[0][1]) 
      width1.pop_back(); //if width1 made, but no corresponding width2

    if (norm!=0) axpt[0]+=step;
    else axpt[1]+=step; 
  } // for (i=0;i<=nslice;i++) 
  

  /* MEASURE DISTANCE SUM OF PERPENDICULAR EXTENSIONS FROM CENTERLINE */
  assert(width1.size()==width2.size());
  
  fprintf(ofn,"\n#width lines\n");
  for (k=0;k<(int)width1.size();k++) {
    fprintf(ofn,"draw line {%9.4f %9.4f %9.4f} {%9.4f %9.4f %9.4f} width 3 \n",width1[k].at(0),
	    width1[k].at(1),zCoord,width2[k].at(0),width2[k].at(1),zCoord);
    epos0.push_back(width1[k]);
    epos1.push_back(width2[k]);
  }
  fprintf(ofn,"\n#width measurements\n");
  for (k=0;k<(int)width1.size();k++) {
    dist=sqrt(pow(width1[k].at(0)-width2[k].at(0),2) +
	      pow(width1[k].at(1)-width2[k].at(1),2) );
    fprintf(ofn,"## %9.4f\n",dist);
    axwdth.push_back(dist); 
  }

  //MARK AND RETURN CENTERPOINTS 
  //getaxis block
  fprintf(ofn,"\n #get axis block\n");
  array<double,2> axarr;
  for (k=0;k<(int)width1.size();k++) {
    axarr.at(0)=(width1[k].at(0)+width2[k].at(0))/2;
    axarr.at(1)=(width1[k].at(1)+width2[k].at(1))/2;
    axvec.push_back(axarr);
  }
  for (i=0;i<(int)axvec.size();i++) {
    fprintf(ofn,"draw sphere {%9.4f %9.4f %9.4f} radius 2\n",axvec[i].at(0),
	    axvec[i].at(1),zCoord);
  }

  //clear memory
  for (i=0;i<2;i++) delete[] pt[i];  
  delete[] pt;

  fclose(ofn);
}

/************************************************************************/
void fnet::moveBeadPos(int p,double *newPos)
{ /* Calls bead::moveBeadPos() to fnet. */ 
  bead::moveBeadPos(p,newPos);
}

/************************************************************************/
void fnet::moveFilSeq(int *filSeq, int filLength, int i_fr, int i_to)
{ /* Move filSeq of length filLength from filament number i_fr to i_to. */
  
  int i,b;
  multimap<int,int>::iterator it;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> it_range;
  
  for (i=0;i<filLength;i++) {
    b=filSeq[i]; // bead number
    it=bead_fil.find(b); assert(it->second==i_fr);
    bead_fil.erase(it);
    it_range=fil_bead.equal_range(i_fr);
    for (it=it_range.first;it!=it_range.second;it++) {
      if (it->second==b) {fil_bead.erase(it); break;}
    }
    putBeadToFilament(b,i_to); // move b to filament i_to
  } // for (i=0;i<filLength;i++) {
}

/************************************************************************/
void fnet::putBeadToFilament(int iBead, int iFil)
{ /* Register iBead to iFil. NOTE: bonds of iBead to other beads in iFil
     have to be formed separately. */
  
  multimap<int,int>::iterator it,jt;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> it_range;

  it=bead_fil.find(iBead);
  if (it !=bead_fil.end() ){ // iBead belongs to a filament
    if (it->second == iFil) return; // iBead already on iFil
    else { // another filament. Delete iBead from it
      it_range=fil_bead.equal_range(it->second);
      for (jt=it_range.first;jt!=it_range.second;jt++) {
	if (jt->second == iBead) {
	  fil_bead.erase(jt); break;
	}
      }
      bead_fil.erase(it);
    }
  }
  bead_fil.insert( pair<int,int>(iBead,iFil) );
  fil_bead.insert( pair<int,int>(iFil,iBead) );

  return;
}

/************************************************************************/
void fnet::readData(string fname, string bd_tag0, ifstream *fin)
{ /* Reads in bead data file to fill data fields of new fnet object. */
  
  stringstream ss;
  string bd_tag, line,count_bnd;
  int iBead,jBead,jBead0; 
  int ivec[1],ivec1[1],nBeads_temp=0,nBonds_temp=0; 
  bool newFile = (fin==NULL); 
  bool b_tag, nbond,flag;
  size_t iCpos, fCpos; 

  multimap<int,int>::iterator it;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> it_range; 
  vector<bead>::iterator it_bead; 
  flag = b_tag = nbond = false ;

  //clear vars for this tag 
  bond.clear(); bead_fil.clear(), fil_bead.clear();
  nBond=0;  

  if (newFile) {
    fin = new ifstream;
    fin->open(fname.c_str(),ios::in);
  }
  if (fin->bad()) error_dodri("dodri: Invalid fnet file"); 
  //read file header 
  while (!flag && fin->good()) {
    flag = b_tag &&nbond ;
    if (flag) break ;

    if (line.find("Bead ")!=std::string::npos)
      error_dodri("READ:: Incomplete bead header");
    getline(*fin,line);

    if (line == "") continue;
    ss.str(line);

    //read nBeads
    if (line.find("nBeads=")!=string::npos) {
      getCommOptInt(ss,"nBeads=",ivec,1,-1);
      nBeads_temp=ivec[0];
      if (ivec[0]<=0) {
    	cout<<"  "<<std::right<<"MESSAGE: file contains 0 beads."<<endl;
    	return;
      }
      else continue; 
    }

    //read bead_tag
    if (!b_tag && line.find("bead_tag=")!=string::npos) {
      if (bd_tag0=="")
      	bd_tag=getCommOptString(ss,"bead_tag=","");
      else bd_tag=bd_tag0;
      if (bd_tag!="") {
    	bead_tag = bd_tag;
    	b_tag = true;
    	continue;
      }
      else error_dodri("READ: missing bead tag");
    }      
    it_bead = findBead(&vbead,bead_tag); //check that bead_tag exists in vbead
    if (it_bead==vbead.end()) error_dodri("READ:: bead tag not found"); 
    
    //read nBonds
    if (line.find("nBonds=")!=string::npos) {
      getCommOptInt(ss,"nBonds=",ivec1,1,-1);
      nBonds_temp=ivec1[0];
      if (ivec1[0]<=0) {
    	cout<<"  "<<std::right<<"WARNING: file contains 0 bonds."<<endl;
	return;
      }
      else nbond=true; 
    }
  } //   while (!flag && fin->good()) {

  assert(it_bead->nBead==nBeads_temp); //check bead objects match

  //read and register bonds
  while (!ss.fail() && getline(*fin,line)) {
    if (!isdigit(fin->peek())) break;     
    if (line=="" || line.find("Bead ")!=std::string::npos) continue;
    ss.str(line);
    ss >> iBead >> count_bnd; 
    if (count_bnd=="0") continue; 

    iCpos = line.find(" "+count_bnd+" ");    
    fCpos = (count_bnd+" ").length(); 
    line.erase(line.begin(),line.begin()+iCpos+fCpos);

    jBead=-1, jBead0=-2; // initiate 
    while (!line.empty()) { //register count_bnd bonds with iBead 
      ss.str(line);
      ss >> jBead;
      iCpos = line.find(" "+to_string(jBead)+" ");
      fCpos = (" "+to_string(jBead)+" ").length();
      line.erase(line.begin(),line.begin()+iCpos+fCpos);
      if (jBead0==jBead) break; 
      jBead0=jBead;

      flag=false; 
      if (bond.find(iBead)!=bond.end()) {   //skip if bond exists
    	it_range=bond.equal_range(iBead);
    	for (it=it_range.first;it!=it_range.second;it++) {
    	  if (it->second==jBead) { flag=true; break; }
    	}
      } //       if (bond.find(iBead)!=bond.end()) {   //skip if bond exists
      if (flag) continue; //bond exists
      registerOneBond(iBead,jBead,-1,-1,&(*findImg(it_bead->img_tag))); //saves 2 bonds
    } //    while (!line.empty()) { //register count_bnd bonds with iBead  
    
    ss.clear(); 
    if (ss.fail()) error_dodri("READ:: error in reading bonds");    
    //    if (!isdigit(fin->peek())) break; 
  } //   while (!ss.fail() && getline(*fin,line)) {
  
  assert(nBonds_temp==nBond);   //assert that multimap size==nbonds to confirm
  
  cout<<"  read:: nBonds= "<< nBond<<endl;
  if (newFile) fin->close(); 
  
  return;
}

/************************************************************************/
int fnet::registerOneBond(int i, int j, int iCut, double rCutSq,img *img)
{ /* Put a bond between bead i and j if their distance squared  
    is less than rCutSq.

    iCut>=0: Do not form a bond if it crosses an area w/ intensity in 
      img <= iCut.
    iCut<0: No intensity cutoff    */
  
  double xi,xj,yi,yj,slope,b;
  int *ipos=new int[2],ix,iy;

  if (getDist(i,j,2)<rCutSq) { 
    if (iCut>0) { // apply pixel intensity cutoff
      xi=cm[i][0];yi=cm[i][1];xj=cm[j][0];yj=cm[j][1];
      slope=(yj-yi)/(xj-xi);
      if (fabs(slope)<1) {// scan horizontal
	ipos[0]=(xi<xj)?(int)xi:(int)xj;
	ipos[1]=(xi<xj)?(int)xj:(int)xi;
	slope=1.*slope;
	b=yi-slope*xi;
	for (ix=ipos[0];ix<=ipos[1];ix++) {
	  iy=(int)(slope*(double)ix)+(int)b;
	  if (img->pxl_gray[ix][iy]<=(unsigned char)iCut) return -1;
	}
      } // if (fabs(slope)<1) {// scan horizontal
      else { // scan vertical
	ipos[0]=(yi<yj)?(int)yi:(int)yj;
	ipos[1]=(yi<yj)?(int)yj:(int)yi;
	slope=1./slope; b=xi-(slope*yi);
	for (iy=ipos[0];iy<=ipos[1];iy++) {
	  ix=(int)(slope*(double)iy)+(int)b;
	  if (img->pxl_gray[ix][iy]<=(unsigned char)iCut) return -1;
	}
      }
    } // if (iCut>0) { // apply pixel intensity cutoff
    bond.insert( pair<int,int>(i,j) );   
    bond.insert( pair<int,int>(j,i) );  nBond+=2;
   } //   if (getDist(i,j,2)<rCutSq) {

  if (rCutSq<0) {   // add bond w/o dist check
    bond.insert( pair<int,int>(i,j) );  
    bond.insert( pair<int,int>(j,i) );  nBond+=2; 
  }
  
  delete[] ipos;
  return 0;
}

/************************************************************************/
void fnet::removeFloatingFil(int length_cut) 
{ /* Delete filaments with length<=length_cut.  */

  int i,j,iBead,nBead0;
  int iFil,filLength,*filSeq=new int[maxFilLength];
  vector<int> markbead;
  
  int nFil0=nFilament, nFil=0;   
  markbead.clear();
  
  // check filament lengths, mark beads for deletion
  for (iFil=1;iFil<=nFilament;iFil++) {
    filLength=findFilamentSeq(iFil,filSeq);
    filLength=abs(filLength); // take care of cyclic
    if (filLength<=length_cut) {
      nFil=nFil+1;
      for (i=0;i<filLength;i++) { // mark beads for deletion
	iBead=filSeq[i];
	markbead.push_back(iBead);
      }
    }
  }
  
  // delete marked beads 
  for (i=0;i<(int)markbead.size();i++)  {
    nBead0=nBead-1;
    for (j=i+1;j<(int)markbead.size();j++) {
      if (markbead[j]==nBead0) {
	// update ibead in markbead since deleteBead will re-index beads
  	markbead[j]=markbead[i]; 
  	break;
      }
    }
    deleteBead(markbead[i],1);
  }

  if (!markbead.empty()) {
    for (iFil=1;iFil<nFilament;iFil++) {
      if ((int)fil_bead.count(iFil)==0) { // replace skippedfil during bead deletion
	filLength=findFilamentSeq(nFilament,filSeq);
	assert (filLength!=0);
	moveFilSeq(filSeq,filLength,nFilament,iFil);
	break;
      }
    }
    if ((int)fil_bead.count(nFilament)==0) --nFilament;

    for (iFil=1;iFil<=nFilament;iFil++) 
      assert((int)fil_bead.count(iFil)!=0); // verify no empty filaments
  }
  
  cout<<"  clean_filaments: segLength= "<<length_cut<<endl;
  cout<<"         filaments deleted= "<<nFil0-nFilament<<
    "; nFilament="<<nFilament<<endl; 
  
  delete[] filSeq; 
  return;
}

/************************************************************************/
void fnet::removeOneBond(int i, int j)
{ /* Remove bonds (i,j) and (j,i) */
  
  multimap<int,int>::iterator it;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> it_range;
  
  assert(i!=j);
  
  it_range=bond.equal_range(i);
  for (it=it_range.first;it!=it_range.second;it++) {
    if (it->second ==j ) {
      bond.erase(it); --nBond;  break;
    }
  }
  it_range=bond.equal_range(j);
  for (it=it_range.first;it!=it_range.second;it++) {
    if (it->second == i ) {
      bond.erase(it); --nBond;  break;
    }
  }
}


/************************************************************************/
void fnet::smoothenFilAverage(int segLength,double rCut)
{ /* Smooth out filaments. Bead is moved to average displacement
     from average fit of segLength number of beads.   */
  
  int iFil, filLength,filLength0; 
  int *filSeq=new int[maxFilLength],*filSeq0=new int[maxFilLength];
  int i,j,k,iBead; 
  double a1,b1,siga1,sigb1,chi2_1; 
  signed char flag; 
  double *x=new double[segLength],*y=new double[segLength],*pos=new double[2];
  double uvec0[2],qvec0[2],vdum[2], vdum1[2];
  double rCut1=rCut*(double)beadDiam, rCutSq=rCut1*rCut1;
  
  if (rCut<0) rCut=(double)rows; // no rcut max; set to rows for safety 
  
  multimap<int,int>::iterator it, it1;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> 
    it_range, it1_range;

  // save sum of pos changes 
  double cm_temp[nBead][2], nCount[nBead];
  for (i=0;i<nBead;i++) {
    for (k=0;k<2;k++) { cm_temp[i][k]=0.; nCount[i]=0.; }
  }
  
  cout<<"  smoothen_average: segLength= "<<segLength<<" rCut= "<<rCut<<endl; 
  
  for (iFil=1;iFil<=nFilament;iFil++) {
    filLength=findFilamentSeq(iFil,filSeq);
    if (abs(filLength)<=segLength*2) continue; 

    if (filLength<0) { // cyclic fil;
      // append first #seglength beads onto end of filSeq
      filLength=abs(filLength);
      filLength0=filLength+segLength; // incr fillength for filSeq0

      for (i=0;i<filLength0;i++) { // populate filSeq0
	if (i>=filLength) {
	  filSeq0[i]=filSeq[i-filLength];
	}
	else filSeq0[i]=filSeq[i];
      }
      
      for (i=0;i<filLength;i++) { 
	flag=1; // flag=1 if y=a+bx is used, flag=-1 if x=a+by is used	
	assignFilSeqPos(&filSeq0[i],x,y,segLength,1);
	fit0(x,y,segLength,&a1, &b1, &siga1, &sigb1, &chi2_1);
	if ((a1!=a1)||(b1!=b1)||(fabs(b1)>1)) { // do x=a+by
	  // a1!=a1 when a1==nan. same for b1
	  fit0(y,x,segLength,&a1, &b1, &siga1, &sigb1, &chi2_1);
	  flag=-1;
	}	
	line2vec(a1,b1,uvec0,qvec0,flag);
	for (j=0;j<segLength;j++) {
	  iBead=filSeq0[i+j]; 
	  for (k=0;k<2;k++) vdum[k]=cm[iBead][k]; 
	  distPoint2Line(vdum,uvec0,qvec0,2,pos,2); 
	  for(k=0;k<2;k++) cm_temp[iBead][k]+=(pos[0]*uvec0[k]+qvec0[k]);
	  ++nCount[iBead];	
	}
      }
      for (i=0;i<filLength;i++) {
	iBead=filSeq[i];
	for (k=0;k<2;k++) vdum1[k]=cm_temp[iBead][k]/nCount[iBead];
	if (getDistPos(cm[iBead],vdum1,2)<rCutSq)  // move bead
	  bead::moveBeadPos(iBead,vdum1);
      }      
    } // if (filLength<0) {
    
    else {    
      // find local fits
      for (i=0;i<filLength-(segLength-1);i++) { 
	flag=1; // flag=1 if y=a+bx is used, flag=-1 if x=a+by is used
	assignFilSeqPos(&filSeq[i],x,y,segLength,1);
	fit0(x,y,segLength,&a1, &b1, &siga1, &sigb1, &chi2_1);
	if ((a1!=a1)||(b1!=b1)||(fabs(b1)>1)) { // do x=a+by
	  // a1!=a1 when a1==nan. same for b1
	  fit0(y,x,segLength,&a1, &b1, &siga1, &sigb1, &chi2_1);
	  flag=-1;
	}	
	line2vec(a1,b1,uvec0,qvec0,flag);      	
	for (j=0;j<segLength;j++) {
	  iBead=filSeq[i+j]; 
	  for (k=0;k<2;k++) vdum[k]=cm[iBead][k]; 
	  distPoint2Line(vdum,uvec0,qvec0,2,pos,2); 
	  for(k=0;k<2;k++) cm_temp[iBead][k]+=(pos[0]*uvec0[k]+qvec0[k]);
	  ++nCount[iBead];	
	}
      }
      for (i=0;i<filLength;i++) {
	iBead=filSeq[i];
	for (k=0;k<2;k++) vdum1[k]=cm_temp[iBead][k]/nCount[iBead];
	if (getDistPos(cm[iBead],vdum1,2)<rCutSq) 
	  bead::moveBeadPos(iBead,vdum1);
      }
    } // else {
  } //   for (iFil=1;iFil<=nFilament;iFil++) {
  
  // free up memory
  delete [] filSeq; delete [] filSeq0;
  delete [] x; delete [] y; delete [] pos;
  
  return; 
}

/************************************************************************/
void fnet::updateFilament(int iBead,int kBead)
{ /* Called after iBead and kBead are joined and respective filaments
     need to be updated. Will split a filament if iBead and 
     kBead belong to the same filament and iBead--kBead are not bonded. 

     Bonds NOT modified here. Any forming/breaking of bonds should
     be done elsewhere.   */
  
  int i,iloc=-1,jloc=-1,maxLoc;
  int iFil,kFil,iFilLength,kFilLength,jFilLength;
  int *iFilSeq=new int[maxFilLength], *kFilSeq=new int[maxFilLength],
    *jFilSeq=new int[maxFilLength]; 
  
  iFil=bead_fil.find(iBead)->second;
  kFil=bead_fil.find(kBead)->second;

  // check if filament needs to be split
  if (iFil==kFil) { // && iFil!=0) {
    if (iFil!=0) {
      if (!isBonded(iBead,kBead)) { // split filament :
	// ensure bond between iBead and kBead was removed
	// cp kFil beads from {maxLoc,end} to new filament jFil
	// temporarily add bond to prevent issues in findFilamentSeq
	registerOneBond(iBead,kBead,-1,-1,&(*findImg(this->img_tag)));
	kFilLength=findFilamentSeq(kFil,kFilSeq);

	if (kFilLength<0) { //cyclic, do nothing
	  removeOneBond(iBead,kBead); // remove temp added bond
	  return; 
	} // else, split
	kFilLength=abs(kFilLength);    
	
	for (i=0;i<kFilLength;i++) { 
	  if (kFilSeq[i]==kBead) iloc=i;
	  else if (kFilSeq[i]==iBead) jloc=i; 
	}
	maxLoc=(iloc>jloc)?(iloc):(jloc); 
	for (i=maxLoc;i<kFilLength;i++)  jFilSeq[i-maxLoc]=kFilSeq[i]; 
	jFilLength=kFilLength-maxLoc;
	
	assert(jFilLength>0); // sanity
	// if single bead on either end: remove that bead from kfil, its fil->0
	if (maxLoc==1) putBeadToFilament(kFilSeq[0],0);
	else if (maxLoc==kFilLength-1) putBeadToFilament(kFilSeq[maxLoc],0); 
	else moveFilSeq(jFilSeq,jFilLength,kFil,++nFilament); 
	removeOneBond(iBead,kBead);     	// remove temp added bond
	return ;
      }
      else return; // do nothing
    }
    else if (iFil==0 && kFil==0) { // make a new filament
      ++nFilament; 
      putBeadToFilament(iBead,nFilament);
      putBeadToFilament(kBead,nFilament);
      return; 
    }
  }
  
  // else: update existing filaments
  if (iFil==0 && kFil!=0)  // insert iFil to kFil
    putBeadToFilament(iBead,kFil); 
  else if (iFil!=0 && kFil==0)  // insert kBead to iFil
    putBeadToFilament(kBead,iFil); 
  else if (iFil==0 && kFil==0) {  // create new filament
    ++nFilament; 
    putBeadToFilament(iBead,nFilament);
    putBeadToFilament(kBead,nFilament); 
  }
  else { // handle existing filaments
    iFilLength=findFilamentSeq(iFil,iFilSeq);
    kFilLength=findFilamentSeq(kFil,kFilSeq);
    
    iFilLength=abs(iFilLength); // cyclic
    kFilLength=abs(kFilLength);    

    if (iFil<kFil) { // move kFil>>iFil,nFil>kFil
      moveFilSeq(kFilSeq,kFilLength,kFil,iFil);
      if (kFil<nFilament) { // reduce # of filaments
	jFilLength=findFilamentSeq(nFilament,jFilSeq);
	jFilLength=abs(jFilLength);
	moveFilSeq(jFilSeq,jFilLength,nFilament,kFil);
      }
      --nFilament;
    } //     if (iFil<kFil) { // move kFil>>iFil,nFil>kFil
    else { // iFil>kFil
      moveFilSeq(iFilSeq,iFilLength,iFil,kFil);
      if (iFil<nFilament) { // reduce # of filaments
	jFilLength=findFilamentSeq(nFilament,jFilSeq);
	jFilLength=abs(jFilLength);	
	moveFilSeq(jFilSeq,jFilLength,nFilament,iFil);
      }
      --nFilament;
    }
  }
  
  return; 
}

/************************************************************************/
void fnet::writeCor(string fname, string fmt, string occ, string outMode)
{ /* fmt=output format: "cor" or "pdb"
     occ=option for occupancy
     "mass": bead mass, "angle": filament angle, "length": filament length  */
  
  FILE *ff;
  string ofname,segname;
  double x1, y1, z1, weightMin, weightMax;
  int imin, imax;
  multimap<int,int>::iterator it;
  int i,resid,*filSeq=new int[maxFilLength]; 
  double r1, *x=new double[maxFilLength], *y=new double[maxFilLength];
  double *weight=new double[nBead];
  
  for (i=0;i<nBead;i++) weight[i]=0.;
  if ((fmt!="cor")&&(fmt!="pdb")) 
    error_dodri("writeCor: unknown output format: "+fmt);
  
  ofname=addFileNameExt(fname,"."+fmt);
  ff=fopen(ofname.c_str(),"w");
  cout<<"  Write fnet coor to "<<ofname<<"   with occ= "
      <<occ<<" outMode= "<<outMode<<endl;
  
  // find min/max for occupancy (col 55-60)
  if (occ=="mass") {
    for (i=0;i<nBead;i++) weight[i]=beadMass[i];
  }
  else if ((occ=="length")||(occ=="angle")) {
  } // else if ((occ=="length")||(occ=="angle")) {
  else {}
  
  MinMax(nBead, weight, &imin, &imax, &weightMin, &weightMax);
  
  if (fmt=="cor") {    // use extended format for more than 99999 atoms
    if (nBead>99999) fprintf(ff,"%10d%s\n",nBead,"  EXT");
    else fprintf(ff,"%d\n",nBead);
  }
  else if (fmt=="pdb") fprintf(ff,"MODEL\n");
  else {}
  
  if (nFilament==0) { segname=tag.substr(0,4); } // tag as 4-letter segid
  resid=1;
  
  for (i=0;i<nBead;i++) {
    if (nFilament>0) { // Use filament number for segid
      resid=bead_fil.find(i)->second; 
      stringstream ss;
      ss<< hex << resid ;// express in hexadecimal 
      segname=ss.str();
    }
    
    if (outMode=="verbose") {x1=cm[i][0]; y1=cm[i][1]; z1=zCoord;} // verbose
    // matches w/ image
    else {y1=((double)rows-cm[i][0]); x1=cm[i][1]; z1=zCoord;} 
    r1=100.*weight[i]/weightMax;
    
    /* From coorio.src:
       fm2='(2I5,1X,A4,1X,A4,3F10.5,1X,A4,1X,A4,F10.5)'
       I,IRES,RES(IRES),ATYPE(I),Xyz(1,I),xYz(2,I),xyZ(3,I),SID,RID,WMAIN(I)
    */
    
    if (fmt=="cor") {
      if (nBead>99999) // use extended format
	fprintf (ff,"%10d%10d  %-8s  %-8s%20.10f%20.10f%20.10f  %-8s  %-8d%20.10f\n",
		 i+1,resid,"FNET","O",x1,y1,z1,segname.c_str(),resid,r1);
      else 
	fprintf (ff,"%5d%5d %-4s %-4s%9.4f %9.4f %9.4f  %4s %-4d%10.5f\n",
		 i+1,resid,"FNET","O",x1,y1,z1,segname.c_str(),resid,r1);
    } // if (fmt=="cor") {
    else if (fmt=="pdb")
      fprintf (ff,"ATOM  %5d  O   FNET %04d    %8.3f%8.3f%8.3f%6.2f 0.00\
    %4s\n",  i+1,resid,x1,y1,z1,r1,segname.c_str());
    else {}
    
  }
  if (fmt=="pdb") fprintf(ff,"ENDMDL\n");
  fclose(ff);
  
  delete filSeq; delete x; delete y; delete weight;
}

/************************************************************************/
void fnet::writeData(string fname, string bd_fname, FILE *dat, FILE *bd_dat)
{ /* Writes fnet data to a dat file with name fname. 
     Dat file contains nBond, bond list. Filament info not saved.   */ 
  
  int iBead,jBead;
  bool newFile = (dat==NULL);
  string datname,bd_tag0,sdum="bd",sdum1,id;
  multimap<int,int>::iterator it;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> it_range; 
  
  vector<bead>::iterator it_bead;
  bool bd_newFile = (bd_dat == NULL);
  string bd_datname; 
  
  std::locale loc;
  id="";
  for (int i=0;i<(int)img_tag.size();i++) {
    if (isdigit(img_tag[i],loc)) {
      sdum1=img_tag[i];
      id=id+sdum1;
    }
  }
  bd_tag0=sdum+id;
  
  datname = addFileNameExt(fname,".dat");   /// write fnet bonds
  if (newFile) {
    dat = fopen(datname.c_str(),"w");
    cout<<"  Write fnet bonds to "<<datname<<endl;
  }

  if (bead_tag!="") 
    fprintf(dat,"bead_tag= %s\nnBeads= %d\ntag= %s\nnBonds= %d\n",
	    bead_tag.c_str(),nBead,tag.c_str(),nBond);
  else  // no bead tag from fnet contour build
    fprintf(dat,"bead_tag= %s\nnBeads= %d\ntag= %s\nnBonds= %d\n",
	    bd_tag0.c_str(),nBead,tag.c_str(),nBond);
  fprintf(dat,"Bead    Count   List \n");
  
  for (iBead=0;iBead<nBead;iBead++) {
    fprintf(dat,"%-8d",iBead);
    fprintf(dat,"%-8d",(int)bond.count(iBead));    
    
    it_range=bond.equal_range(iBead); //print all bonds to this bead 
    for (it=it_range.first;it!=it_range.second;it++) {
      jBead=it->second;
      fprintf(dat,"%-8d",jBead);
    }
    fprintf(dat,"\n"); //print newline at end of this iBead's bond list
  }
  
  // Write bead file also
  bd_datname = addFileNameExt(bd_fname,".dat");
  if (bd_newFile) {
    bd_dat = fopen(bd_datname.c_str(),"w");
    cout<<"  Write bead data to "<<bd_datname<<endl;
  }
  
  fprintf(bd_dat, "img_tag= %s\ntag= %s\nmaxBead= %d\nnBeads= %d\nbeadDiam= %d\nrows= %d\ncolumns= %d\ndRow= %.2f\ndCol= %.2f\n",
	  img_tag.c_str(), bd_tag0.c_str(), maxBead, nBead,
	  beadDiam, rows, columns, dRow, dCol); 
  
  fprintf(bd_dat, "Bead     Mass               X                  Y                  Z\n");
  for (iBead = 0; iBead < nBead; iBead++) {
    fprintf(bd_dat, "%-8d %-15.12e %-15.12e %-15.12e %-15.12e \n",
	    iBead, beadMass[iBead], cm[iBead][0], cm[iBead][1], zCoord);
  }
  
  if (newFile) fclose(dat);
  if (bd_newFile) fclose(bd_dat); 
  
  return; 
}

/************************************************************************/
void fnet::writeImage(int beadL, int bondL, int nColor, double *colorBg,
		      double *colorBead, double *colorBond,
		      string fname, string fmt)
{ /* Write bead coordinates and fnet lines to an image file.
     beadL,bondL: Sizes of beads and bonds. If <=0, do not draw the 
     corresponding element. 
     nColor=1: grayscale, nColor=3: rgb
     color_{bg,bead,bond}: corresponding colors.
     fname,fmt: output file name.  */
  
  int i,flag;
  string ofname;
  
  Image img0(Geometry(0,0),"black");
  list<Magick::Drawable> drawBead,drawBond; //construct drawing list
  
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
    cout<<"  "<<"WARNING: nBead= 0. No image written"<<endl;
    return;
  }
  if (beadL<=0) {
    if (bondL<=0) {
      cout<<"  "<<"WARNING: Zero bead/bond sizes. No image written"<<endl;
      return;
    }
    cout<<"  Beads not rendered."<<endl;
  }
  else {
    if (bondL<=0) cout<<"  Bonds not rendered."<<endl;
  }
  
  ///////////////////////////////////////////////////////
  // Write parameters and render
  ///////////////////////////////////////////////////////
  
  cout<<"  Write image to "<<ofname<<endl;
  cout<<"  Bead/Bond sizes: "<<setw(5)<<beadL<<"  "<<bondL<<endl;
  
  if (nColor==3) { // rgb
    cout<<endl<<"    RGB colors:"<<endl;
    cout<<"    Background: ";
    for (i=0;i<3;i++) cout<<setw(5)<<setprecision(2)<<colorBg[i]<<" ";
    cout<<endl;
    if (beadL>0) {
      cout<<"    Bead: ";
      for (i=0;i<3;i++) cout<<setw(5)<<setprecision(2)<<colorBead[i]<<" ";
      cout<<endl;
    }
    if (bondL>0) {
      cout<<"    Bond: ";
      for (i=0;i<3;i++) cout<<setw(5)<<setprecision(2)<<colorBond[i]<<" ";
      cout<<endl;
    }
    img0.extent(Geometry(columns,rows),ColorRGB(colorBg[0],colorBg[1],
						colorBg[2]));
    img0.modifyImage(); //ensure there is only 1 reference to underlying image
  } //   if (nColor==3) { // rgb
  else { // gray: only use the first color
    cout<<endl<<"    Grayscale colors:"<<endl;
    cout<<"    Background: "<<setw(5)<<setprecision(2)<< *colorBg <<endl;
    if (beadL>0) {
      cout<<"    Bead: "<<setw(5)<<setprecision(2)<< *colorBead <<endl;
    }
    if (bondL>0) {
      cout<<"    Bond: "<<setw(5)<<setprecision(2)<< *colorBond <<endl;
    }
    img0.extent (Geometry(columns,rows), ColorGray(*colorBg)); //base image
    img0.modifyImage(); //ensure there is only 1 reference to underlying image
  } //   else { // gray: only use the first color
  
  img0.matte(true);   
  writeImageBead(drawBead,cm,nBead,beadL, nColor, colorBead);
  writeImageBond(drawBond,cm,nBead,bondL, bond, nColor, colorBond);
  img0.draw(drawBond); //draw all objects on list
  img0.draw(drawBead); //draw all objects on list
  img0.depth(8); //write to 8-bit image
  img0.transparent(ColorGray(*colorBg)); //add transparency 
  img0.write(ofname);  

  return; 
}

/************************************************************************/
void fnet::writePsf(string fname)
{
  int j,n,resid, **btemp;
  FILE *psf;
  string segname="",ofname=fname, fmt="psf";
  stringstream ss;
  multimap<int,int>::iterator it;

  // determine filename extension
  std::size_t pos=fname.find(fmt);
  j=(int)pos;
  if (j<0) ofname=fname+"."+fmt;
  else ofname=fname;

  psf=fopen(ofname.c_str(),"w");

//  std::size_t pos = fname.find(".psf");
//  if (fname.length()- pos != 4) ofname=fname+".psf";
//  else ofname=fname;
  cout<<"  Write fnet psf to "<<ofname<<endl;
  
  //fprintf(psf,"PSF CMAP CHEQ\n\n !NTITLE\n");
  fprintf(psf,"PSF \n\n !NTITLE\n");
  fprintf(psf,"* PSF for %s \n* \n\n",fname.c_str());

  fprintf(psf,"%6d !NATOM\n",nBead);
  btemp=new int*[nBond/2]; //initialize w/ or w/o bonds, cpp warning
  if (nBond>0) { // case w/ bonds
    //    btemp=new int*[nBond/2];
    for (j=0;j<(nBond/2);j++) btemp[j]=new int[2];
  }
  for (j=0;j<(nBond/2);j++) { btemp[j][0]=-9999;btemp[j][1]=-9999;}

  if (nFilament==0) { segname=tag.substr(0,4); } // tag as 4-letter segid
  resid=1;

  for (j=0;j<nBead;j++) { 
    if (nFilament>0) { // Use filament number for segid
      resid=bead_fil.find(j)->second; 
      stringstream ss;
      ss<< hex << resid ;// express in hexadecimal 
      segname=ss.str();
    }
    fprintf(psf,"%8d %4s %04d FNET O       0 %14.6f %14.6f %d\n",
	    j+1,segname.c_str(),resid,0.0,beadMass[j],0);
    /* from psfres.src: fmt01=(I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8)
                   I(8),SEGID(ISEG)(4),RESID(IRES)(4),RES(IRES)(4),
		   ATYPE(I)(4),IAC(I)(4; atomID),CG(I)(14.6;charge),
		   AMASS(I)(14.6;mass),IMOVE(I)(8)
    
        fprintf(psf,"%8d %4s %04d FNET  O     0   0.000000        1.0000 \
               0   0.00000      0.000000E-02\n",j+1,segname.c_str(),k );
        //    fprintf(psf,"%8d %04d %04d NET O 0 0.000000 1.0000 0 0.00000
        //    0.000000E-02\n",j+1,k,k );
    */
  }
  fprintf(psf,"\n%8d !NBOND\n",nBond/2);
  if (nBond==0) goto WRITEPSF0;

  n=0;
  for (it=bond.begin();it!=bond.end();it++) {
    for (j=0;j<n;j++) { // check for duplicate bond
      if ((btemp[j][0]==it->first)&&(btemp[j][1]==it->second)) goto WRITEPSF1;
      if ((btemp[j][1]==it->first)&&(btemp[j][0]==it->second)) goto WRITEPSF1;
    }
    btemp[n][0]=it->first; btemp[n++][1]=it->second;
  WRITEPSF1: continue;
  }
  assert (n==(nBond/2));
  for (j=0;j<n;j++) {
    fprintf(psf,"%8d%8d",btemp[j][0]+1,btemp[j][1]+1);
    if(j%4==3) fprintf(psf,"\n");
  }

 WRITEPSF0: fclose(psf);
  // free up memory
  if (nBond>0) {
    for (j=0;j<nBond/2;j++) delete [] btemp[j];
  }
  delete [] btemp;

  return;
}

/************************************************************************/
// Non-member functions
/************************************************************************/

void checkDuplicateFnetTag(vector<fnet> * vfnet,string tag)
{ /* In vfnet, find iterator for the element with tag */
  
  vector<fnet>::iterator it_fnet;
  for (it_fnet=(*vfnet).begin();it_fnet!=(*vfnet).end();++it_fnet) {
    if ( (*it_fnet).tag==tag) error_dodri("Duplicate fnet tag: " + tag);
  }
  return;
}

/************************************************************************/
vector<fnet>::iterator findFnet(vector<fnet> * vfnet,string tag)
{ /* In vfnet, find iterator for the element with tag */
  
  vector<fnet>::iterator it_fnet;
  for (it_fnet=(*vfnet).begin();it_fnet!=(*vfnet).end();++it_fnet) {
    if ( (*it_fnet).tag==tag)  break;
  }
  return it_fnet;
}
