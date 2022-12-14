#include "dodri.h"

/************************************************************************/
// Constructors
/************************************************************************/

fnet3d::fnet3d(img3d *img3d0, int iDiam,unsigned char pxlCut,string tag0)
{ /* constructor directly from img0. */

  int ndigit;
  string img_tag,fnet_tag; 
  stringstream ss; 
  vector<img>::iterator it_img;
  vector<fnet>::iterator it_fnet; 
   
  nStack=img3d0->nStack;
  img3d_tag=img3d0->tag; 
  tag=tag0; bead3d_tag="";
  
  // fnet  
  bond_a=new multimap<int,int>[nStack];
  bond_b=new multimap<int,int>[nStack];

  ndigit=nDigit(nStack);
  
  for (int i=0;i<nStack;i++) {
    ss.str(""); ss.clear();
    ss<< tag<<i; fnet_tag=ss.str();
    img_tag=img3d0->img_stack[i].tag;
    it_img=findImg(img_tag); 
    cout<<"  ID= " <<setw(ndigit)<<i<<" tag= "<<fnet_tag;
    vfnet.push_back( fnet( &(*it_img),iDiam,pxlCut,fnet_tag));
    it_fnet=findFnet(&vfnet,fnet_tag);
    fnet_stack.push_back(*it_fnet);  
  }
  cout<<"  "<<nStack<<" images used for assigning beads."<<endl;            
}

/************************************************************************/
fnet3d::fnet3d(bead3d *bd3d, int *id_vec, string tag0)
{ // constructor using existing bead3d. id_vec: id range in bead3d

  int i,ndigit;
  vector<fnet>::iterator it_fnet;
  vector<bead>::iterator it_bead;
  string fnet_tag,bead_tag,sdum;
  stringstream ss;

  tag=tag0; bead3d_tag=bd3d->tag;
  nStack=id_vec[1]-id_vec[0]+1;

  bond_a=new multimap<int,int>[nStack];
  bond_b=new multimap<int,int>[nStack];

  ndigit=nDigit(nStack);
  
  for (i=id_vec[0];i<=id_vec[1];i++) {
    ss.str(""); ss.clear();
    ss<< tag<<i; fnet_tag=ss.str();
    cout<<"  ID= " <<setw(ndigit)<<i;
    vfnet.push_back( fnet( &(bd3d->bead_stack[i]), fnet_tag) );
    it_fnet=findFnet(&vfnet,fnet_tag);
    fnet_stack.push_back(*it_fnet);
  }
}

/************************************************************************/
fnet3d::fnet3d(string tag,string ifname,string bead3d_tag0)
{ /* Construct fnet from the data in a .dat file with name ifname. 
     Does not build bead. Bead should be built/loaded separately.   */
  this->tag = tag;
  if (findFnet3d(&vfnet3d,tag)==vfnet3d.end())
    error_dodri("Fnet3d tag needs to be built first."); 
  this->readData(ifname,bead3d_tag0); 
}

/************************************************************************/
// Member functions
/************************************************************************/

void fnet3d::addMissingZBond(int istk,int jstk,
			     vector<int> iSeq, vector<int> jSeq)
{ /* Add z-bonds to beads in iSeq and jSeq. 
     Sequence lists should already be ordered.
     If equal number of beads, add z-bond in sequential order.  */ 
  
  int i,minLength;
  if ((int)iSeq.size()==(int)jSeq.size()) {  // verify seq order
    for (i=0;i<(int)iSeq.size();i++)
      registerOneBondZ(iSeq[i],jSeq[i],-1,istk); 
    return;
  } //   if ((int)iSeq.size()==(int)jSeq.size()) {  
  else {
    minLength=((int)iSeq.size()<(int)jSeq.size())?((int)iSeq.size()):((int)jSeq.size());
    for (i=0;i<minLength;i++) 
      registerOneBondZ(iSeq[i],jSeq[i],-1,istk); // rough z-bonding.
    // adjust with subsequent pruning 
  }
  
  return; 
}

/************************************************************************/
void fnet3d::assignBondZ(double rCut, int flag,int *id_vec)
{ /* Assign bonds in the z-direction. For each bead on a plane, scan beads
     in nearby planes, and assign bond_a, bond_b if their xy-distance
     is less than rCut*beadDiam. First assign bond_a, then accordingly
     set bond_b for the same pair of beads in the reverse way. 
     
     if flag==1: skip adding zbonds to any beads already (from a previous 
     execution of assignbondZ, NOT this one) z-bonded.   */
  
  int ib,kb,istk,kstk,icell; // ib,kb: bead index, istk,kstk: stack index
  int ndigit=nDigit(id_vec[1]-id_vec[0]+1);
  double rCutSq,rdum;
  
  multimap<int,int>::iterator it,kt;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> 
     it_range, kt_range;

  multimap<int,int> *bond_a0,*bond_b0; // temp copies of bond_a0 and bond_b0
  bond_a0=new multimap<int,int>[nStack];
  bond_b0=new multimap<int,int>[nStack];
  for (istk=id_vec[0];istk<id_vec[1];istk++) { // to compare existing zbonds
    bond_a0[istk]=bond_a[istk]; bond_b0[istk]=bond_b[istk];
  }
  
  cout<<"  assign_bond_z: rCut= "<<rCut<<endl;
  if (flag) cout<<"    only assigning to beads without zbonds to both stacks."<<endl; 

  for (istk=id_vec[0];istk<id_vec[1];istk++) {
    kstk=istk+1; // index for the next plane
    rdum=(double)(fnet_stack[istk].beadDiam);
    rCutSq=(rdum*rCut); rCutSq=rCutSq*rCutSq;
    // Look for neighboring cells
    for (ib=0;ib<fnet_stack[istk].nBead;ib++) {
      if (flag) {
	if (bond_a0[istk].find(ib)!=bond_a0[istk].end()) continue; 
      }
      
      it=fnet_stack[istk].bead_cell.find(ib); 
      icell=it->second; // current cell number
      // Look for neighbor in the current cell
      kt_range=fnet_stack[kstk].cell_bead.equal_range(icell);
      for (kt=kt_range.first; kt!=kt_range.second; kt++) {
	kb=kt->second; 
	registerOneBondZ(ib,kb,rCutSq,istk);
      }
      // cell on the right
      if (icell%nCell!=(nCell-1)) {// not the last column
	kt_range=fnet_stack[kstk].cell_bead.equal_range(icell+1);
	for (kt=kt_range.first; kt!=kt_range.second; kt++) {
	  kb=kt->second;
	  registerOneBondZ(ib,kb,rCutSq,istk);
	}
      }
      // cell on the left
      if (icell%nCell!=0) {// not the first column
	kt_range=fnet_stack[kstk].cell_bead.equal_range(icell-1);
	for (kt=kt_range.first; kt!=kt_range.second; kt++) {
	  kb=kt->second;
	  registerOneBondZ(ib,kb,rCutSq,istk);
	}
      }
      //  below
      if (icell/nCell!=(nCell-1)) {// not the last row
	kt_range=fnet_stack[kstk].cell_bead.equal_range(icell+nCell);
	for (kt=kt_range.first; kt!=kt_range.second; kt++) {
	  kb=kt->second;
	  registerOneBondZ(ib,kb,rCutSq,istk);
	}
	// lower right
	if (icell%nCell!=(nCell-1)) {// not the last column
	  kt_range=fnet_stack[kstk].cell_bead.equal_range(icell+nCell+1);
	  for (kt=kt_range.first; kt!=kt_range.second; kt++) {
	    kb=kt->second;
	    registerOneBondZ(ib,kb,rCutSq,istk);
	  }
	}
	// lower left
	if (icell%nCell!=0) {// not the first column
	  kt_range=fnet_stack[kstk].cell_bead.equal_range(icell+nCell-1);
	  for (kt=kt_range.first; kt!=kt_range.second; kt++) {
	    kb=kt->second;
	    registerOneBondZ(ib,kb,rCutSq,istk);
	  }
	}
      } // if (icell/nCell!=(nCell-1)) {// not the last row
      // up
      if (icell>=nCell) {// not the first row
	kt_range=fnet_stack[kstk].cell_bead.equal_range(icell-nCell);
	for (kt=kt_range.first; kt!=kt_range.second; kt++) {
	  kb=kt->second;
	  registerOneBondZ(ib,kb,rCutSq,istk);
	}
	// upper right
	if (icell%nCell!=(nCell-1)) {// not the last column
	  kt_range=fnet_stack[kstk].cell_bead.equal_range(icell-nCell+1);
	  for (kt=kt_range.first; kt!=kt_range.second; kt++) {
	    kb=kt->second;
	    registerOneBondZ(ib,kb,rCutSq,istk);
	  }
	}
	// upper left
	if (icell%nCell!=0) {// not the first column
	  kt_range=fnet_stack[kstk].cell_bead.equal_range(icell-nCell-1);
	  for (kt=kt_range.first; kt!=kt_range.second; kt++) {
	    kb=kt->second;
	    registerOneBondZ(ib,kb,rCutSq,istk);
	  }
	}
      } // if (icell>=nCell) {// not the first row
    } // for (ib=0;ib<fnet_stack[istk].nBead;ib++) {

    cout<<"  "<<setw(ndigit)<<std::right<<istk<<": "<<
      bond_a[istk].size()-bond_a0[istk].size()<<" z-bonds assigned."<<endl;
  } // for (istk=0;istk<(nStack-1);istk++) {

  // clear memory
  for (istk=id_vec[0];istk<id_vec[1];istk++) { // to compare existing zbonds
    bond_a0[istk].clear(); bond_b0[istk].clear(); 
  }
  delete[] bond_a0; delete[] bond_b0;   
}

/************************************************************************/
void fnet3d::clearBondZ(int *id_vec) 
{ /* Delete all z-bonds */

  int istk; 
  for (istk=id_vec[0];istk<=id_vec[1];istk++) {
    bond_a[istk].erase(bond_a[istk].begin(),bond_a[istk].end());
    bond_b[istk].erase(bond_b[istk].begin(),bond_b[istk].end());
    assert(bond_a[istk].empty() && bond_b[istk].empty()); 
  }
  return; 
}

/************************************************************************/
void fnet3d::equalizeFilamentBond(double rCut, int *id_vec)
{
  int i,ndigit;
  ndigit=nDigit(nStack);
  for (i=id_vec[0];i<=id_vec[1];i++) {
    cout<<"  ID= "<<setw(ndigit)<<i<<" tag= "<<fnet_stack[i].tag;
    fnet_stack[i].equalizeFilamentBond(rCut);
  }
}

/************************************************************************/
void fnet3d::fillEnclosedHole(int *id_vec)
{ /* Add beads/bonds to close holes. 
     Interpolate beads if slice is skipped. 
     First test done with linear interpolation.    */
  
  int i,j,k,istk,jstk,kstk,ih,del_stk,bdIndx;
  int tot_bead,nHole,iBead,kBead;
  int ib,jb,kb,iloc=-1,iFil,jFil,kFil,filLength,*filSeq; 
  double rdum[2],rdum1[2],avgBdMass,avgDist;

  vector<int> iSeq,jSeq,*edge,*delBead;
  vector<array<int,2>> *hole;
  vector<array<double,2>> pos;
  multimap<int,int> *addBnd,*delBnd; // save bond and bead modifications
  multimap<int,int>::iterator it;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> it_range; 
  
  tot_bead=0; 
  for (istk=id_vec[0];istk<=id_vec[1];istk++) tot_bead+=fnet_stack[istk].nBead;
  hole=new vector<array<int,2>>[tot_bead];
  edge=new vector<int>[nStack]; 
  addBnd=new multimap<int,int>[nStack];
  delBnd=new multimap<int,int>[nStack];
  delBead=new vector<int>[nStack];  

  // hole[ih][0,1=bot, 2,3=top, else=edges ].at(0=stk,1=bead index):
  // there are 6 beads for "no slice" holes: the last two are the same as the "top" beads
  nHole=locateEnclosedHole(hole,id_vec); 

  /* Classify holes to determine how to fill 
     : no slice skip, one slice skip, multi slice skip  */
  for (ih=0;ih<nHole;ih++) {
    istk=hole[ih][0].at(0);
    // bottom  : vals saved this way in locateEnclosedHole()
    jstk=hole[ih][3].at(0); // top 
    del_stk=jstk-istk; // # slices in hole
    
    if (del_stk==1) {
      // just add z-bonds ; add beads (or delete) if count does not match
      iSeq.clear(); jSeq.clear();
      // list of beads in sequence at hole top/bottom      
      trackFilSeq(istk,hole[ih][0].at(1),hole[ih][1].at(1),iSeq);
      trackFilSeq(jstk,hole[ih][2].at(1),hole[ih][3].at(1),jSeq);
      // hole too small      
      if ((int)iSeq.size()<=2 && (int)jSeq.size()<=2) continue;
      // hole not enclosed       
      if ((int)iSeq.size()<=1 || (int)jSeq.size()<=1) continue; 

      // simplest case - just needs z-bonds      
      if ((int)iSeq.size()==(int)jSeq.size()) { 
	// first/last beads of *Seq should be z-bonded.
	//else, flip order (necessary for cyclic).
	assert(jstk>istk); //sanity
	if (bond_a[istk].find(iSeq.front())->second!=jSeq.front()) {
	  reverse(jSeq.begin(),jSeq.end());
	}
	iSeq.pop_back(); jSeq.pop_back(); // remove last element 
	iSeq.erase(iSeq.begin());  jSeq.erase(jSeq.begin());
	// remove first element
	addMissingZBond(istk,jstk,iSeq,jSeq); 
      }
      else { // will need to add new beads, update bond/filament/zbond
	assert(hole[ih][4].at(0)==hole[ih][5].at(0)); //sanity check
	kstk=hole[ih][4].at(0); // stack where bead will be added
	
	// only need to modify 2D bonds for pentagons 
	if (!fnet_stack[kstk].isBonded(hole[ih][4].at(1),hole[ih][5].at(1))) 
	  goto checkZBnd;
	// check/skip duplicates	
	it_range=delBnd[kstk].equal_range(hole[ih][4].at(1)); 
	for (it=it_range.first;it!=it_range.second;it++) {
	  if (it->second==hole[ih][5].at(1)) goto checkZBnd;
	}
	it_range=delBnd[kstk].equal_range(hole[ih][5].at(1));
	for (it=it_range.first;it!=it_range.second;it++) {
	  if (it->second==hole[ih][4].at(1))  goto checkZBnd; 
	}
	
	// use midpoint of "edge" beads for new bead pos
	for (j=0;j<2;j++) 
	  rdum[j]=(fnet_stack[kstk].cm[hole[ih][4].at(1)][j] +
		   fnet_stack[kstk].cm[hole[ih][5].at(1)][j])/2.;
	// use average mass for new bead mass
	avgBdMass=(fnet_stack[kstk].beadMass[hole[ih][4].at(1)]+
		   fnet_stack[kstk].beadMass[hole[ih][5].at(1)])/2.;
	bdIndx=fnet_stack[kstk].nBead;
	fnet_stack[kstk].addBead2Fn(avgBdMass,rdum[0],rdum[1]);

	// bond update: add 2 new bonds from new bead (bdIndx) to each edge,
	// delete existing bond b/en edge beads
	// update to bond object done outside of this loop
	addBnd[kstk].insert(pair<int,int>(bdIndx,hole[ih][4].at(1)));
	addBnd[kstk].insert(pair<int,int>(bdIndx,hole[ih][5].at(1)));
	delBnd[kstk].insert(pair<int,int>(hole[ih][4].at(1),
					  hole[ih][5].at(1))); 
	
      checkZBnd: // add zbonds to non-equal sequences and newly added beads 
	// list of beads in sequence at hole top/bottom
	iSeq.clear(); jSeq.clear(); 
	// seq w/ new bead	
	trackFilSeq(istk,hole[ih][0].at(1),hole[ih][1].at(1),iSeq); 
	trackFilSeq(jstk,hole[ih][2].at(1),hole[ih][3].at(1),jSeq);

	// do not adjust pent	
	if ((int)iSeq.size()!=2 && (int)jSeq.size()!=2) { 
	  assert(jstk>istk); //sanity
	  if (bond_a[istk].find(iSeq.front())->second!=jSeq.front())
	    // flip order
	    reverse(jSeq.begin(),jSeq.end());
	  iSeq.pop_back(); jSeq.pop_back(); // remove last element
	  // remove first element	  
	  iSeq.erase(iSeq.begin());  jSeq.erase(jSeq.begin());
	  addMissingZBond(istk,jstk,iSeq,jSeq); 
	}
      }
    } //     if (del_stk==1) { 
    
    else if (del_stk>1) {
      // In this pass, just add beads to missing segments
      // add based on distance of *seq beads.
      // set start pos at mdpt and move toward ends to avoid skewing placement
      // this will be fine regardless if *seq are the same size      
      iSeq.clear(); jSeq.clear();
      // list of beads in sequence at hole top/bottom      
      trackFilSeq(istk,hole[ih][0].at(1),hole[ih][1].at(1),iSeq);
      trackFilSeq(jstk,hole[ih][2].at(1),hole[ih][3].at(1),jSeq);

      /******/
      /* Get details of top/bottom boundaries: avg bond length, linear fit */
      avgDist=0; // calc avg bond length
      avgBdMass=fnet_stack[istk].beadMass[iSeq[0]]; // calc avg bead mass
      for (i=1;i<(int)iSeq.size();i++) {
	for (k=0;k<2;k++) {
	  rdum[k]=fnet_stack[istk].cm[iSeq[i-1]][k];
	  rdum1[k]=fnet_stack[istk].cm[iSeq[i]][k];
	}
	avgDist+=getDistPos(rdum,rdum1,0);
	avgBdMass+=fnet_stack[istk].beadMass[iSeq[i]]; 
      }
      avgBdMass+=fnet_stack[jstk].beadMass[jSeq[0]]; // calc avg bead mass      
      for (i=1;i<(int)jSeq.size();i++) {
	for (k=0;k<2;k++) {
	  rdum[k]=fnet_stack[jstk].cm[jSeq[i-1]][k];
	  rdum1[k]=fnet_stack[jstk].cm[jSeq[i]][k];
	}
	avgDist+=getDistPos(rdum,rdum1,0);
	avgBdMass+=fnet_stack[jstk].beadMass[jSeq[i]]; 	
      }
      // no. bonds= sequence size -1      
      avgDist/=(((int)iSeq.size()-1)+((int)jSeq.size()-1)); 
      avgBdMass/=((double)iSeq.size()+(double)jSeq.size());
      /*****/
      
      // add beads in b/en ends from hole[ih][4-end].at(1) per matching stack
      for (kstk=id_vec[0];kstk<id_vec[1];kstk++) edge[kstk].clear(); 
      for (i=4;i<(int)hole[ih].size();i++) {
      	kstk=hole[ih][i].at(0);
      	edge[kstk].push_back(hole[ih][i].at(1));  // edge[stk][bd];
      }
      for (kstk=istk+1;kstk<jstk;kstk++) { // check all slices
      	pos.clear(); 
      	assert(edge[kstk].size()==2); // sanity check
	//	flag=1; // change flag==0 if hole is not proper

	// check for "improper" hole assignments: 
	// two edge beads are bonded - do nothing, skip adding beads/bonds 
	if (fnet_stack[kstk].isBonded(edge[kstk][0],edge[kstk][1])) continue; 
	else { // check beads are not at the ends of a filament 
	  //  for now, sever "dangling" part of filament, bond the two edge beads (no new bead)
	  // later, implement a smoothen to smoothen kstk based on hole top/bot
	  ib=edge[kstk][0]; jb=edge[kstk][1]; 
	  iFil=fnet_stack[kstk].bead_fil.find(ib)->second;
	  jFil=fnet_stack[kstk].bead_fil.find(jb)->second;

	  // allow if both beads are not part of a filament 	  
	  if (iFil==0 && jFil==0) continue;  
	  
	  filSeq=new int[fnet_stack[kstk].maxFilLength]; 	  
	  // check iff just one or both edges need to be severed
	  for (k=0;k<2;k++) { 
	    if (k==0) {kFil=iFil;kb=ib;}
	    else {kFil=jFil; ; kb=jb;}
	    filLength=fnet_stack[kstk].findFilamentSeq(kFil,filSeq);
	    filLength=abs(filLength); 
	    if (kb!=filSeq[0] && kb!=filSeq[filLength-1]) {
	      for (i=0;i<filLength;i++) {
		if (filSeq[i]==kb) { iloc=i; break; }
	      }
		// delete from iloc to end	      
	      if (iloc>floor((double)filLength/2.)) {
		for (i=iloc;i<filLength;i++)  // save beads for deletion
		  delBead[kstk].push_back(filSeq[i]); 
	      }
	      else { // delete from 0 to iloc
		for (i=0;i<iloc;i++)  // save beads for deletion
		  delBead[kstk].push_back(filSeq[i]); 
	      }
	    } // 	  if (ib!=filSeq[0] || ib!=filSeq[nVec-1]) {
	  } // pair<int,int>(kb+1,kb0)
	  delete[] filSeq; 	  
	}	
      }
    }
  } //    for (ih=0;ih<nHole;ih++) {
  
  /*****  Add/remove bonds, update filament *****/ 
  for (istk=id_vec[0];istk<=id_vec[1];istk++) {
    // delete bonds between ends if necessary
    for (it=delBnd[istk].begin();it!=delBnd[istk].end();++it) {
      iBead=it->first; kBead=it->second;       
      if (!(fnet_stack[istk].isBonded(iBead,kBead))) continue;
      // bond already deleted
      fnet_stack[istk].removeOneBond(iBead,kBead);
      fnet_stack[istk].updateFilament(iBead,kBead); // filament should be split
      assert(!(fnet_stack[istk].isBonded(iBead,kBead))); //sanity check
    }
    
    // add bonds b/en new beads
    for (it=addBnd[istk].begin(); it!= addBnd[istk].end();++it) {
      fnet_stack[istk].registerOneBond(it->first,it->second,-1,-1,
				       &(*findImg(fnet_stack[istk].img_tag)));
      // temporarily add any new beads to filament=++nFil
      if (fnet_stack[istk].bead_fil.find(it->first)==
  	  fnet_stack[istk].bead_fil.end()) {
  	fnet_stack[istk].bead_fil.insert(pair<int,int>(it->first,0));
  	fnet_stack[istk].fil_bead.insert(pair<int,int>(0,it->first)); 	
      }
      if (fnet_stack[istk].bead_fil.find(it->second)==
  	  fnet_stack[istk].bead_fil.end()) {
  	fnet_stack[istk].bead_fil.insert(pair<int,int>(it->second,0));
  	fnet_stack[istk].fil_bead.insert(pair<int,int>(0,it->second)); 	
      }
      
      // update filament
      // iBead and kBead should be in the same fil      
      iBead=it->first; kBead=it->second; 
      fnet_stack[istk].updateFilament(iBead,kBead);
    }
  }  
  return; 
}

/************************************************************************/
void fnet3d::findFilament(int flag,int *id_vec)
{
  int i,ndigit;
  ndigit=nDigit(nStack);
  for (i=id_vec[0];i<=id_vec[1];i++) {
    cout<<"  ID= "<<setw(ndigit)<<i<<" tag= "<<fnet_stack[i].tag;
    fnet_stack[i].findFilament(flag);
  }
}

/************************************************************************/
void fnet3d::findNextZ(int istk,int iBead, int iFil,int kstk,
		       vector<int> &z_bonded)
{ /* Search filament for flanking z-bonds from reference iBead. 
     Returns bead z-bonded to kstk after "empty" region.   */

  if (!iFil) {
    cout<<" WARNING: "<< iBead<<" does not belong to a filament."<<endl; 
    return;
  }
  
  int i,iloc,jloc,jBead;
  int cyclic,filLength, *filSeq=new int[fnet_stack[istk].maxFilLength];
  multimap<int,int> temp_bond; 

  // determine whether to use bond_a or bond_b for zbond check 
  if (kstk>istk) temp_bond=bond_a[istk];
  else temp_bond=bond_b[istk]; 

  filLength=fnet_stack[istk].findFilamentSeq(iFil,filSeq);
  cyclic=(filLength<0)?1:0;  // record if filament is cyclic
  
  filLength=abs(filLength); // abs value of filLength
  iloc=-1; 
  for (i=0;i<filLength;i++) { //locate iBead position in filament
    if (filSeq[i]==iBead) { iloc=i; break; }
  }
  if (iloc<0)
    error_dodri("  fnet3d::findNextZ: invalid position in filament.");
  
  if (!cyclic) { // take care of ends
    if (iloc<filLength) {   // check to the right (+iloc) of iBead
      for (i=iloc+1;i<filLength;i++) {
	jBead=filSeq[i];
	if (temp_bond.find(jBead)!=temp_bond.end()) {
	  if (temp_bond.find(filSeq[i-1])==temp_bond.end())
	    z_bonded.push_back(jBead);
	  break; // end search 
	}
      }
    }
    if (iloc-1>=0) {     // check to the left (-iloc) of iBead 
      for (i=iloc-1;i>=0;i--) {
	jBead=filSeq[i];
	if (temp_bond.find(jBead)!=temp_bond.end()) {
	  if (temp_bond.find(filSeq[i+1])==temp_bond.end())
	    z_bonded.push_back(jBead);
	  break; // end search
	}
      }
    } //   if (iloc-2>=0) {
  } //   if (!cyclic) { // take care of ends
  else { //cyclic bounds
    for (i=1;i<filLength;i++) { // check to the right (+i)
      jloc=iloc+i;
      if (jloc>=filLength) jloc=jloc-filLength; //i-(jloc-iloc); // position for cyclic
      jBead=filSeq[jloc];
      if (temp_bond.find(jBead)!=temp_bond.end()) {
	z_bonded.push_back(jBead);
	break; // end search 
      }
    }
    for (i=1;i<filLength;i++) { // check to the left (-i)
      jloc=iloc-i;
      if (jloc<0) jloc=filLength+jloc;
      jBead=filSeq[jloc];
      if (temp_bond.find(jBead)!=temp_bond.end()) {
	z_bonded.push_back(jBead);
	break; // end search 
      }
    }
  }
  
  delete[] filSeq;   
  return; 
}

/************************************************************************/
void fnet3d::findSeqZ(int iBead0,int istk,multimap<int,int> &zlist,
		      int *id_vec)
{ /* Analogous to fnet::findFilament() but for z-bonds. Saves 
     z-bond sequence from ibead. Checks prior and next slices. 
     Assume each bead is z-bonded once (ie no triangles).   */
  
  int jstk,iBead,jBead;
  
  assert(istk<=id_vec[1] && istk>=id_vec[0]); // sanity check
  zlist.clear();
  
  iBead=iBead0;       // z-bond to next slice; stop if does not exist
  for (jstk=istk;jstk<id_vec[1];jstk++) {
    if (bond_a[jstk].find(iBead)==bond_a[jstk].end()) break; 
    jBead=bond_a[jstk].find(iBead)->second;
    zlist.insert(pair<int,int>(jstk+1,jBead));
    iBead=jBead;       
  }
  
  iBead=iBead0;
  for (jstk=istk;jstk>=id_vec[0];jstk--) { // from istk to id_vec[0]
    if (bond_b[jstk].find(iBead)==bond_b[jstk].end()) break; 
    jBead=bond_b[jstk].find(iBead)->second;
    zlist.insert(pair<int,int>(jstk-1,jBead));
    iBead=jBead;   
  }
  
  return; 
}

/************************************************************************/
int fnet3d::fnetTag2Id(string sdum)
{ /* Find id number for the fnet with tag sdum.
     Returns: id number or -1 if tag not found */
  int i,j;
  vector<fnet>::iterator it_fnet;
  j=-1;
  for (i=0;i<nStack;i++) {
    if (fnet_stack[i].tag==sdum) {j=i; break;}
  }
  return j;
}

/************************************************************************/
int fnet3d::getBondVecZ(int p, int istk, int iSkip, double **bondVec,
			int *bondVecList, int dir)
{ /* Get z-dir bond vectors relative to bead p.
     Save bond vectors to bondVec, bead indexes to bondVecList
     dir>=0: use bond_a, dir<0: use bond_b
     In both dir's, bondVec points away from p.
     Do not include bead iSkip in considering bond vector (iSkip=-1 for no skip)
     Returns: number of bond vectors  */
  
  int q, n=0, idir=0;
  double r0,r1;
  multimap<int,int>::iterator it;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> it_range;  
  
  if (dir<0) {
    idir=-1; it_range=bond_b[istk].equal_range(p);
  }
  else {
    idir=1;
    it_range=bond_a[istk].equal_range(p);
  }
  
  if (istk>=nStack) {
    cout<<"  Warning: Went beyond the last stack. No z-bond assigned."<<endl;
    return 0;
  }
  for (it=it_range.first;it!=it_range.second;it++) {
    q=it->second; // p-q form a bond
    if (q!=iSkip) {
      r0=fnet_stack[istk+idir].cm[q][0]-fnet_stack[istk].cm[p][0];     
      r1=fnet_stack[istk+idir].cm[q][1]-fnet_stack[istk].cm[p][1];
      bondVec[n][0]=r0;
      bondVec[n][1]=r1;
      bondVec[n][2]=fnet_stack[istk+idir].zCoord-fnet_stack[istk].zCoord;
      bondVecList[n]=q;
      ++n; assert(n<maxVec);
    }
  }  
  return n;
}

/************************************************************************/
void fnet3d::joinCyclicFil(double rCut, int *id_vec)
{
  int i,ndigit;
  ndigit=nDigit(nStack);
  for (i=id_vec[0];i<=id_vec[1];i++) {
    cout<<"  ID= "<<setw(ndigit)<<i<<" tag= "<<fnet_stack[i].tag;
    fnet_stack[i].joinCyclicFil(rCut);
  }
}

/************************************************************************/
int fnet3d::locateEnclosedHole(vector<array<int,2>> *hole,int *id_vec)
{ /* Locate bounded surface holes. This method just uses fnet3d BOC.  */ 
  
  int i,j,k,istk,jstk,nHole,flag,ifil;
  int ib,jb,ib_zba,jb_zba;
  double top_match[2]; 

  array<int,2> addarr,temparr;
  vector<int> add,*seed, *delBot; // seed beads for holes
  vector<array<int,2>> *bot,*top, edge, cpBot; // beads in presumptive hole outline 
  array<double,2> matcharr;
  vector<array<double,2>> match; // to match top/bottom: [stk][index]
  
  bot=new vector<array<int,2>>[nStack];
  top=new vector<array<int,2>>[nStack];    
  seed=new vector<int>[nStack];
  delBot=new vector<int>[nStack];  

  //  hole[ihole].at(0)==stack; .at(1)==ibead
  // top[ihole (presumptive}.at(0)==nstack; .at(1)==ibead
  // bot[ihole (presumptive}.at(0)==nstack; .at(1)==ibead
  
  // seed beads: beads w/ 2 2d bonds and onely one z-bond
  for (istk=id_vec[0];istk<=id_vec[1];istk++) { 
    seed[istk].clear(); 
    for (ib=0;ib<fnet_stack[istk].nBead;ib++) {
      if (fnet_stack[istk].bond.count(ib)==2) {
	// locate corners of presumptive hole "top"
	if (istk==id_vec[0]) goto checkBot; // first slice can't be a "top"
	if (bond_b[istk].find(ib)==bond_b[istk].end()) {
	  add.clear();
	  ifil=fnet_stack[istk].bead_fil.find(ib)->second; 
	  findNextZ(istk,ib,ifil,istk-1,add);
	  if ((int)add.size()!=2) goto checkBot ; // to save ENCLOSED holes
	  flag=0;
	  for (i=0;i<(int)top[istk].size();i++) { // to avoid duplicates 
	    if ( (top[istk][i].at(0)==add[0] && top[istk][i].at(1)==add[1]) ||
		 (top[istk][i].at(0)==add[1] && top[istk][i].at(1)==add[0]) ) 
	      flag=1; // this pairing already exists
	  }	  
	  if (!flag) { // else, add to top vector
	    addarr.at(0)=add[0]; addarr.at(1)=add[1];
	    top[istk].push_back(addarr); 
	  }
	}
      checkBot: 	// locate corners of presumptive hole "bottom"
	if (istk==id_vec[1]) continue; // last slice can't be a "bottom"
	if (bond_a[istk].find(ib)==bond_a[istk].end() ) { //&&
	  add.clear();
	  ifil=fnet_stack[istk].bead_fil.find(ib)->second; 	  	  
	  findNextZ(istk,ib,ifil,istk+1,add);
	  if ((int)add.size()!=2) continue; // to save ENCLOSED holes
	  flag=0;
	  for (i=0;i<(int)bot[istk].size();i++) { //  to avoid duplicates
	    if ( (bot[istk][i].at(0)==add[0] && bot[istk][i].at(1)==add[1]) ||
		 (bot[istk][i].at(0)==add[1] && bot[istk][i].at(1)==add[0]) ) 
	      flag=1; // this pairing already exists
	  }
	  if (!flag) { // else, add to bot vector
	    addarr.at(0)=add[0]; addarr.at(1)=add[1];
	    bot[istk].push_back(addarr); 
	  }
	}
      } //       if (fnet_stack[istk].bond.count(ib)==2) {
    }
  } //   for (istk=id_vec[0];istk<id_vec[1];istk++) {
  
  // check for no-slice skip pentagons (most common "hole")
  nHole=0;
  for (istk=id_vec[0];istk<id_vec[1];istk++) { 
    for (i=0;i<(int)bot[istk].size();i++) { // for each presumptive hole bottom 
      ib=bot[istk][i].at(0);  jb=bot[istk][i].at(1);
      if (bond_a[istk].find(ib)!=bond_a[istk].end() &&
  	  bond_a[istk].find(jb)!=bond_a[istk].end()) {
  	ib_zba=bond_a[istk].find(ib)->second;	    
  	jb_zba=bond_a[istk].find(jb)->second;
  	jstk=istk+1;	// *_zba beads are in istk+1;
  	if (fnet_stack[jstk].isBonded(ib_zba,jb_zba)) { // pentagon hole
  	  // record hole outline ("top" will be the same as edge here)
  	  temparr.at(0)=istk; temparr.at(1)=bot[istk][i].at(0); // pt 1
  	  hole[nHole].push_back(temparr);	
  	  temparr.at(0)=istk; temparr.at(1)=bot[istk][i].at(1); // pt 2
  	  hole[nHole].push_back(temparr);
  	  temparr.at(0)=jstk; temparr.at(1)=ib_zba; // pt 3
  	  hole[nHole].push_back(temparr);		
  	  temparr.at(0)=jstk; temparr.at(1)=jb_zba; // pt 4
  	  hole[nHole].push_back(temparr);
  	  // record edge where new bead is needed
  	  // to fit with hole filling, add edge beads to hole[istk][4 and 5]
  	  temparr.at(0)=jstk; temparr.at(1)=ib_zba; 
  	  hole[nHole].push_back(temparr); 
  	  temparr.at(0)=jstk; temparr.at(1)=jb_zba; 
  	  hole[nHole].push_back(temparr);
  	  // incr number of holes 
  	  nHole=nHole+1;
  	  //save curr bot combo for deletion to avoid being recorded again in general search	  
  	  delBot[istk].push_back(i); 
  	} // 	if (fnet_stack[jstk].isBonded(ib_zba,jb_zba)) { // pentagon hole
      }
    } //  for (i=0;i<(int)bot[istk].size();i++) { // for each presumptive hole bottom 
    // do the same from any top bead but checking stack below
    for (i=0;i<(int)top[istk].size();i++) { // for each presumptive hole bottom 
      ib=top[istk][i].at(0);  jb=top[istk][i].at(1);
      if (bond_b[istk].find(ib)!=bond_b[istk].end() &&
  	  bond_b[istk].find(jb)!=bond_b[istk].end()) {
  	ib_zba=bond_b[istk].find(ib)->second;	    
  	jb_zba=bond_b[istk].find(jb)->second;
  	jstk=istk-1;	// *_zba beads are in istk-1;	
    	if (fnet_stack[jstk].isBonded(ib_zba,jb_zba)) { // pentagon hole
  	  // record hole outline ("bottom" will be the same as edge here)
  	  temparr.at(0)=jstk; temparr.at(1)=ib_zba; // pt 1
  	  hole[nHole].push_back(temparr);	
  	  temparr.at(0)=jstk; temparr.at(1)=jb_zba; // pt 2
  	  hole[nHole].push_back(temparr);	
  	  temparr.at(0)=istk; temparr.at(1)=top[istk][i].at(0); // pt 3
  	  hole[nHole].push_back(temparr);	
  	  temparr.at(0)=istk; temparr.at(1)=top[istk][i].at(1); // pt 4
  	  hole[nHole].push_back(temparr);
  	  // record edge where new bead is needed
  	  temparr.at(0)=jstk; temparr.at(1)=ib_zba; 
  	  hole[nHole].push_back(temparr); 
  	  temparr.at(0)=jstk; temparr.at(1)=jb_zba; 
  	  hole[nHole].push_back(temparr);
  	  // incr number of holes
  	  nHole=nHole+1;
  	}
      }
    }
  }
  
  // update bot vector to remove bead combos found as pentagon holes.
  // no need to do this to top vector (for now), since search is bottom-up  
  for (istk=id_vec[0];istk<id_vec[1];istk++) {
    cpBot.clear();
    for (i=0;i<(int)bot[istk].size();i++) {
      if (find(delBot[istk].begin(),delBot[istk].end(),i)==delBot[istk].end()) {
	// current bot combo was not marked for deletion, save to temp bot 
	temparr.at(0)=bot[istk][i].at(0);  temparr.at(1)=bot[istk][i].at(1);
	cpBot.push_back(temparr); 
      }
    } //     for (i=0;i<(int)bot[istk].size();i++) {
    bot[istk].clear(); // remove all entries from bot[istk]
    bot[istk]=cpBot; // replace with cpBot
  } //   for (istk=id_vec[0];istk<id_vec[1];istk++) { 
  
  /**  search for non-pentagon holes **/ 
  // nHole=0; 
  for (istk=id_vec[0];istk<id_vec[1];istk++) {
    for (i=0;i<(int)bot[istk].size();i++) { // for each presumptive hole bottom             
      match.clear();  // to save matches
      edge.clear(); // to save presumptive edges during search 
      ib=bot[istk][i].at(0);  // search one bot end to find a match in top end
      
      //// save vec from .at(0) to .at(1) to track hole "direction" 
      //ibnd[0]=fnet_stack[istk].cm[bot[istk][i].at(1)][0] - fnet_stack[istk].cm[ib][0];
      //ibnd[1]=fnet_stack[istk].cm[bot[istk][i].at(1)][1] - fnet_stack[istk].cm[ib][1];  
      
      flag=0; 
      jb=ib; jstk=istk; // to start search
      while (!flag && jstk+1<=id_vec[1]) {
	// update jb for next iteration of while loop
	if (bond_a[jstk].find(jb)!=bond_a[jstk].end()) {
	  jb_zba=bond_a[jstk].find(jb)->second;	    
	  jb=jb_zba; jstk=jstk+1;	  
	  temparr.at(0)=jstk; temparr.at(1)=jb; 
	  edge.push_back(temparr); 
	}
	else break; // for now jstk=jstk+1; 
	
	for (j=0;j<(int)top[jstk].size();j++) { // check if this bead is in top
	  if (top[jstk][j].at(0)==jb || top[jstk][j].at(1)==jb) {
	    flag=1; //match[0][0]=jstk; match[0][1]=j;
	    matcharr.at(0)=jstk; matcharr.at(1)=j;
	    match.push_back(matcharr); // keep track of ALL matches
	  }
	}
      } //       while (!flag && jstk+1<=id_vec[1]) {
      
      if (match.empty()) continue; // no need to chck other bot end
      ib=bot[istk][i].at(1);

      flag=0; 
      jb=ib; jstk=istk; // to start search
      while (!flag && jstk+1<=id_vec[1]) {
	// update jb for next iteration of while loop
	if (bond_a[jstk].find(jb)!=bond_a[jstk].end()) {
	  jb_zba=bond_a[jstk].find(jb)->second;	    
	  jb=jb_zba; jstk=jstk+1;
	  temparr.at(0)=jstk; temparr.at(1)=jb; 
	  edge.push_back(temparr); 
	}
	else break; //for now jstk=jstk+1; 
	for (j=0;j<(int)top[jstk].size();j++) { // check if this bead is in top
	  if (top[jstk][j].at(0)==jb || top[jstk][j].at(1)==jb) {
	    for (k=0;k<(int)match.size();k++) {
	      if (match[k][0]==jstk && match[k][1]==j) {
		// this top matches that found in "bot" search
		flag=1; top_match[0]=jstk; top_match[1]=j; 
		goto insert; 
	      }
	    }
	  }
	}
      } //       while (!flag && jstk+1<=id_vec[1]) {
      
    insert: 
      // top index matches for both bot beads
      if (flag) {	// save top and bottom 
	temparr.at(0)=istk; temparr.at(1)=bot[istk][i].at(0); // pt 1
	hole[nHole].push_back(temparr);	
	temparr.at(0)=istk; temparr.at(1)=bot[istk][i].at(1); // pt 2
	hole[nHole].push_back(temparr);
	jstk=top_match[0]; jb=top_match[1]; 
	temparr.at(0)=jstk; temparr.at(1)=top[jstk][jb].at(0); // pt 3
	hole[nHole].push_back(temparr);		
	temparr.at(0)=jstk; temparr.at(1)=top[jstk][jb].at(1); // pt 4
	hole[nHole].push_back(temparr);	
	
	for (k=0;k<(int)edge.size();k++) { // save edges to hole
	  temparr.at(0)=edge[k][0]; temparr.at(1)=edge[k][1];
	  hole[nHole].push_back(temparr);
	}  	
	nHole=nHole+1;
      }
    } //     for (i=0;i<(int)bot[istk].size();i++) { 
  } //    for (istk=id_vec[0];istk<id_vec[1];istk++) { // for "bottom" of hole
  
  return nHole;     
}

/************************************************************************/
void fnet3d::measureWidth3D(int N, double arr[], int nslice,string ofname,
			    int *id_vec,vector<array<double,2>> *axis_pos,
			    vector<double> *axis_width,
			    vector<array<double,2>> *end_pos0,
			    vector<array<double,2>> *end_pos1,
			    string ax0)
{ /* Measure the width of an image based on making a plane from 3 points. 
     Call fitPlane() to calc plane equation.  

    --> arr[] will be found here if ==-1 based on boc orientation. 
    ax0:: should match w/ orientation axis used with bead3d::orient   */
  
  int istk,i,j,k,istart,iend,flag,ndigit;
  double maxZ,minZ,dz,zdum;
  double *pvec=new double[3];
  double **plpts=new double*[nStack]; //bounds for plane points for each z
  array<double,3> ptarr;
  vector<array<double,3>> pttemp;
  vector<int> range; // non-empty slices 

  // for locating reference values
  int ib,totBead,imin,imax,zStk0,zStk1=-1;
  double min,max,sig;
  double *bm=new double[2],**pos, *bd_mass;

  // arbitrary scaling values: necessary when plane is auto-located
  double delt=10., a0=delt/2.;
  
  if (arr[0]!=-1) goto measure; 

  /*********************/
  /* LOCATE reference points: boc centroid, max&min of x or y */
  // beads should prefereably be oriented along either x or y axis 
  cout<<"  MESSAGE:: 3 reference points will be automatically located."<<endl;

  // get total bead count for pos sizing
  totBead=0; zStk0=id_vec[1]; 
  for (istk=id_vec[0];istk<=id_vec[1];istk++) {
    totBead+=fnet_stack[istk].nBead; 
    // get range of slices that contain beads 
    if (fnet_stack[istk].nBead>0 && istk > zStk0) zStk1=istk;
    else if (fnet_stack[istk].nBead>0) zStk0=istk; 
  }

  pos=new double*[2];
  for (j=0;j<2;j++) pos[j]=new double[totBead];
  bd_mass=new double[totBead]; 

  for (i=0;i<totBead;i++) bd_mass[i]=1.; 
  
  ib=0; 
  for (istk=id_vec[0];istk<=id_vec[1];istk++) {
    if (!fnet_stack[istk].nBead) continue; // skip if no beads    
    for (i=0;i<fnet_stack[istk].nBead;i++) {
      for (j=0;j<2;j++) pos[j][ib]=fnet_stack[istk].cm[i][j];
      ++ib;
    }
  }
  for (j=0;j<2;j++)  {
    getavgw(&bm[j],&sig,pos[j],bd_mass,totBead);  // centroid
    cout<<"getavgw:"<<bm[0]<<" " <<bm[1]<<endl; 
  }

  
  if (ax0=="x") MinMax(totBead,pos[0],&imin,&imax,&min,&max); 
  else if (ax0=="y")  MinMax(totBead,pos[1],&imin,&imax,&min,&max);
  else error_dodri("fnet3d::measureWidth3D: invalid axis."); 
  
  // save min/max pos to arr (which is linear)
  // i{min,max}Stk: to get bead pos at imin and imax
  // zStk{0,1}: zcoords of first/last stack with beads  
  // arr[] order: min, bcm, max 
  
  if (ax0=="x") {
    arr[0]=pos[0][imin];
    arr[1]=bm[1]-delt; // bmy with offset    
    arr[6]=pos[0][imax];
    arr[7]=bm[1]+delt; // bmy with offset    
  }
  else if (ax0=="y") {
    arr[0]=bm[0]-delt; // bmx with offset
    arr[1]=pos[1][imin];    
    arr[6]=bm[0]+delt; // bmx with offset 
    arr[7]=pos[1][imax];
  }
  arr[2]=fnet_stack[zStk0].zCoord;
  arr[8]=fnet_stack[zStk1].zCoord;    

  // 3rd point: centroid 
  arr[3]=bm[0]; arr[4]=bm[1];
  arr[5]= (arr[2]+arr[8])/2. ;
  
  /*********************/
  
 measure:
  
  if (ofname=="") ofname="w_dist";   
  cout<<"  measure_width: npts= "<<N/3<<" nslice= "<<nslice<<endl;
  cout<<"    reference points (x y z)=";
  for (i=0;i<N;i++) {
    if (i%3==0)   cout <<" #"<<i/3+1<<": ";
    cout<< arr[i]<<" ";
  }
  cout<<endl;

  //ensure all output vectors are clear
  for (istk=id_vec[0];istk<=id_vec[1];istk++) {       
    axis_pos[istk].clear();
    axis_width[istk].clear();
    end_pos0[istk].clear();
    end_pos1[istk].clear();
  }
  
  //calc plane equation format: z=pvec[0]*x+pvec[1]*y+pvec[2]
  fitPlane(N,arr,pvec);
  
  //calc plane vals & save to plpts[nStack][5];
  //4: (x0,y0 and x1,y1) *top and bottom bounds
  for (i=0;i<nStack;i++) plpts[i]=new double[4];
  
  //take care of different setz dz and origins
  for (istk=id_vec[0];istk<=id_vec[1];istk++) {
    if (fnet_stack[istk].nBead==0) continue;
    range.push_back(istk); //save istk of non-empty slices
  }
  
  istart=(range.front()>range.back())?(range.back()):(range.front()); 
  iend=(range.front()>range.back())?(range.front()):(range.back());
 
  dz=abs(fnet_stack[istart].zCoord-fnet_stack[istart+1].zCoord);
  minZ=(fnet_stack[range.front()].zCoord<fnet_stack[range.back()].zCoord)?
    (fnet_stack[range.front()].zCoord):(fnet_stack[range.back()].zCoord); 
  maxZ=(fnet_stack[range.front()].zCoord>fnet_stack[range.back()].zCoord)?
    (fnet_stack[range.front()].zCoord):(fnet_stack[range.back()].zCoord); 

  for (i=0;i<fnet_stack[range.front()].rows;i++) {
    for (j=0;j<fnet_stack[range.front()].rows;j++) {
      zdum=pvec[0]*i+pvec[1]*j+pvec[2];
      if (zdum>minZ-dz && zdum<maxZ+dz) {
  	ptarr.at(0)=i;ptarr.at(1)=j;ptarr.at(2)=zdum;
  	pttemp.push_back(ptarr); //save i and j as well
      }
    }
  } //   for (i=0;i<fnet_stack[range.front()].rows;i++) { 

  for (i=istart;i<=iend;i++) { // initialize
    plpts[i][0]=plpts[i][1]=plpts[i][2]=plpts[i][3]=-1;
  }
  
  for (i=istart;i<=iend;i++) {
    zdum=fnet_stack[i].zCoord;       flag=0;
    for (k=0;k<(int)pttemp.size();k++) {
      //save first and last vals in range
      if ( (pttemp[k][2] > zdum-dz*a0 && pttemp[k][2] < zdum+dz*a0) ||
  	   (pttemp[k][2] < zdum-dz*a0 && pttemp[k][2] > zdum+dz*a0) )  {
  	//save only first and last vals in set
  	if (!flag) { 	  flag=1;
  	  plpts[i][0]=pttemp[k][0]; // first x
  	  plpts[i][1]=pttemp[k][1]; //first y
  	}
  	plpts[i][2]=pttemp[k][0];  //last x
  	plpts[i][3]=pttemp[k][1];  //last y
      }
    }
  }
  
  ndigit=nDigit(nStack);
  for (i=istart;i<=iend;i++) {
    cout<<"  ID= "<<setw(ndigit)<<i<<" tag= "<<fnet_stack[i].tag ;
    if (plpts[i][0]==-1)
      error_dodri("  fnet3d::measureWidth: invalid reference value");
    fnet_stack[i].measureWidth(plpts[i],nslice,ofname+to_string(i),
			       axis_pos[i],axis_width[i],
			       end_pos0[i],end_pos1[i],ax0,a0);
  }
  
  for (i=0;i<nStack;i++) delete[] plpts[i];
  delete pvec; delete[] plpts;
}

/************************************************************************/
void fnet3d::pruneBondZ(int *id_vec)
{ /* Removes multiple bondZ's from the same bead */
  
  int ib, istk, nRemoved, pdum0,pdum1;
  int ndigit=nDigit(id_vec[1]-id_vec[0]+1);
  pair<multimap<int, int>::iterator, multimap<int, int>::iterator> it_range;
  multimap<int, int>::iterator it;
  for (istk=id_vec[0];istk<=id_vec[1];istk++) {
    nRemoved = 0;
    for (ib = 0; ib<fnet_stack[istk].nBead; ib++) {
      pdum0=pdum1=0;
      if (istk< id_vec[1]) pdum0 = pruneBondSubZ(ib, istk,1);
      if (istk>0) pdum1 = pruneBondSubZ(ib, istk,-1);
      nRemoved += pdum0+pdum1;
    }
    cout<<"  ID= "<<setw(ndigit)<<istk<<" tag= "<<fnet_stack[istk].tag ; 
    cout<<"  "<<nRemoved<<" z-bonds removed"<<endl;
  }
}

/************************************************************************/
int fnet3d::pruneBondSubZ(int ib, int istk, int dir)
{  /* For multiple bondZ for bead ib on plane istk, keep only the shortest
   bondZ and delete the rest. 
   dir>=0: use bond_a, dir<0: use bond_b
   Two dirs are needed since bond_a and bond_b may separately have 
   multiple bonds.

   returns number of bonds removed. 
*/
  
  int iv,imin;
  int  nVec, *bondVecList=new int[maxVec];
  double dr,dr_min, **bondVec=new double*[maxVec];
  multimap<int,int>::iterator jt;
  for (iv=0;iv<maxVec;iv++) bondVec[iv]=new double[3];
  
  nVec=getBondVecZ(ib,istk,-1,bondVec,bondVecList,dir);
  if (nVec<2) return 0 ; // no need to check a single bond

  imin=0; // index for min bond distance
  dr_min=dotprod(bondVec[0],bondVec[0],3);
  
  for (iv=1;iv<nVec;iv++) { // find minimum distance bond
    dr=dotprod(bondVec[iv],bondVec[iv],3);
    if (dr<dr_min) { imin=iv; dr_min=dr; }
  }
  for (iv=0;iv<nVec;iv++) {
    if (iv!=imin) {
      if (dir>=0) removeOneBondZ(ib,bondVecList[iv],istk);
      else removeOneBondZ(bondVecList[iv],ib,istk-1);
    }
  }
  for (iv = 0; iv < maxVec; iv++) delete[] bondVec[iv];
  delete[] bondVec;  delete[] bondVecList;
  return nVec-1;
}

/************************************************************************/
void fnet3d::readData(string fname,string bd3d_tag0)
{ /* Read a .dat file made from fnet3d::writeData containing 
     bond (2D and 3D) info. Beads should be created separately. */
    
  int i, j, ivec[1], istk=0, stk=0, ndigit, nStack0=0;
  int bnd_temp[8], iBead, jBead;
  
  stringstream ss;
  string bd3d_tag,bd_tag,line,ntag,ofname;
  ifstream fin;
  vector<fnet>::iterator it_fnet; 
  bool b3d_tag, nstk, flag;
  vector<bead3d>::iterator it_bead3d; 
  flag = b3d_tag = nstk = false;

  multimap<int,int>::iterator it;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> it_range;                              
  ofname=addFileNameExt(fname,".dat"); 
  cout<<"  read:: Read data from: " << ofname<<endl;
  fin.open(ofname.c_str(),ios::in);
  if (fin.bad()) error_dodri("dodri: Invalid fnet3d file"); 
  
  // read file header
  while (!flag && fin.good()) {
    flag = b3d_tag && nstk;
    if (flag) break;
    
    if (line.find("Stack: ")!=std::string::npos)
      error_dodri("READ: Incomplete fnet3d header");
    getline(fin,line);
    
    if (line=="") continue;
    ss.str(line); 
    
    //read bead3d tag
    if (!b3d_tag && line.find("bead3d_tag=")!=string::npos) {
      if (bd3d_tag0=="")
	bd3d_tag=getCommOptString(ss,"bead3d_tag=","");
      else bd3d_tag=bd3d_tag0;
      if (bd3d_tag!="") {
	bead3d_tag = bd3d_tag;
	b3d_tag = true;
	continue;
      }
      else error_dodri("READ: missing bead3d tag");
    }
    
    //read nstack
    if (!nstk && line.find("nStack=")!=string::npos) {
      getCommOptInt(ss, "nStack=",ivec,1,-1);
      if (ivec[0] != -1) {
	nStack0 = ivec[0];
	nstk = true;
	continue;
      }
      else error_dodri("READ: invalid or missing nStack");
    }
  } //  while (!flag && fin.good()) {

  assert(nStack0==nStack); //sanity check
  ndigit = nDigit(nStack);
  it_bead3d=findBead3d(&vbead3d,bead3d_tag);
  if (it_bead3d==vbead3d.end()) error_dodri("READ: bead3d tag not found");
  
  while (fin.good()) {
    getline(fin,line);
    if (line == "") continue;
    ss.str(line);
    if (line.find("Stack:")!=string::npos) {
      ntag = tag.c_str() + to_string(istk);
      bd_tag=(*it_bead3d).bead_stack[istk].tag;
      cout<<"  ID= "<<setw(ndigit)<<std::right<<istk<<" tag= "<<ntag<<" ";//endl; 
      fnet_stack[istk].readData(ntag,bd_tag,&fin);
      stk=istk; //for z-bond assingmnet below
      istk++;
    }
    if (line.find("zbonds")!=string::npos)  {//zbonds
      while (!ss.fail() && getline(fin,line)) {
	if (line=="") break;
	ss.str(line);
	for (i=0;i<8;i++) bnd_temp[i]=-1; // clear just in case	 
	ss >> bnd_temp[0] >> bnd_temp[1] >> bnd_temp[2] >> bnd_temp[3] >> bnd_temp[4] >>
	  bnd_temp[5] >> bnd_temp[6] >> bnd_temp[7];  //4 z-bonds per line, bond_a recorded
	
	//if bond not found in in bond_a nor bond_b, call registerbondZ
	for (i=0;i<8;i=i+2) {
	  j=i+1;  iBead=bnd_temp[i]; jBead=bnd_temp[j];
	  if (iBead==-1 || jBead==-1) continue;
	  
	  flag=false; 
	  if (bond_a[stk].find(iBead)!=bond_a[stk].end()) { //skip if bond exists
	    it_range=bond_a[stk].equal_range(iBead);
	    for (it=it_range.first; it!=it_range.second;it++) {
	      if (it->second==jBead) {flag=true; break; }
	    }
	  } //	  if (bond_a[stk].find(iBead)!=bond_a[stk].end()) { //skip if bond exists
	  if (flag) continue; //bond exists
	  registerOneBondZ(iBead,jBead,-1,stk); //no dist cutoff, based on read
	}
	ss.clear();
	if (ss.fail()) error_dodri("READ: error in reading bonds."); 
      }   //  while (!ss.fail() && getline(fin,line)) {
      cout<<"  "<<setw(ndigit)<<std::right<<stk<<": "<<bond_a[stk].size()
	  <<" z-bonds assigned."<<endl;
    } //     if (line.find("zbonds")!=string::npos)  {//zbonds
  } //   while (fin.good())   
  return; 
}

/************************************************************************/
void fnet3d::registerOneBondZ(int ib, int kb, double rCutSq, int istk)
{ /* For bead ib in fnet_stack[istk], if the squared xy-distance with kb in
     fnet_stack[istk+1] is less than rCutSq, put (ib,kb) in bond_a[istk] and 
    (kb,ib) in bond_b[istk+1]. */

   double rdum;

   if (rCutSq>0)
     rdum=getDistPos(fnet_stack[istk].cm[ib],fnet_stack[istk+1].cm[kb],1);
   else
     rdum=-1./TOL; // no distance cutoff 
   
   if (rdum<rCutSq) {
     bond_a[istk].insert(pair<int,int>(ib,kb));
     bond_b[istk+1].insert(pair<int,int>(kb,ib));
   }
}

/************************************************************************/
void fnet3d::removeFloatingFil(int segLength,int *id_vec) 
{ /* Remove fils with length < segLength. */ 

  int i,ndigit;
  ndigit=nDigit(nStack); 
  
  for (i=id_vec[0];i<=id_vec[1];i++) {
    cout<<"  ID= "<<setw(ndigit)<<i<<" tag= "<<fnet_stack[i].tag ;
    fnet_stack[i].removeFloatingFil(segLength);
  }
}

/************************************************************************/
void fnet3d::removeOneBondZ(int ib, int kb, int istk)
{ /* remove bond between bead ib on istk and kb on istk+1 */
  pair<multimap<int, int>::iterator, multimap<int, int>::iterator> it_range;
  multimap<int, int>::iterator it;
  it_range = bond_a[istk].equal_range(ib);
  for (it=it_range.first; it!=it_range.second; it++) {
    if (it->second == kb) {
      bond_a[istk].erase(it);
      break;
    }
  }
  it_range = bond_b[istk+1].equal_range(kb);
  for (it = it_range.first; it!=it_range.second; it++) {
    if (it->second == ib) {
      bond_b[istk+1].erase(it);
      break;
    }
  }
  return;
}

/************************************************************************/
void fnet3d::removeUnZbondedBeads(int *id_vec)
{ /* Remove beads that have 0 z-bonds. */ 
  
  int istk, iBead,i,j,nBead0,iFil;
  vector<int> markbead; 
  multimap<int,int>::iterator it,zt;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> it_range,zt_range;
  vector<int>::iterator vt;
  
  int delBead,updateBead; 
  multimap<int,int> renumber;
  
  for (istk=id_vec[0];istk<=id_vec[1];istk++) {
    markbead.clear(); renumber.clear();
    for (iFil=1;iFil<=fnet_stack[istk].nFilament;iFil++) {
      it_range=fnet_stack[istk].fil_bead.equal_range(iFil);
      for (it=it_range.first;it!=it_range.second;it++) {
	iBead=it->second;
	if (bond_a[istk].find(iBead)==bond_a[istk].end() &&
	    bond_b[istk].find(iBead)==bond_b[istk].end())
	  markbead.push_back(iBead);
      }
    }

    for (i=0;i<(int)markbead.size();i++)  {
      nBead0=fnet_stack[istk].nBead-1;
      for (j=i+1;j<(int)markbead.size();j++) {
     	if (markbead[j]==nBead0) { // once bead in pos markbead[i] is deleted,
     	  // bead in nBead0 will be moved to pos i; update pos in markbead if found
     	  markbead[j]=markbead[i];
     	  break;
     	}
      }
      renumber.insert(pair<int,int>(nBead0,markbead[i]));                  
      fnet_stack[istk].deleteBead(markbead[i],1);
    }
    
    // take care of z-bonds
    for (zt=renumber.begin();zt!=renumber.end();zt++) {
      // delete all z-bonds of deleted (first) entry; update second entry
      delBead=zt->first; updateBead=zt->second;
      
      while (bond_a[istk].find(delBead)!=bond_a[istk].end()) { // bond_a
	it=bond_a[istk].find(delBead);
	removeOneBondZ(delBead,it->second,istk);
	if (renumber.find(updateBead)==renumber.end()) { 
	  // if updateBead is also a deleteBead, renumbering deleted this bead
	  // thus its z-bonds should not be registered
	  registerOneBondZ(updateBead,it->second,-1,istk); 
	}
      }
      while (bond_b[istk].find(delBead)!=bond_b[istk].end()) { // bond_b
	it=bond_b[istk].find(delBead);
	removeOneBondZ(it->second,delBead,istk-1);
	registerOneBondZ(it->second,updateBead,-1,istk-1); 
      }
    } //     for (zt=renumber.begin();zt!=renumber.end();zt++) {
  } //   for (istk=id_vec[0];istk<=id_vec[1];istk++) {
  
  return;       
}

/************************************************************************/
void fnet3d::setz(double dz, string origin)
{ /* Set z-plane coordinates for each slice.
     origin: zCoord=0 for 0th bead ("ini"); (nStack-1)-th bead ("fin");
     at (nStack-1)/2) ("mid")
     dz: adds to previous bead in the stack. So, for "fin" and dz>0,
     0-th bead has negative zCoord.  */
  
  int i;
  if (origin=="ini") {
    for (i=0;i<nStack;i++) fnet_stack[i].zCoord=dz*(double)i;
  }
  else if (origin=="fin") {
    for (i=(nStack-1);i>=0;i--) fnet_stack[i].zCoord=-dz*(double)i;
  }
  else if (origin=="mid") {
    for (i=0;i<nStack;i++) 
      fnet_stack[i].zCoord=dz*((double)i-0.5*(double)(nStack-1));
  }
  else 
    error_dodri("  Unrecognized origin keyword: "+origin);
  
  cout<<"  setz: dz= "<<dz<<" origin= "<<origin<<
    " zCoord range: "<<fnet_stack[0].zCoord<<" to "
      <<fnet_stack[nStack-1].zCoord<<endl;
}

/************************************************************************/
void fnet3d::smoothenFilAverage(int segLength,double rCut,int *id_vec)
{
  int i,ndigit;
  ndigit=nDigit(nStack);
  for (i=id_vec[0];i<=id_vec[1];i++) {
    cout<<"  ID= " <<setw(ndigit)<<i<<" tag= "<<fnet_stack[i].tag;
    fnet_stack[i].smoothenFilAverage(segLength,rCut);
  }
}

/************************************************************************/
void fnet3d::smoothenZ(int segLength,int delZ, double rCut,int *id_vec)
{ /* Between delZ slices, check if any 
     seglength # of beads in a filament are sequentiallly not z-bonded to
     previous slice. If so, adjust bead position of next slice based on 
     fit from last z-bonded beads betwen slices. 
     
     This DOES NOT add/delete beads/bonds.

     Compares zbonds between slices istk-->istk+1 
     Modifies bead positions in slice istk+1

     segLength: minimum length of non-zbonded sequence to be considered 
     delZ: # slices to check for sequential smoothening. This should be 
           based on the expected # of slices affected by a kink/protrusion           
     rCut: distance cutoff for bead re-position. If moving distance>rCut, 
           bead will not be repositioned. This will not affect other beads 
	   in sequence. Scaled by beadDiam.  
  */

  int i,k,istk,jstk,kstk; 
  int jb,jb0,jb1,ifil;
  double *pos=new double[2],pos1[2];    
  vector<int> znghb,iSeq,jSeq,zSeq; // flanking beads z-bonded from jstk--istk
  vector<int> checked; // bead ends already checked
  vector<pair<int,int>> jbz; // store zbonded beads [stk,bead index]
  vector<pair<int,int>>:: iterator it;
  multimap<int,int> jbz0,jbz1,ztemp;

  // for fitQuad0
  int iN,nMove,zk,zb0,zb1;
  double dist,dist1,mdPt[2],rdum[2],rdum1[2],bLeng,disp,fract; 
  double *vel=new double[2], *acc=new double[2]; 
  vector<array<double,2>> newPos;
  
  cout<<"  smoothen_z: segLength="<<segLength<<" delZ="<<
    delZ<<" rCut="<<rCut<<endl; 

  if (segLength<3)
    error_dodri("fnet3d::smoothenZ: minimum segLength required is 3.");
  
  for (istk=id_vec[0];istk<id_vec[1];istk++) {
    checked.clear();     
    jstk=istk+1;    znghb.clear(); 
    for (jb=0;jb<fnet_stack[jstk].nBead;jb++) {
      if (bond_b[jstk].find(jb)==bond_b[jstk].end()) { // no z-bond, find nearest zbond
	znghb.clear();
	ifil=fnet_stack[jstk].bead_fil.find(jb)->second; 
	findNextZ(jstk,jb,ifil,istk,znghb); // get flanking zbonded beads
	if ((int)znghb.size()!=2) continue; // avoid "hanging" segments
	if (find(checked.begin(),checked.end(),znghb[0])!=checked.end() &&
	    find(checked.begin(),checked.end(),znghb[1])!=checked.end()) continue;
	checked.push_back(znghb[0]);  checked.push_back(znghb[1]); // save checked beads
	jSeq.clear();	
	trackFilSeq(jstk,znghb[0],znghb[1],jSeq); //consecutive seq of non-zbonded beads
	if ((int)jSeq.size()<segLength)  continue; // do not adjust
	
	/* CALC new position */
	jb0=jSeq.front(); jb1=jSeq.back(); // save bead indeces of contour end to fill
	nMove=(int)jSeq.size()-2; // adjust position of beads between two ends
	for (k=0;k<2;k++) { // reference positions
	  rdum[k]=fnet_stack[jstk].cm[jb0][k];
	  rdum1[k]=fnet_stack[jstk].cm[jb1][k];
	  mdPt[k]=(rdum[k]+rdum1[k])/2.;	          
	}
	dist=getDistPos(rdum,rdum1,0); // total contour length
	fract=(1./(nMove+1)); // fraction for new bead placement along endBd vector
	fitQuad0(rdum,mdPt,rdum1,vel,acc,2);
	
	for (iN=1;iN<=nMove;iN++) {
	  bLeng=iN*fract*dist;
	  for (i=0;i<2;i++) {
	    disp=fract*acc[i]*bLeng*bLeng;
	    pos[i]=rdum[i]+vel[i]*bLeng+disp;
	  }

	  for (k=0;k<2;k++)  // use midpoint to smoothen transition 
	    pos1[k]=(pos[k]+fnet_stack[jstk].cm[jSeq[iN]][k])/2.;
	  if (dist>TOL) 
	    fnet_stack[jstk].moveBeadPos(jSeq[iN],pos1); // avoid beads w/ same pos
	}
	/*********/

	/* Adjust bond_a upstream bead positions. Do this by filament sequence restricted 
	   by zbonds to original jb0 and jb1.  	*/
	jbz0.clear(); jbz1.clear();
	findSeqZ(jb0,jstk,jbz0,id_vec);
	findSeqZ(jb1,jstk,jbz1,id_vec);	
	
	kstk=jstk+delZ;
	for (zk=jstk+1;zk<=kstk;zk++) {
	  // if bond_a sequence ends before delZ slices, stop adjusting	  
	  if (jbz0.find(zk)==jbz0.end() || jbz1.find(zk)==jbz1.end()) break; 
	  zb0=jbz0.find(zk)->second;   zb1=jbz1.find(zk)->second; // udpate ref pos
	  zSeq.clear(); 
	  trackFilSeq(zk,zb0,zb1,zSeq);
	  if ((int)zSeq.size()>segLength*3) continue; 
	  
	  /*********/
	  /* CALC new position (same procedure as above) */
	  nMove=(int)zSeq.size()-2; // adjust all positions
	  for (k=0;k<2;k++) { // reference positions
	    rdum[k]=fnet_stack[zk].cm[zb0][k];
	    rdum1[k]=fnet_stack[zk].cm[zb1][k];
	    mdPt[k]=(rdum[k]+rdum1[k])/2.;	  
	  }
	  dist=getDistPos(rdum,rdum1,0); // total contour length
	  fract=(1./(nMove+1)); 
	  fitQuad0(rdum,mdPt,rdum1,vel,acc,2);
	  
	  for (iN=1;iN<=nMove;iN++) {
	    bLeng=iN*fract*dist;
	    for (i=0;i<2;i++) {
	      disp=fract*acc[i]*bLeng*bLeng;
	      pos[i]=rdum[i]+vel[i]*bLeng+disp;
	    }
	    
	    // for the bead about to be moved:
	    // if it is not directly z-bonded to original moving slice (jstk), do not move
	    ztemp.clear(); 
	    findSeqZ(zSeq[iN],zk,ztemp,id_vec); // bead to move: zSeq[iN]
	    if (ztemp.find(jstk)==ztemp.end()) continue;
	    
	    // distance cutoff for moving. Used to prevent "new" segments from being moved
	    dist1=getDistPos(fnet_stack[zk].cm[zSeq[iN]],pos,0);
	    
	    for (k=0;k<2;k++) 
	      pos1[k]=(pos[k]+fnet_stack[zk].cm[zSeq[iN]][k])/2.;
	    
	    if (dist1<rCut*fnet_stack[zk].beadDiam && dist>TOL) 
	      fnet_stack[zk].moveBeadPos(zSeq[iN],pos1);
	  }
	} //	for (zk=0;zk<delZ;zk++) { 
      } //       if (bond_b[jstk].find(jb)==bond_b[jstk].end()) { 
    } //      for (jb=0;jb<fnet_stack[jstk].nBead;jb++) {
  }

  return; 
}

/************************************************************************/
void fnet3d::trackFilSeq(int istk, int b0, int b1, vector<int> &list)
{ /* Get list of beads in between b0 and b1 in the same filament/
     For cyclic, check both sequences, return shorter one. */ 

  int i,j=-2,ifil,jfil,iloc=-1,jloc=-1,temploc,cyclic; // vars initialized w/ dummy vals
  int filLength, *filSeq=new int[fnet_stack[istk].maxFilLength];
  vector<int> list1; 
  
  list.clear();

  ifil=fnet_stack[istk].bead_fil.find(b0)->second;
  jfil=fnet_stack[istk].bead_fil.find(b1)->second;
  if (ifil!=jfil) return;   // verify ifil==jfil

  filLength=fnet_stack[istk].findFilamentSeq(ifil,filSeq);
  cyclic=(filLength<0)?1:0;  // record if filament is cyclic

  filLength=abs(filLength); // abs value of filLength
  for (i=0;i<filLength;i++) {
    if (filSeq[i]==b0) iloc=i;
    if (filSeq[i]==b1) jloc=i;
  }

  // record list in direction from increasing fil position (iloc<jloc; iloc-->jloc)
  temploc=jloc;
  if (iloc>jloc) { jloc=iloc; iloc=temploc; }
  for (i=iloc;i<=jloc;i++) list.push_back(filSeq[i]);

  // cyclic: save list1 with fil beads in other search direction (jloc-->iloc)
  if (cyclic) {
    list1.clear();  i=iloc; 
    while (j!=jloc) {
      if (i<0) j=filLength+i;
      else j=i;
      list1.push_back(filSeq[j]);
      i--;       
    }
    if ((int)list1.size()<(int)list.size()) list=list1; // return list1;
  }
  
  delete[] filSeq; 
  return; 
}

/************************************************************************/
void fnet3d::writeCor(string fname, string outMode, int *id_vec)
{ /* Write coordinate file of fnet data */
  
  FILE *ff;
  string ofname,segname,sdum;
  stringstream ss;

  double x1, y1, z1, massMin, massMax, occupancy;
  double massMin0,massMax0;
  int i,j,imin, imax,iBead=0,nBead0=0,ndigit;
  multimap<int,int>::iterator it;

  // determine filename extension
  ofname = addFileNameExt(fname,".cor");
  cout<<"  Write fnet3d cor to "<<ofname<<" outMode= "<<outMode<<endl;
  ff=fopen(ofname.c_str(),"w");

  // find min/max mass: normalized mass entered as occupancy (col 55-60)
  massMin=1e308; massMax=-1;
  for (i=id_vec[0];i<=id_vec[1];i++) {
    MinMax(fnet_stack[i].nBead, fnet_stack[i].beadMass, 
	   &imin, &imax, &massMin0, &massMax0);
    massMin=(massMin0<massMin)?massMin0:massMin;
    massMax=(massMax0>massMax)?massMax0:massMax;
  }

  for (i=id_vec[0];i<=id_vec[1];i++) nBead0+=fnet_stack[i].nBead;
  // use extended format for more than 99999 atoms
  if (nBead0>99999) fprintf(ff,"%10d%s\n",nBead0,"  EXT");
  else fprintf(ff,"%d\n",nBead0);

  ndigit=nDigit(nStack);

  for (i=id_vec[0];i<=id_vec[1];i++) {
    ss.str(""); ss.clear();
    if (ndigit<5) ss<<i;
    else ss<<hex<<i;
    segname=ss.str();
    for (j=0;j<fnet_stack[i].nBead;j++) {
      if (outMode=="verbose") {
	x1=fnet_stack[i].cm[j][0]; y1=fnet_stack[i].cm[j][1]; 
	z1=fnet_stack[i].zCoord;
      }
      else { // outMode=="image"
	y1=((double)(fnet_stack[i].rows)-fnet_stack[i].cm[j][0]); 
	x1=fnet_stack[i].cm[j][1]; z1=fnet_stack[i].zCoord;
      } 
      occupancy=fnet_stack[i].beadMass[j]/massMax;
      
      /* From coorio.src:
	 fm2='(2I5,1X,A4,1X,A4,3F10.5,1X,A4,1X,A4,F10.5)'
	 I,IRES,RES(IRES),ATYPE(I),Xyz(1,I),xYz(2,I),xyZ(3,I),SID,RID,WMAIN(I)
      */
      if (nBead0>99999) { // use extended format
	fprintf (ff,"%10d%10d  %-8s  %-8s%20.10f%20.10f%20.10f  %-8s  %-8d%20.10f\n",
		 ++iBead,j,"FN3D","O",x1,y1,z1,segname.c_str(),j,occupancy);
      }
      else {
	//    fprintf (ff,"%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %4s %-4d%10.5f\n",
	fprintf (ff,"%5d%5d %-4s %-4s%9.4f %9.4f %9.4f  %4s %-4d%10.5f\n",
		 ++iBead,j,"FN3D","O",x1,y1,z1,segname.c_str(),j,occupancy);
      }
    } // for (j=0;j<fnet_stack[i].nBead;j++) {
  } // for (i=id_vec[0];i<=id_vec[1];i++) {
  fclose(ff);
}

/************************************************************************/
void fnet3d::writeData(string fname, string bd_fname, int *id_vec)
{ /* Writes fnet data to a dat file with name fname. 
     Also writes bead cm from fnet-modified functions. 
     Filament info not saved.   */
  
  FILE *dat, *bd_dat;
  string datname, bd_datname,bd3d_tag0="";
  int i,ndigit,bcnt;
  //  std::size_t pos; 
  multimap<int,int>::iterator it; 

  // fnet fopen
  datname = addFileNameExt(fname,".dat");
  dat = fopen(datname.c_str(), "w");
  cout<<"  Write fnet3d data to "<<datname<<endl;
  
  // bead fopen 
  bd_datname = addFileNameExt(bd_fname,".dat");
  bd_dat = fopen(bd_datname.c_str(), "w");
  cout<<"  Write bead3d data to "<<bd_datname<<endl;
  
  // if bead3d tag exists
  if (findBead3d(&vbead3d,bead3d_tag)!=vbead3d.end()) {
    fprintf(dat, "nStack= %d\nbead3d_tag= %s\n", id_vec[1]-id_vec[0]+1,
	    bead3d_tag.c_str());
    fprintf(bd_dat, "nStack= %d\nimg3d_tag= %s\n", id_vec[1]-id_vec[0]+1,
	    (findBead3d(&vbead3d,bead3d_tag)->img3d_tag).c_str());
  }
  else { // no bead tag, use bead file name for bead tag
    bd3d_tag0="bd";
    fprintf(dat, "nStack= %d\nbead3d_tag= %s\n", id_vec[1]-id_vec[0]+1,bd3d_tag0.c_str());
    fprintf(bd_dat, "nStack= %d\nimg3d_tag= %s\n", id_vec[1]-id_vec[0]+1,img3d_tag.c_str());
  }    
  
  ndigit=nDigit(nStack); 
  for (i=id_vec[0]; i<=id_vec[1];i++) {
    cout << "  ID= " <<setw(ndigit)<<i<<" tag= " << fnet_stack[i].tag<<endl; 
    fprintf(dat, "\nStack: %d\n\n", i-id_vec[0]);
    fprintf(bd_dat, "\nStack: %d\n\n", i-id_vec[0]);
    fnet_stack[i].writeData(fname, bd_fname, dat, bd_dat); //writes 2 files: fnet and bead .dat

    //add z-bonds; use bond_a
    //print bond[istk] bond[istk+1], in list format (like psf) 
    //    fprintf(dat, "\n");
    fprintf(dat, "\nzbonds\n");
    bcnt=0; 
    for (it=bond_a[i].begin();it!=bond_a[i].end();it++) {
      fprintf(dat,"%8d%8d",it->first,it->second);
      if (bcnt%4==3) fprintf(dat,"\n");
      ++bcnt; 
    }
    fprintf(dat,"\n");     
  }

  fclose(dat);
  fclose(bd_dat);
}

/************************************************************************/
void fnet3d::writeImage(int beadL, int bondL, int nColor, double *col_bg,
			double *col_bead, double *col_bond,
			string fname, string fmt, int *id_vec)
{ /* Write fnet3d into a set of images. */
  
  string ofname,sdum; stringstream ss;
  int i,ndigit,k;
  
  // count number of digits to use
  ss << id_vec[1]; sdum=ss.str(); ndigit=sdum.length(); 

  for (i=id_vec[0];i<=id_vec[1];i++) {
    ss.str(""); ss.clear();
    ss<<fname<<"_"<<setfill('0')<<setw(ndigit)<<i; //<<"."<<fmt;
    ofname=ss.str();
    cout << "  ID= " <<setw(ndigit)<<i <<"  tag= " << fnet_stack[i].tag;
    // when beadL=0, use beadDiam for each fnet_stack
    k=(beadL==0)?(fnet_stack[i].beadDiam):beadL; 
    fnet_stack[i].writeImage(k, bondL, nColor,col_bg,col_bead, col_bond,
			     ofname,fmt);
  }
  return;
}

/************************************************************************/
void fnet3d::writePsf(string fname, int *id_vec)
{ /* Writes psf file.
     id_vec[2]: begin/end stack numbers
     Use stack number as SEGID. If number of stacks is greater than 1e4,
     use hexadecimal expression to fit it into 4-letter SEGID. */
  
  int istk,j,ndigit,iBead=0,nBead0=0, **btemp; 
  int nb,nbz[nStack],nBond0=0,nStack0, bcnt;
  double rdum;
  string ofname,sdum;  stringstream ss;
  string segname,psfname;
  FILE *psf;
  multimap<int,int>::iterator it,zt;

  ofname = addFileNameExt(fname,".psf"); 
  cout<<"  Write fnet3d psf to "<<ofname<<endl;

  psf=fopen(ofname.c_str(),"w");
  fprintf(psf,"PSF \n\n !NTITLE\n");
  fprintf(psf,"* PSF for %s \n* \n\n",fname.c_str());

  // total number of beads & bonds in selection
  for (istk=id_vec[0];istk<=id_vec[1];istk++) {
    nBead0+=fnet_stack[istk].nBead;
    if (istk<id_vec[1]) {
      nBond0+=fnet_stack[istk].nBond/2+bond_a[istk].size();
    }
    else  // don't count z-bond for the last plane
      nBond0+=fnet_stack[istk].nBond/2;
  }
  
  fprintf(psf,"%6d !NATOM\n",nBead0);
  
  nStack0=id_vec[1]-id_vec[0]+1;
  ndigit=nDigit(nStack0);
  // cumulative number of beads before the current plane
  for (istk=id_vec[0];istk<=id_vec[1];istk++) {
    nbz[istk]=0;
    if (istk>id_vec[0]) {
      for (j=id_vec[0];j<istk;j++) nbz[istk]+=fnet_stack[j].nBead; 
    }
  }

  for (istk=id_vec[0];istk<=id_vec[1];istk++) {
    ss.str(""); ss.clear();
    if (ndigit<5) ss<<istk;
    else ss<<hex<<istk;
    segname=ss.str();
    for (j=0;j<fnet_stack[istk].nBead;j++) {
      rdum=fnet_stack[istk].beadMass[j];
      fprintf(psf,"%8d %4s %04d BD3D O       0 %14.6f %14.6f %d\n",
	      ++iBead,segname.c_str(),j,0.0,rdum,0);
      /* from psfres.src: fmt01=(I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8)
	 I(8),SEGID(ISEG)(4),RESID(IRES)(4),RES(IRES)(4),
	 ATYPE(I)(4),IAC(I)(4; atomID),CG(I)(14.6;charge),
	 AMASS(I)(14.6;mass),IMOVE(I)(8)
      */
    } // for (j=0;j<fnet_stack[istk].nBead;j++) {
  } // for (istk=id_vec[0];istk<=id_vec[1];istk++) {
  fprintf(psf,"\n%8d !NBOND\n",nBond0);
  
  btemp=new int*[nBond0]; nb=0; bcnt=0; // bcnt: bond number counter
  if (nBond0>0) { // case w/ bonds
    //    btemp=new int*[nBond0]; nb=0; bcnt=0; // bcnt: bond number counter
    for (j=0;j<nBond0;j++) btemp[j]=new int[2];
    for (istk=id_vec[0];istk<=id_vec[1];istk++) {
      // btemp: stores bead pairs for in-plane bonds
      for (j=0;j<fnet_stack[istk].nBond/2;j++) btemp[j][0]=btemp[j][1]=0;
      nb=0;
      for (it=fnet_stack[istk].bond.begin();
	   it!=fnet_stack[istk].bond.end(); it++) {
	for (j=0;j<nb;j++) { // check for duplicate bond
	  if ((btemp[j][0]==it->first)&&(btemp[j][1]==it->second)) 
	    goto WRITEPSF1;
	  if ((btemp[j][1]==it->first)&&(btemp[j][0]==it->second)) 
	    goto WRITEPSF1;
	}
	btemp[nb][0]=it->first;	btemp[nb++][1]=it->second;
      WRITEPSF1: continue;
      }
      assert (nb==(fnet_stack[istk].nBond/2));
      for (j=0;j<nb;j++) {
	fprintf(psf,"%8d%8d",btemp[j][0]+1+nbz[istk],btemp[j][1]+1+nbz[istk]);
	if(bcnt%4==3) fprintf(psf,"\n"); 
	++bcnt; 
      }
    } // for (istk=id_vec[0];istk<=id_vec[1];istk++) {
    
    for (istk=id_vec[0];istk<id_vec[1];istk++) { // write z-bond
      for (zt=bond_a[istk].begin();zt!=bond_a[istk].end();zt++) {
	fprintf(psf,"%8d%8d",(zt->first)+nbz[istk]+1,
		(zt->second)+nbz[istk+1]+1);
	if(bcnt%4==3) fprintf(psf,"\n"); 
	++bcnt; 
      }
    } // for (istk=id_vec[0];istk<id_vec[1];istk++) { // write z-bond
  } // if (nBond0>0) { // case w/ bonds
  
  
  fclose(psf);
  
  // free up memory
  if (nBond0>0) {
    for (j=0;j<nBond0;j++) delete [] btemp[j];
    delete [] btemp;
  }
  return;
}



/************************************************************************/
// Non-member functions
/************************************************************************/

void checkDuplicateFnet3dTag(vector<fnet3d> *vfnet3d, string tag)
{ /* In vfnet3d, find iterator for the element with tag for duplicate */
  vector<fnet3d>::iterator it_fnet3d;
  for (it_fnet3d=(*vfnet3d).begin();it_fnet3d!=(*vfnet3d).end();++it_fnet3d) {
    if ( (*it_fnet3d).tag==tag) error_dodri("Duplicate fnet3d tag: " + tag);
  }
  return;
}

/************************************************************************/
vector<fnet3d>::iterator findFnet3d(vector<fnet3d> *vfnet3d,string tag)
{ /* In vfnet3d, find iterator for the element with tag */
  vector<fnet3d>::iterator it_fnet3d;
  for (it_fnet3d=(*vfnet3d).begin();it_fnet3d!=(*vfnet3d).end();++it_fnet3d)
    if ( (*it_fnet3d).tag==tag) break;
  return it_fnet3d;
}
