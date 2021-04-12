#include "dodri.h"

int fnet3d_cmd(stringstream& ss, ifstream *fin0) {
  int flag,bflag,idum,ivec[2],ivec1[2], ivec2[2];
  double avec[10],avec1[3],avec2[3];
  string s_action, tag,tag1,sdum,sdum1,sdum2,ifname,ofname;
  vector<fnet3d>::iterator it_fnet3d;
  vector<bead3d>::iterator it_bead3d;
  vector<img3d>::iterator it_img3d; 
  
  tag=getCommOptString(ss,"tag","");
  if (tag=="") error_dodri("Missing tag.");

  bflag=0; //check for do block
  if (checkCommOptKey(ss,"do",0)==0) goto FNET3D0; //no do block
  checkCommand(ss);
  
  bflag=1;
  while (bflag==1) {
    sdum=""; 
    processCommand(ss, *fin0);
    if (ss.str().size() >0)  cout<<"  "<<tag<<": "<< ss.str() <<endl;
    else continue;
    
  FNET3D0: it_fnet3d=findFnet3d(&vfnet3d,tag); // find fnet3d with tag
    while(ss.good()) {
      s_action=getCommNext(ss); if (s_action.length()==0) continue;
      
      // check if tag exists
      if ((s_action!="stop")&&(s_action!="done")&&(s_action!="build")
	  &&(s_action!="read")) {
	if (it_fnet3d==vfnet3d.end())
	  error_dodri("Cannot find fnet3d "+ tag);
      }
      if (s_action=="stop") {
	cout<< "dodri:: Program finished."<<endl; exit(-1);
      }
      else if (s_action=="done") {
	bflag=0; checkCommand(ss); break;
      }

      /************ Constructor ************/
      
      else if (s_action=="build") {
	checkDuplicateFnet3dTag(&vfnet3d,tag);
	flag=checkCommOptKey(ss,"img3d",1); //flag or using image
	if (!flag) { // build fnet from bead
	  tag1=getCommOptString(ss,"bead3d",""); // tag of bead3d to use
	  it_bead3d=findBead3d(&vbead3d,tag1); // find bead3d with tag
	  if (it_bead3d==vbead3d.end()) 
	    error_dodri("findBead3d: cannot find bead3d "+ tag1);
	  getBead3dIdRange(ss,ivec,it_bead3d);
	  vfnet3d.push_back(fnet3d(&(*it_bead3d),ivec,tag));
	}
	else {  // build from image
	  if (checkCommOptKey(ss,"contour",0)) { 
	    tag1=getCommOptString(ss,"img3d",""); // tag of image to use
	    it_img3d=findImg3d(tag1);
	    if (it_img3d==vimg3d.end()) error_dodri("Cannot find image "+tag1);
	    getCommOptInt(ss,"size",&ivec[1],1,1); // iDiam
	    getCommOptInt(ss,"pxl_cut",ivec,1,1); 
	    vfnet3d.push_back(fnet3d(&(*it_img3d),ivec[1],
				     (unsigned char)ivec[0],tag)); 
	  }
	  else { error_dodri("Fnet build option not found."); }
	}
	checkCommand(ss);break;
      } // if (s_action=="build") {
      
      /************ Member Function ************/

      else if (s_action=="assign_bond_z") {
	getCommOptDouble(ss,"rcut",avec,1,1.5);
	getFnet3dIdRange(ss,ivec1,it_fnet3d);
	flag=checkCommOptKey(ss,"new_only",0); // restrict new zbonds 
	it_fnet3d->assignBondZ(avec[0],flag,ivec1);
	checkCommand(ss);break;
      } //       else if (s_action=="assign_bond_z") {
      
      else if (s_action=="clean_filaments") {
	getCommOptInt(ss,"length_cut",ivec,1,4); 
	getFnet3dIdRange(ss,ivec1,it_fnet3d);
	it_fnet3d->removeFloatingFil(ivec[0],ivec1);
	checkCommand(ss);break;
      } //       else if (s_action=="clean_filaments") {

      else if (s_action=="clear_bond_z") {
	getFnet3dIdRange(ss,ivec1,it_fnet3d);
	it_fnet3d->clearBondZ(ivec1);
	checkCommand(ss);break;
      } // else if (s_action=="clear_bondz_z") {
      
      else if (s_action=="equalize_bond") {
	getCommOptDouble(ss,"rcut",&avec[0],1,1.0);
	getFnet3dIdRange(ss,ivec1,it_fnet3d);
	it_fnet3d->equalizeFilamentBond(avec[0],ivec1);
	checkCommand(ss);break;
      } //       else if (s_action=="equalize_bond") {

      else if (s_action=="fill_enclosed_hole") {
	getFnet3dIdRange(ss,ivec,it_fnet3d);
	it_fnet3d->fillEnclosedHole(ivec); 
	checkCommand(ss);break;	
      }

      else if (s_action=="find_filament") {
	getFnet3dIdRange(ss,ivec1,it_fnet3d);
	flag=checkCommOptKey(ss,"cyclic",0);
	it_fnet3d->findFilament(!flag,ivec1);
	checkCommand(ss);break;
      } // else if (s_action=="find_filament")
      
      else if (s_action=="join_cyclic_fil") {
	getCommOptDouble(ss,"rcut",&avec[0],1,1.0);
	it_fnet3d->joinCyclicFil(avec[0],ivec1);
	checkCommand(ss);break;
      } //       else if (s_action=="join_cyclic_fil") {
      
      else if (s_action=="prune_bond_z") {
	getFnet3dIdRange(ss,ivec1,it_fnet3d);
	it_fnet3d->pruneBondZ(ivec1);
	checkCommand(ss);break;
      } //       else if (s_action=="prune_bond_z") {
      
      else if (s_action=="read") {
      	ifname=getCommOptString(ss,"name","");
	sdum=getCommOptString(ss,"bead3d",""); //tag of bead to use 
      	if (ifname == "") error_dodri("fnet3d_cmd: No filename given.");
      	if (it_fnet3d==vfnet3d.end())
	  vfnet3d.push_back(fnet3d(tag,ifname,sdum));
      	else it_fnet3d->readData(ifname,sdum);
      	checkCommand(ss); break;
      } //      else if (s_action=="read") {
      
      else if (s_action=="remove_unzbonded_beads") {
	getFnet3dIdRange(ss,ivec1,it_fnet3d);
	it_fnet3d->removeUnZbondedBeads(ivec1);
	checkCommand(ss);break;
      } //       else if (s_action=="remove_unzbonded_beads") {
      
      else if (s_action=="setz") { // set z coords of beads
	sdum1=getCommOptString(ss,"origin","ini");
	getCommOptDouble(ss, "dz", avec, 1, 1);
	it_fnet3d->setz(avec[0],sdum1);
	checkCommand(ss);break;
      } //       else if (s_action=="setz") {

      else if (s_action=="smoothen_average") {
	getCommOptInt(ss,"seglength",ivec,1,4); 
	getCommOptDouble(ss,"rcut",&avec[0],1,-1.0);
	getFnet3dIdRange(ss,ivec1,it_fnet3d);
	it_fnet3d->smoothenFilAverage(ivec[0],avec[0],ivec1);
	checkCommand(ss);break;
      } //       else if (s_action=="smoothen_average") {
      
      else if (s_action=="smoothen_z") {
	getCommOptInt(ss,"seglength",ivec,1,0);
	getCommOptInt(ss,"delz",&ivec[1],1,1); 
	getCommOptDouble(ss,"rcut",&avec[0],1,0.5); 
	getFnet3dIdRange(ss,ivec1,it_fnet3d);
	it_fnet3d->smoothenZ(ivec[0],ivec[1],avec[0],ivec1);
	checkCommand(ss);break;
      } //       else if (s_action=="smoothen_z") {
      
      else if (s_action=="write") {
	ofname=getCommOptString(ss,"name",tag); // output filename
	sdum1=getCommOptString(ss,"type",""); // output type
	sdum2="";
	getFnet3dIdRange(ss,ivec,it_fnet3d);
	if (sdum1=="") error_dodri("No output type given");
	else if (sdum1=="psf") it_fnet3d->writePsf(ofname,ivec);
	else if ((sdum1=="cor")||(sdum1=="cor+psf")||(sdum1=="psf+cor")) {
	  sdum2=getCommOptString(ss,"mode","image"); // output type
	  it_fnet3d->writeCor(ofname,sdum2,ivec); // write cor
	  if (sdum1!="cor") it_fnet3d->writePsf(ofname,ivec); // psf too
	}
	else if (sdum1=="dat") {
	  sdum=getCommOptString(ss,"bd_name",""); // output bead filename
	  if (sdum=="") error_dodri("No output name for bead tag given."); 
	  it_fnet3d->writeData(ofname,sdum,ivec);
	}
	else if (sdum1=="image") {
	  sdum=getCommOptString(ss,"format",""); // no format given here
	  getCommOptInt(ss,"bead_size",ivec1,1,0);
	  getCommOptInt(ss,"bond_size",&ivec2[1],1,1);
	  if (checkCommOptKey(ss,"rgb",0)==0) { idum=1; } //grayscale
	  else { idum=3; } //color
	  getCommOptDouble(ss,"color_bg",avec,idum,1.0);
	  getCommOptDouble(ss,"color_bead",avec1,idum,0.);
	  getCommOptDouble(ss,"color_bond",avec2,idum,0.);
	  for (int i=0;i<idum;i++) {
	    if ((avec[i]>1.0)||(avec1[i]>1.0)||(avec2[i]>1.0))
	      error_dodri("img_cmd:RGB values must be <= 1.0");
	    if ((avec[i]<0.0)||(avec1[i]<0.0)||(avec2[i]<0.0))
	      error_dodri("img_cmd: RGB values must be >= 1.0");
	  }
	  it_fnet3d->writeImage(ivec1[0],ivec2[1],idum,avec,avec1,avec2,
			      ofname,sdum,ivec);
	} //	else if (sdum1=="image") 
	else error_dodri("Unrecognized output type: "+sdum1);
	checkCommand(ss);break;
      } //       else if (s_action=="write") {
      

      /************************/
      
      else { // no proper action keyword
	cout <<"WARNING: Unrecognized action keyword: "<<s_action<<endl;
	checkCommand(ss);break;
      }
    } //while(ss.good())  
    ss.clear();
    cout<<endl;
  } //while (bflag==1) 
  return 0;
}
/************************************************************************/
void getFnet3dIdRange(stringstream& ss, int *id_vec,
		      vector<fnet3d>::iterator it_fnet3d)
{/* Find index range and save to id_vec :
    No range given (default): 0-(ns-1) (full range in it_fnet3d)
    fnet3d_id: work with id range within it_fnet3d 
    fnet: work with fnet tag range
    fnet_id: work with fnet id range */
  
  int ns=(*it_fnet3d).nStack,i,j,k,i0,j0, flag0,flag1,flag2;
  string sdum,sdum0,sdum1;
  stringstream ss0;
  size_t found; char *cdum;
  
  id_vec[0]=0; id_vec[1]=ns-1; // default
  flag0=checkCommOptKey(ss,"fnet3d_id",1); // priority ordered
  flag1=checkCommOptKey(ss,"fnet",1);
  flag2=checkCommOptKey(ss,"fnet_id",1);

  if (flag0==1) { // fnet3d_id
    if (flag1==1) {
      cout<<"  WARNING: Both fnet3d_id and fnet given.";
      cout<<"  Ignoring fnet."<<endl;
    }
    if (flag2==1) {
      cout<<"  WARNING: Both fnet3d_id and fnet_id given.";
      cout<<"  Ignoring fnet_id."<<endl;
    }
    sdum=getCommOptString(ss,"fnet3d_id","");
    found=sdum.find(":");
    if (found!=string::npos) { 
      sdum0=sdum.substr(0,found);
      sdum1=sdum.substr(found+1);
      i=(sdum0=="")?0:(int)strtol(sdum0.c_str(),&cdum,10);
      if (*cdum) error_dodri(sdum0+" is not a number.");
      j=(sdum1=="")?ns-1:(int)strtol(sdum1.c_str(),&cdum,10);
      if (*cdum) error_dodri(sdum0+" is not a number.");
      if (i>j) {
	cout<<"  WARNING: id_ini larger than id_fin. Swapping."<<endl;
	k=j; j=i; i=k;
      }
      if ((i<0)||(i>(ns-1))) {
	ss0<<"id_ini out of range. Must be between 0 and "<<(ns-1);
	error_dodri(ss0.str());
      }
      if ((j<0)||(j>(ns-1))) {
	ss0<<"id_fin out of range. Must be between 0 and "<<(ns-1);
	error_dodri(ss0.str());
      }
      id_vec[0]=i; id_vec[1]=j;
    }
    else error_dodri("fnet3d_id keyword must be followed by id_ini:id_fin.");
  }
   else { // fnet3d_id absent
    if (flag1==1) { // fnet
      if (flag2==1)
	cout<<"  Both fnet and fnet_id given. Ignoring fnet_id."<<endl;
      sdum=getCommOptString(ss,"fnet","");
      found=sdum.find(":");
      if (found!=string::npos) { 
	sdum0=sdum.substr(0,found);
	sdum1=sdum.substr(found+1);
	// only process provided tags
	if (sdum0!="") {
	  i=it_fnet3d->fnetTag2Id(sdum0);
	  if (i==-1) error_dodri("fnet tag "+sdum1+" not found in fnet3d "
				+(*it_fnet3d).tag);
	  else id_vec[0]=i;
	}
	if (sdum1!="") {
	  i=it_fnet3d->fnetTag2Id(sdum1);
	  if (i==-1) error_dodri("fnet tag "+sdum1+" not found in fnet3d "
				+(*it_fnet3d).tag);
	  else id_vec[1]=i;
	}
      }
      else error_dodri("fnet keyword must be followed by tag_ini:tag_fin.");
    } // if (flag1==1) { // fnet
    else { // fnet not given
      if (flag2==1) { // fnet_id
	sdum=getCommOptString(ss,"fnet_id","");
	found=sdum.find(":");
	if (found!=string::npos) { 
	  sdum0=sdum.substr(0,found);
	  sdum1=sdum.substr(found+1);
	  i=(sdum0=="")?0:(int)strtol(sdum0.c_str(),&cdum,10);
	  if (*cdum) error_dodri(sdum0+" is not a number.");
	  j=(sdum1=="")?(int)vfnet.size()-1:(int)strtol(sdum1.c_str(),&cdum,10);
	  if (*cdum) error_dodri(sdum0+" is not a number.");
	  if (i>j) {
	    cout<<"  WARNING: id_ini larger than id_fin. Swapping."<<endl;
	    k=j; j=i; i=k;
	  }
	  k=(int)vfnet.size()-1;
	  if ((i<0)||(i>k)) {
	    ss0<<"id_ini out of range. Must be between 0 and "<<k;
	    error_dodri(ss0.str());
	  }
	  if ((j<0)||(j>k)) {
	    ss0<<"id_fin out of range. Must be between 0 and "<<k;
	    error_dodri(ss0.str());
	  }
	  // find corresponding fnet tags
	  sdum0=vfnet[i].tag;
	  sdum1=vfnet[j].tag;
	  i0=it_fnet3d->fnetTag2Id(sdum0);
	  j0=it_fnet3d->fnetTag2Id(sdum1);
	  if (i0==-1) {
	    ss0<<"fnet_id "<<i<< " not found in fnet3d "+(*it_fnet3d).tag;
	    error_dodri(ss0.str());
	  }
	  else id_vec[0]=i0;
	  if (j0==-1) {
	    ss0<<"fnet_id "<<j<< " not found in fnet3d "+(*it_fnet3d).tag;
	    error_dodri(ss0.str());
	  }
	  else id_vec[1]=j0;
	} // if (found!=string::npos) {
	else error_dodri("fnet_id keyword must be followed by id_ini:id_fin.");
      } // if (flag2==1) { // fnet_id
    } // else { // fnet not given
   } // else { // fnet3d_id absent
  cout<<"  fnet3d_id selected: "<<id_vec[0]<<" - "<<id_vec[1]<<endl;
  return;
}


