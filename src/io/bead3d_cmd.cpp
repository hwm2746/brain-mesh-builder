#include "dodri.h"

int bead3d_cmd(stringstream& ss, ifstream *fin0) {
  int bflag,idum,ivec[2],ivec1[2],ivec3[2],ivec4[2];
  double avec[10],avec1[3]; // rdum,avec[10],avec1[3],avec2[3];
  string s_action, tag,tag1,ifname,ofname,sdum,sdum1,sdum2;
  //  string svec[2];
  vector<bead3d>::iterator it_bead3d;
  vector<img3d>::iterator it_img3d;
  
  tag=getCommOptString(ss,"tag","");  
  if (tag=="") error_dodri("Missing tag.");
  
  bflag=0; // check for do block
  if (checkCommOptKey(ss,"do",0)==0) goto BEAD3D0; // no do block
  checkCommand(ss);
  
  bflag=1;
  while (bflag==1) {
    sdum=""; 
    processCommand(ss, *fin0);
    if (ss.str().size() >0)  cout<<"  "<<tag<<": "<< ss.str() <<endl;
    else continue;

  BEAD3D0: it_bead3d=findBead3d(&vbead3d,tag); // find bead with tag
    while(ss.good()) {
      s_action=getCommNext(ss); if (s_action.length()==0) continue; 

      // check if tag exists
      if ((s_action!="stop")&&(s_action!="done")&&(s_action!="build")
	  &&(s_action!="read")) 
	if (it_bead3d==vbead3d.end())
	  error_dodri("Cannot find bead3d "+ tag);
      if (s_action=="stop") {
	cout<< "dodri:: Program finished."<<endl; exit(-1);
      }
      else if (s_action=="done") {
	bflag=0; checkCommand(ss); break;
      }

      /************ Constructor ************/
      else if (s_action=="build") {
	checkDuplicateBead3dTag(&vbead3d,tag);
	
	// tag range of image to use
	tag1=getCommOptString(ss,"img3d",""); 
	it_img3d=findImg3d(tag1);
	if (it_img3d==vimg3d.end()) 
	  error_dodri("Cannot find image "+ tag1);
	
	getImg3dIdRange(ss,ivec3,it_img3d);
	
	sdum1=getCommOptString(ss,"mode","area"); 
	ivec[0]=0; getCommOptInt(ss,"size",ivec,1,0); // iDiam
	if (ivec[0]==0) error_dodri("Zero bead diameter.");
	getCommOptInt(ss, "dir", &ivec1[1], 1, 0);
	getCommOptInt(ss,"verbosity",ivec4,1,1); // print level
	
	sdum=getCommOptString(ss, "save", "");
	if (sdum != "") checkDuplicateImg3dTag(sdum);
	if (sdum1=="area") {
	  if ((ivec1[1]>=0)&&(ivec1[1]<4))  //one direction
	    vbead3d.push_back(bead3d( &(*it_img3d),tag,ivec[0],"",ivec1[1]));
	  else //combine all directions
	    vbead3d.push_back(bead3d(&(*it_img3d),tag,ivec[0],sdum,ivec4[0],ivec3));
	} //	  if (sdum1=="area") {
	checkCommand(ss); break;
      } //  if (s_action=="build")
      
      /************ Member functions ************/

      else if (s_action=="orient") {
 	getCommOptDouble(ss, "tol", avec, 1, 0.1); // tolerance in radians
	sdum1=getCommOptString(ss,"to","x"); // to which axis to orient
	getCommOptInt(ss,"slice",ivec1,1,-1);
	getBead3dIdRange(ss,ivec,it_bead3d);	
	it_bead3d->orient(avec[0],sdum1,ivec1[0],ivec);
	checkCommand(ss);break; 	
      }
      
      else if (s_action=="read") {
 	ifname=getCommOptString(ss,"name","");
	sdum=getCommOptString(ss,"img3d",""); //tag of image to use 
 	if (ifname == "") error_dodri("bead3d_cmd: No filename given.");
 	if (it_bead3d==vbead3d.end())
	  vbead3d.push_back(bead3d(tag, ifname,sdum));
 	else it_bead3d->readData(ifname,sdum);
	checkCommand(ss); break;
      } //       else if (s_action=="read") {

      else if (s_action=="rotate") {
	getBead3dIdRange(ss,ivec,it_bead3d);
	getCommOptDouble(ss, "angle", avec, 1, 0);
	it_bead3d->rotate(avec[0],ivec);
	checkCommand(ss); break;
      } //       else if (s_action=="rotate") {

      else if (s_action=="setz") { // set z coords of beads
	sdum1=getCommOptString(ss,"origin","ini");
	getCommOptDouble(ss, "dz", avec, 1, 1);
	it_bead3d->setz(avec[0],sdum1);
	checkCommand(ss); break;
      } //       else if (s_action=="setz") { // set z coords of beads
      
      else if (s_action=="write") {
	ofname=getCommOptString(ss,"name",tag); // output filename
	sdum1=getCommOptString(ss,"type",""); // output type
	getBead3dIdRange(ss,ivec,it_bead3d);
	if (sdum1=="") error_dodri("No output type given");
	else if (sdum1=="psf") it_bead3d->writePsf(ofname,ivec);
	else if ((sdum1=="cor")||(sdum1=="cor+psf")||(sdum1=="psf+cor")) {
	  sdum2=getCommOptString(ss,"mode","image"); // output type
	  it_bead3d->writeCor(ofname,sdum2,ivec); // write cor
	  if (sdum1!="cor") it_bead3d->writePsf(ofname,ivec); // psf too
	}
	else if (sdum1=="dat") it_bead3d->writeData(ofname, ivec);
	else if (sdum1=="image") {
	  sdum=getCommOptString(ss,"format",""); // no format given here
	  getCommOptInt(ss,"bead_size",ivec1,1,0);
	  if (checkCommOptKey(ss,"rgb",0)==0) { idum=1; } //grayscale
	  else { idum=3; } //color
	  getCommOptDouble(ss,"color_bg",avec,idum,1.0);
	  getCommOptDouble(ss,"color_bead",avec1,idum,0.);
	  for (int i=0;i<idum;i++) {
	    if ((avec[i]>1.0)||(avec1[i]>1.0))
	      error_dodri("bead3d_cmd: RGB values must be <= 1.0");
	    if ((avec[i]<0.0)||(avec1[i]<0.0))
	      error_dodri("bead3d_cmd: RGB values must be >= 1.0");
	  }
	  it_bead3d->writeImage(ivec1[0],idum,avec,avec1,ofname,sdum,ivec);
	} // else if (sdum1=="image") {
	else error_dodri("Unrecognized output type: "+sdum1);
	checkCommand(ss); break;
      } //       else if (s_action=="write") {
      
      /************************/      
      else { // no proper action keyword
	cout <<"WARNING: Unrecognized action keyword: "<<s_action<<endl;
	checkCommand(ss); break;
      }
    } //while(ss.good())  
    ss.clear();
    cout<<endl; //blank line after each command 
  } //while (bflag==1) 
  return 0;
}

/************************************************************************/
void getBead3dIdRange(stringstream& ss, int *id_vec,
		      vector<bead3d>::iterator it_bead3d)
{/* Find index range and save to id_vec :
    No range given (default): 0-(ns-1) (full range in it_bead3d)
    bead3d_id: work with id range within it_bead3d 
    bead: work with bead tag range
    bead_id: work with bead id range */
  
  int ns=(*it_bead3d).nStack,i,j,k,i0,j0, flag0,flag1,flag2;
  string sdum,sdum0,sdum1;
  stringstream ss0;
  size_t found; char *cdum;
  
  id_vec[0]=0; id_vec[1]=ns-1; // default
  flag0=checkCommOptKey(ss,"bead3d_id",1); // priority ordered
  flag1=checkCommOptKey(ss,"bead",1);
  flag2=checkCommOptKey(ss,"bead_id",1);

  if (flag0==1) { // bead3d_id
    if (flag1==1)
      cout<<"  WARNING: Both bead3d_id and bead given. Ignoring bead."<<endl;
    if (flag2==1)
      cout<<"  WARNING: Both bead3d_id and bead_id given. Ignoring bead_id."<<endl;
    sdum=getCommOptString(ss,"bead3d_id","");
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
    else error_dodri("bead3d_id keyword must be followed by id_ini:id_fin.");
  }
   else { // bead3d_id absent
    if (flag1==1) { // bead
      if (flag2==1)
	cout<<"  Both bead and bead_id given. Ignoring bead_id."<<endl;
      sdum=getCommOptString(ss,"bead","");
      found=sdum.find(":");
      if (found!=string::npos) { 
	sdum0=sdum.substr(0,found);
	sdum1=sdum.substr(found+1);
	// only process provided tags
	if (sdum0!="") {
	  i=it_bead3d->beadTag2Id(sdum0);
	  if (i==-1) error_dodri("bead tag "+sdum1+" not found in bead3d "
				+(*it_bead3d).tag);
	  else id_vec[0]=i;
	}
	if (sdum1!="") {
	  i=it_bead3d->beadTag2Id(sdum1);
	  if (i==-1) error_dodri("bead tag "+sdum1+" not found in bead3d "
				+(*it_bead3d).tag);
	  else id_vec[1]=i;
	}
      }
      else error_dodri("bead keyword must be followed by tag_ini:tag_fin.");
    } // if (flag1==1) { // bead
    else { // bead not given
      if (flag2==1) { // bead_id
	sdum=getCommOptString(ss,"bead_id","");
	found=sdum.find(":");
	if (found!=string::npos) { 
	  sdum0=sdum.substr(0,found);
	  sdum1=sdum.substr(found+1);
	  i=(sdum0=="")?0:(int)strtol(sdum0.c_str(),&cdum,10);
	  if (*cdum) error_dodri(sdum0+" is not a number.");
	  j=(sdum1=="")?(int)vbead.size()-1:(int)strtol(sdum1.c_str(),&cdum,10);
	  if (*cdum) error_dodri(sdum0+" is not a number.");
	  if (i>j) {
	    cout<<"  WARNING: id_ini larger than id_fin. Swapping."<<endl;
	    k=j; j=i; i=k;
	  }
	  k=(int)vbead.size()-1;
	  if ((i<0)||(i>k)) {
	    ss0<<"id_ini out of range. Must be between 0 and "<<k;
	    error_dodri(ss0.str());
	  }
	  if ((j<0)||(j>k)) {
	    ss0<<"id_fin out of range. Must be between 0 and "<<k;
	    error_dodri(ss0.str());
	  }
	  // find corresponding bead tags
	  sdum0=vbead[i].tag;
	  sdum1=vbead[j].tag;
	  i0=it_bead3d->beadTag2Id(sdum0);
	  j0=it_bead3d->beadTag2Id(sdum1);
	  if (i0==-1) {
	    ss0<<"bead_id "<<i<< " not found in bead3d "+(*it_bead3d).tag;
	    error_dodri(ss0.str());
	  }
	  else id_vec[0]=i0;
	  if (j0==-1) {
	    ss0<<"bead_id "<<j<< " not found in bead3d "+(*it_bead3d).tag;
	    error_dodri(ss0.str());
	  }
	  else id_vec[1]=j0;
	} // if (found!=string::npos) {
	else error_dodri("bead_id keyword must be followed by id_ini:id_fin.");
      } // if (flag2==1) { // bead_id
    } // else { // bead not given
   } // else { // bead3d_id absent
  cout<<"  bead3d_id selected: "<<id_vec[0]<<" - "<<id_vec[1]<<endl;
  return;
}

