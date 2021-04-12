#include "dodri.h"

int fnet_cmd(stringstream& ss, ifstream *fin0) {
  int bflag,flag,ivec[2],idum; //,ivec1[2],ivec2[2];
  double avec[4],avec1[10],avec2[3];
  string s_action, tag,tag1,ifname,ofname,sdum,sdum1,sdum2; 
  vector<bead>::iterator it_bead;
  vector<fnet>::iterator it_fnet;
  vector<img>::iterator it_img;   
  
  tag=getCommOptString(ss,"tag","");
  if (tag=="") error_dodri("Missing tag.");

  bflag=0; //check for do block
  if (checkCommOptKey(ss,"do",0)==0) goto FNET0;
  checkCommand(ss);
  
  bflag=1; //do block
  while (bflag==1) {
    sdum=""; 
    processCommand(ss, *fin0);
    if (ss.str().size() >0)  cout<<"  "<<tag<<": "<< ss.str() <<endl;
    else continue;
	  
  FNET0: it_fnet=findFnet(&vfnet,tag);
    while (ss.good()) {
      s_action=getCommNext(ss); if (s_action.length()==0) continue;

      // check if tag exists
      if ((s_action!="stop")&&(s_action!="done")&&(s_action!="build")
	  &&(s_action!="read")) 
	if (it_fnet==vfnet.end()) error_dodri("cannot find fnet "+ tag);
      if (s_action=="stop") {
	cout<< "Program finished."<<endl; exit(-1);
      }
      else if (s_action=="done") {
	bflag=0; checkCommand(ss);break;
      }      
      
      /************ Constructor ************/
      
      else if (s_action=="build") { 
	checkDuplicateFnetTag(&vfnet,tag);
	flag=checkCommOptKey(ss,"img",1); // flag for using image or not
	if (!flag) { // build fnet from bead
	  tag1=getCommOptString(ss,"bead",""); // tag of bead to use
	  if (tag1=="") error_dodri("Missing bead tag ");	  	  
	  it_bead=findBead(&vbead,tag1); // find bead with tag1
	  if (it_bead==vbead.end()) 
	    error_dodri("findBead: cannot find bead "+ tag1);
	  vfnet.push_back(fnet(&(*it_bead),tag));
	}
	else { // build from image
	  tag1=getCommOptString(ss,"img",""); // tag of image to use
	  if (tag1=="") error_dodri("Missing img tag. ");
	  it_img=findImg(tag1); // find img with tag1
	  if (it_img==vimg.end()) 
	    error_dodri("findImg: cannot find img "+ tag1);	    
	  
	  if (checkCommOptKey(ss,"contour",0)) { // accg2010
	    getCommOptInt(ss,"size",&ivec[1],1,1); // iDiam
	    getCommOptInt(ss,"pxl_cut",ivec,1,1);	    
	    vfnet.push_back(fnet(&(*it_img),ivec[1],
				 (unsigned char)ivec[0],tag));
	  }
	  else error_dodri("Invalid fnet build option."); 
	}
	checkCommand(ss);break;
      } // if (s_action=="build") { 
      
      /************ Member functions ************/
      
      else if (s_action=="clean_filaments") {
	getCommOptInt(ss,"length_cut",ivec,1,4); 
	it_fnet->removeFloatingFil(ivec[0]);
	checkCommand(ss);break;
      } //      else if (s_action=="clean_beads") {
      
      else if (s_action=="equalize_bond") {
	getCommOptDouble(ss,"rcut",&avec[0],1,1.0);	
	it_fnet->equalizeFilamentBond(avec[0]);
	checkCommand(ss);break;
      } //      else if (s_action=="equalize_bond") {
      
      else if (s_action=="find_filament") {
	flag=checkCommOptKey(ss,"cyclic",0); 
	it_fnet->findFilament(!flag);
	checkCommand(ss);break;
      } //      else if (s_action=="find_filament") {
      
      else if (s_action=="join_cyclic_fil") {
	getCommOptDouble(ss,"rcut",&avec[0],1,1.);
	it_fnet->joinCyclicFil(avec[0]); 
	checkCommand(ss);break;
      } //      else if (s_action=="close_fil_gaps") {

      else if (s_action=="read") {
	ifname=getCommOptString(ss,"name","");
	sdum=getCommOptString(ss,"bead",""); // tag of bead to use  
	if (ifname == "") error_dodri("fnet_cmd: No filename given.");
	if (it_fnet==vfnet.end()) vfnet.push_back(fnet(tag,ifname,sdum)); 
	else it_fnet->readData(ifname,sdum); 
	checkCommand(ss); break;
      } // else if (s_action=="read") {

      else if (s_action=="smoothen_average") {
	getCommOptInt(ss,"seglength",ivec,1,4); 
	getCommOptDouble(ss,"rcut",&avec[0],1,-1.0); 
	it_fnet->smoothenFilAverage(ivec[0],avec[0]);
	checkCommand(ss);break;
      } //       else if (s_action=="smoothen_average") {

      else if (s_action=="write") {
	ofname=getCommOptString(ss,"name",it_fnet->tag); // output filename
	sdum1=getCommOptString(ss,"type",""); // output type
	if (sdum1=="") error_dodri("  No output type given");
	else if (sdum1=="psf") it_fnet->writePsf(ofname);
	else if ((sdum1=="cor")||(sdum1=="pdb")) { // write coordinates
	  sdum=getCommOptString(ss,"occ","mass");
	  sdum2=getCommOptString(ss,"mode","image");
	  it_fnet->writeCor(ofname,sdum1,sdum,sdum2);
	}
	else if ((sdum1=="psf+cor")||(sdum1=="cor+psf")) {
	  sdum=getCommOptString(ss,"occ","mass");
	  sdum2=getCommOptString(ss,"mode","image");
	  it_fnet->writePsf(ofname);
	  it_fnet->writeCor(ofname,"cor",sdum,sdum2);
	}
	else if ((sdum1=="psf+pdb")||(sdum1=="pdb+psf")) {
	  sdum=getCommOptString(ss,"occ","mass");
	  sdum2=getCommOptString(ss,"mode","image");
	  it_fnet->writePsf(ofname);
	  it_fnet->writeCor(ofname,"pdb",sdum,sdum2);
	}
	else if (sdum1=="image") {
	  sdum=getCommOptString(ss,"format",""); // no format given here
	  getCommOptInt(ss,"bead_size",ivec,1,it_fnet->beadDiam);
	  idum=((it_fnet->beadDiam)>1)?((it_fnet->beadDiam)/2):1;
	  getCommOptInt(ss,"bond_size",&ivec[1],1,idum);
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
	  it_fnet->writeImage(ivec[0],ivec[1],idum,avec,avec1,avec2,
			     ofname,sdum);
	} // 	else if (sdum1=="image") {
	else if (sdum1=="dat") {
	  sdum=getCommOptString(ss,"bd_name",""); // output filename
	  if (sdum=="") error_dodri("No output name for bead tag given."); 
	  it_fnet->writeData(ofname,sdum);
	}
	else error_dodri("Unrecognized output type: "+sdum1);
	checkCommand(ss);break;
      } // if (s_action=="write") {
      
      else { // no proper action keyword
	cout <<"WARNING: Unrecognized action keyword: "<<s_action<<endl;
	checkCommand(ss); break;
      }
    } // while (ss.good()) {
    ss.clear();
    cout<<endl; // add a blank line between command execution    
  } // whille (bflag==1) {
  return 0;
}
