#include "dodri.h"

int bead_cmd(stringstream& ss, ifstream *fin0) {
  int bflag,idum,ivec[2], ivec1[2], ivec2[2];
  double avec[3],avec1[3]; //rdum ,avec2[3];
  string s_action, tag,tag1,ifname,ofname,sdum,sdum1,sdum2; 
  vector<bead>::iterator it_bead;
  vector<img>::iterator it_img;
    
  tag=getCommOptString(ss,"tag","");
  if (tag=="") error_dodri("Missing tag.");
  
  bflag=0; // check for do block
  if (checkCommOptKey(ss,"do",0)==0) goto BEAD0; // no do block
  checkCommand(ss);
  
  bflag=1; // do block
  while (bflag==1) {
    sdum=""; 
    processCommand(ss, *fin0);
    if (ss.str().size() >0)  cout<<"  "<<tag<<": "<< ss.str() <<endl;
    else continue;

  BEAD0: it_bead=findBead(&vbead,tag);
    while(ss.good()) {
      s_action=getCommNext(ss); if (s_action.length()==0) continue;

      // check if tag exists
      if ((s_action!="stop")&&(s_action!="done")&&(s_action!="build")
	  && (s_action!="read")) 
	if (it_bead==vbead.end()) error_dodri("cannot find bead "+ tag);
      if (s_action=="stop") {
	cout<< "Program finished."<<endl; exit(-1);
      }
      else if (s_action=="done") {
	bflag=0; checkCommand(ss); break;
      }
      
      /************ Constructor ************/
      
      else if (s_action=="build") {
	checkDuplicateBeadTag(&vbead,tag);
	sdum1=getCommOptString(ss,"mode","area"); 
	tag1=getCommOptString(ss,"img",""); // tag of image to use
	ivec[0]=0; getCommOptInt(ss,"size",ivec,1,0); // iDiam
	if (ivec[0]==0) error_dodri("dodri: Zero bead diameter.");
	getCommOptInt(ss, "dir", &ivec1[1], 1, 0);
	getCommOptInt(ss,"verbosity",ivec2,1,0); // print level		
	it_img=findImg(tag1); // find image with tag
	if (it_img==vimg.end()) 
	  error_dodri("Cannot find image "+ tag1);

	// flag for keeping original pxl_gray
	sdum=getCommOptString(ss, "save", "");
	checkDuplicateImgTag(&vimg, sdum);
	if (sdum1=="area") {
	  if ((ivec1[1]>=0)&&(ivec1[1]<4)) { //one direction
	    vbead.push_back(bead( &(*it_img),tag,ivec[0],sdum,ivec1[1]));
	  }
	  else  //combine all directions, -1 for unneeded vars
	    vbead.push_back(bead(&(*it_img),tag,ivec[0], sdum,
				 -1.,-1.,-1,-1,-1,-1.,ivec2[0],sdum1));
	} //	  if (sdum1=="area") {
	checkCommand(ss); break;
      } // if (s_action=="build")
      
      /************ Member functions ************/

      else if (s_action=="read") {
	ifname=getCommOptString(ss,"name","");
	sdum=getCommOptString(ss,"img",""); //tag of image to use 
	if (ifname == "") error_dodri("bead_cmd: No filename given.");
	if (it_bead==vbead.end()) vbead.push_back(bead(tag, ifname,sdum));
	else it_bead->readData(ifname,sdum);
	checkCommand(ss); break;
      } // else if (s_action=="read") {
      
      else if (s_action=="write") {
	ofname=getCommOptString(ss,"name",tag); // output filename
	sdum1=getCommOptString(ss,"type",""); // output type
	if (sdum1=="") error_dodri("bead_cmd: No output type given");
	else if (sdum1=="psf") it_bead->writePsf(ofname);
	else if (sdum1=="cor") {
	  sdum2=getCommOptString(ss,"mode","image"); // output type
	  it_bead->writeCor(ofname,sdum2);
	}
	else if ((sdum1=="psf+cor")||(sdum1=="cor+psf")) {
	  sdum2=getCommOptString(ss,"mode","image");
	  it_bead->writePsf(ofname);
	  it_bead->writeCor(ofname,sdum2);
	}
	else if (sdum1=="dat") it_bead->writeData(ofname);
	else if (sdum1=="image") {
	  sdum=getCommOptString(ss,"format",""); // no format given here
	  getCommOptInt(ss,"bead_size",ivec,1,it_bead->beadDiam);
	  if (checkCommOptKey(ss,"rgb",0)==0) { idum=1; } //grayscale
	  else { idum=3; } //color
	  getCommOptDouble(ss,"color_bg",avec,idum,1.0);
	  getCommOptDouble(ss,"color_bead",avec1,idum,0.);
	  for (int i=0;i<idum;i++) {
	    if ((avec[i]>1.0)||(avec1[i]>1.0))
	      error_dodri("bead_cmd: RGB values must be <= 1.0");
	    if ((avec[i]<0.0)||(avec1[i]<0.0))
	      error_dodri("bead_cmd: RGB values must be >= 1.0");
	  }
	  it_bead->writeImage(ivec[0],idum,avec,avec1, 
			      ofname,sdum);
	}
	else error_dodri(" Unrecognized output type: "+sdum1);
	checkCommand(ss); break;
      } //        else if (s_action=="write") {
      
      /************ Non-member functions ************/
      
      else if (s_action=="copy") {
	tag1 = getCommOptString(ss,"save",""); //output tag
	if (tag1=="")
	  error_dodri("bead_cmd: Copy was not specified a new tag name.");
	copyBead( (*it_bead), it_bead->cm, it_bead->beadMass, tag1);
	checkCommand(ss); break;
      } //       else if (s_action=="copy") {

      /************************/      
      
      else { // no proper action keyword
	cout <<"WARNING: Unrecognized action keyword: "<<s_action<<endl;
	checkCommand(ss); break;
      }
    } //while(ss.good())  
    ss.clear();
    cout<<endl; // add a blank line between command execution    
  } //while (flag==1) 
  return 0;
}

