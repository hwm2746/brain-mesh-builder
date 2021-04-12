#include "dodri.h"

int img_cmd(stringstream& ss, ifstream *fin0) { 
  int flag,bflag,idum,ivec[2],ivec1[4];
  double avec[3],avec1[3]; 
  string s_action,tag,tag1,ifname,ofname,sdum,sdum1,sdum2;
  vector<img>::iterator it_img;
  ifstream ifile;
  
  tag=getCommOptString(ss,"tag","");
  if (tag=="") error_dodri("Missing tag.");
  
  bflag=0; // check for do block
  if (checkCommOptKey(ss,"do",0)==0) goto IMG0; // no "do" block
  checkCommand(ss);
  
  bflag=1;  // "do" block
  while (bflag==1) {
    sdum=""; 
    processCommand(ss, *fin0);
    if (ss.str().size() >0)  // output tag info
      cout<<"  "<<tag<<": "<< ss.str() <<endl;
    else continue;
    
  IMG0: it_img=findImg(tag);
    while(ss.good()) {
      s_action=getCommNext(ss); if (s_action.length()==0) continue;
      
      // check if tag exists
      if ((s_action!="stop")&&(s_action!="done")&&(s_action!="read")) 
	if (it_img==vimg.end()) error_dodri("cannot find image "+ tag);
      if (s_action=="stop") {
	cout<< "Program finished."<<endl; exit(-1);
      }
      else if (s_action=="done") {
	bflag=0; checkCommand(ss); break;
      }
      
      /************ Constructor ************/
      
      else if (s_action=="read") {
	getCommOptInt(ss,"ncell",ivec,1,5); nCell=ivec[0];
	if (nCell<=0) error_dodri("nCell must be >=0"); 
	checkDuplicateImgTag(&vimg,tag); 	
	sdum1=getCommOptString(ss,"type","single"); // input type 
	ifname=getCommOptString(ss,"name",""); // input filename
	if (ifname=="") error_dodri("img_cmd: No input file given.");
	ifile.open(ifname);
	if (ifile.fail()) error_dodri("img_cmd: Image or stack file does not exist.");
	if (sdum1=="single") { 
	  assert(ifname!="");
	  vimg.push_back(img(ifname,tag)); 
	}
	else if (sdum1=="stack")  // image stack
	  readImageStack(&vimg,ifname,sdum1,tag); 
 	else 
	  error_dodri("img read: " + sdum1 + " : Unrecognized read option.");
	ifile.close();
	checkCommand(ss); break;
      } // if (s_action=="read") {
      
      /************ Member functions ************/
      
      else if (s_action=="fill_region") { // position-based region fill
	getCommOptInt(ss,"pixel_ini",ivec,2,-1);   // reference pixel
	if (ivec[0]<0 || ivec[1]<0) error_dodri("Reference pixel missing.");
	getCommOptInt(ss,"wsize",&ivec1[0],1,1);   // window size
	if (ivec1[0]<=0) error_dodri("wsize must be greater than zero.");
	if ((ivec1[0]> (int)(it_img->rows -1)) || (ivec1[0]> (int)(it_img->columns -1)))
	  error_dodri("wsize must be smaller than image size.");
	getCommOptInt(ss,"foreground",&ivec1[1],1,255);   // foreground fill value
	if (ivec1[1]<0) {ivec1[1]=-1;} if (ivec1[1]>255) {ivec1[1]=255;}
	if (ivec1[1]<0) {ivec1[1]=-1;} if (ivec1[1]>255) {ivec1[1]=255;}
	sdum2=getCommOptString(ss,"mode","region"); 
	getCommOptInt(ss,"background",&ivec1[2],1,0);   // background fill value
	if (ivec1[2]<0) {ivec1[2]=-1;} if (ivec1[2]>255) {ivec1[2]=255;}
	if (ivec1[2]<0) {ivec1[2]=-1;} if (ivec1[2]>255) {ivec1[2]=255;}
	
	tag1=getCommOptString(ss,"save",""); // Saving option	
	it_img->fillRegion(ivec,ivec1[0],ivec1[1],ivec1[2],sdum2,tag1);	
	checkCommand(ss); break;
      } // else if (s_action=="fill_region") { // position-based region fill

      else if (s_action=="filter") {
	sdum1=getCommOptString(ss,"kind","highpass"); // filter kind
	flag=checkCommOptKey(ss,"binary",0); // flag for binary or not 
	if (sdum1=="highpass"||sdum1=="bandpass"||sdum1=="lowpass") {
	  sdum2=getCommOptString(ss,"mode","relative"); // filter mode
	}
	tag1=getCommOptString(ss,"save",""); // Saving option
	if (sdum1=="highpass"||sdum1=="lowpass") { // single cutoff
	  getCommOptDouble(ss,"pxl_cut",avec,1,0);
	  it_img->threshold(sdum1,sdum2,avec,flag,tag1);
	} // if (sdum1=="highpass"||sdum1=="lowpass") { // single cutoff
	else if (sdum1=="bandpass") { // bandpass: two cutoffs
	  getCommOptDouble(ss,"pxl_cut",avec,2,0);
	  it_img->threshold(sdum1,sdum2,avec,flag,tag1);
	}
	else if (sdum1=="gaussian") { // blur filter, gaussBlur function
	  // Only one stdev argument is given
	  if (ss.str().find("sx")!=string::npos
	      && ss.str().find("sy")==string::npos ) {
	    getCommOptDouble(ss,"sx",avec,1,1.5);  // std of gaussian
	    avec[1]=avec[0];  }
	  else if (ss.str().find("sx")==string::npos
		   && ss.str().find("sy")!=string::npos ) {
	    getCommOptDouble(ss,"sy",&avec[1],1,1.5);
	    avec[0]=avec[1];	  }
	  else {
	    getCommOptDouble(ss,"sx",avec,1,1.5);  // std of gaussian
	    getCommOptDouble(ss,"sy",&avec[1],1,1.5);	  }
	  
	  if (avec[0]<0. || avec[1]<0.)
	    error_dodri("img_cmd: sx, sy must be positive.");
	  
	  // Only one width argument given
	  if (ss.str().find("wx")!=string::npos
	      && ss.str().find("wy")==string::npos ) {
	    getCommOptInt(ss,"wx",ivec,1,3);   // window size
	    ivec[1]=ivec[0];	  }
	  else if (ss.str().find("wx")==string::npos
		   && ss.str().find("wy")!=string::npos ) {
	    getCommOptInt(ss,"wy",&ivec[1],1,3);
	    ivec[0]=ivec[1];	  }
	  else {
	    getCommOptInt(ss,"wx",ivec,1,3);   // window size
	    getCommOptInt(ss,"wy",&ivec[1],1,3);	  }
	  
	  //check window sizes
	  if (ivec[0]%2==0 || ivec[1]%2==0 || ivec[0]<3 || ivec[1] <3)
	    error_dodri("img_cmd: wx, wy must be an odd integer 3 or above.");
	  getCommOptInt(ss,"pxl_cut",ivec1,1,0); // 0: include all pxls.
	  it_img->gaussBlur(avec,ivec,ivec1[0],tag1);
	} // else if (sdum1=="gaussian");
	else if (sdum1=="auto") {
	  getCommOptDouble(ss,"cutoff",avec,1,0.75);
	  getCommOptInt(ss,"nbin",ivec,1,32);
	  if (avec[0]<0. || avec[0]>1.)
	    error_dodri("cutoff must be between 0.0 and 1.0.");
	  it_img->thresholdAuto(avec[0],ivec[0],tag1);
	}
	else if (sdum1=="histo") {
	  getCommOptDouble(ss,"locut",avec,1,-1.);
	  getCommOptDouble(ss,"hicut",&avec[1],1,-1.);
	  getCommOptInt(ss,"nbin",ivec,1,32);
	  it_img->thresholdHisto(avec,ivec[0],tag1);
	}
	else {
	  error_dodri("Unknown filter kind: "+sdum1);
	}
	checkCommand(ss); break;
      } // else if (s_action=="filter") {

      else if (s_action=="invert") {//invert image from max intensity
	tag1 = getCommOptString(ss,"save",""); //output tag
  	it_img->invert(tag1);
	checkCommand(ss); break;
      } //       else if (s_action=="invert") {//invert image from max intensity

      else if (s_action=="merge") {
	sdum=getCommOptString(ss,"with","");
	sdum1=getCommOptString(ss,"mode","add");
	if (sdum=="") error_dodri("Destination image name missing.");
	tag1 = getCommOptString(ss,"save",""); //output tag
	mergeImg( &(*it_img), &(*findImg(sdum)),sdum1,tag1);
	checkCommand(ss); break;
      } //      else if (s_action=="merge") {

      else if (s_action=="shift") {//shift image 
	getCommOptInt(ss,"dx",ivec,1,0);
	getCommOptInt(ss,"dy",&ivec[1],1,0);
	getCommOptInt(ss,"fill",ivec1,1,0);
	tag1 = getCommOptString(ss,"save",""); //output tag
	sdum = getCommOptString(ss,"tag_ref",""); //reference tag for filling
	if (ivec1[0]<0) ivec1[0]=0;
	else if (ivec1[0]>255) ivec1[0]=255; 
  	it_img->shift(ivec,(unsigned char)ivec1[0],tag1,sdum);
	checkCommand(ss); break;
      } //else if (s_action=="shift"
      
      else if (s_action=="write") {
	sdum1=getCommOptString(ss,"format",""); // output format
	ofname=getCommOptString(ss,"name",it_img->tag); // output filename
	if (checkCommOptKey(ss,"rgb",0)==0) idum=1;
	else idum=3; //color
	getCommOptDouble(ss,"color_bg",avec,idum,0.0);
	getCommOptDouble(ss,"scale",avec1,idum,1.0);	
	it_img->writeImage(ofname,sdum1,idum,avec,avec1);
	checkCommand(ss); break;
      } // else if (s_action=="write") {

      /******************************************/
      else { // no proper action keyword
	cout <<"  WARNING: Unrecognized action keyword: "<<s_action<<endl;
	checkCommand(ss); break;
      }
    } // while(ss.good()) {
    ss.clear();
    cout<<endl; // add a blank line between command execution
  } //while (bflag==1)
  return 0;
}
