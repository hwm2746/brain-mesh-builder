#include "dodri.h"

int img3d_cmd(stringstream& ss, ifstream *fin0) {
  int flag,bflag, id_vec[2],ivec[2], ivec1[4],ivec2[2];
  double avec[3],avec1[3]; //,rdum;
  unsigned char pvec[2];   
  string s_action,tag,tag1,ifname,ofname,sdum,sdum1,sdum2; 
  string svec[2];
  vector<img3d>::iterator it_img3d;
  ifstream ifile;
  
  tag=getCommOptString(ss,"tag","");
  if (tag=="") error_dodri("Missing tag.");
  
  bflag=0; // check for do block
  if (checkCommOptKey(ss,"do",0)==0) goto IMG3D0; // no do block
  checkCommand(ss);
  
  bflag=1;
  while (bflag==1) {
    sdum=""; 
    processCommand(ss, *fin0);
    if (ss.str().size() >0)  cout<<"  "<<tag<<": "<< ss.str() <<endl;
    else continue;
    
  IMG3D0: it_img3d=findImg3d(tag);
    while(ss.good()) {
      s_action=getCommNext(ss); if (s_action.length()==0) continue;
      
      // check if tag exists
      if ((s_action!="stop")&&(s_action!="done")&&(s_action!="build")&&
	  (s_action!="read")) {
	if (it_img3d==vimg3d.end()) 
	  error_dodri("cannot find img3d "+ tag);
      }
      if (s_action=="stop") {
	cout<< "Program finished."<<endl; exit(-1);
      }
      else if (s_action=="done") {
	bflag=0; checkCommand(ss); break;
      }
      
      /************ Constructor ************/
      
      else if (s_action=="build") {
	checkDuplicateImg3dTag(tag); // check for duplicate tag
	getImgTagRange(ss,svec);
	vimg3d.push_back(img3d(svec[0],svec[1],tag));
	checkCommand(ss); break;
      }
      else if (s_action=="copy") {
	tag1 = getCommOptString(ss,"save",""); //output tag
	if (tag1=="")
	  error_dodri("img3d_cmd: Copy was not specified a new tag name.");
	getImg3dIdRange(ss,id_vec,it_img3d); 
	vimg3d.push_back(img3d(&(*it_img3d),tag1,id_vec));
	checkCommand(ss); break;
      } //else if (s_action=="copy")
      
      /************ Member functions ************/
      
      else if (s_action=="fill_region") {
	//if only one argument is given
	if (ss.str().find("pixel_fin")==string::npos 
	    && ss.str().find("pixel_ini")!=string::npos) {
	  getCommOptInt(ss,"pixel_ini",ivec,2,-1); //reference pixel
	  ivec2[0]=ivec[0]; ivec2[1]=ivec[1];	}
	else if (ss.str().find("pixel_ini")==string::npos
		 && ss.str().find("pixel_fin")!=string::npos) {
	  getCommOptInt(ss,"pixel_fin",ivec2,2,-1);   // reference pixel
	  ivec[0]=ivec2[0]; ivec[1]=ivec2[1];	}
	else {  //both given 
	  getCommOptInt(ss,"pixel_ini",ivec,2,-1); //reference pixel
	  getCommOptInt(ss,"pixel_fin",ivec2,2,-1);   // reference pixel
	}
	if (ivec[0]<0 || ivec[1]<0 || ivec2[0]<0 || ivec2[1]<0)
	  error_dodri("Reference pixel missing.");	
	
	getImg3dIdRange(ss,id_vec,it_img3d); //use full range
	getCommOptInt(ss,"wsize",&ivec1[0],1,1);   // window size
	if (ivec1[0]<=0) error_dodri("wsize must be greater than zero.");
	//wSize check moved to function	
	getCommOptInt(ss,"foreground",&ivec1[1],1,255);   // fill value
	//if <0, will keep original pixel value
	if (ivec1[1]<0) {ivec1[1]=-1;}	if (ivec1[1]>255) {ivec1[1]=255;} 
	if (ivec1[1]<0) {ivec1[1]=-1;}	if (ivec1[1]>255) {ivec1[1]=255;}
	getCommOptInt(ss,"background",&ivec1[2],1,0);   // fill value
	if (ivec1[2]<0) {ivec1[2]=-1;}	if (ivec1[2]>255) {ivec1[2]=255;} 
	if (ivec1[2]<0) {ivec1[2]=-1;}	if (ivec1[2]>255) {ivec1[2]=255;}
	sdum2=getCommOptString(ss,"mode","region"); // fill mode
	it_img3d->fillRegion2D(ivec,ivec2,ivec1[0],
			       ivec1[1],ivec1[2],sdum2,id_vec);
	checkCommand(ss); break;
      } //else if (s_action=="fill_region")
      
      else if (s_action=="filter") {
	getImg3dIdRange(ss,id_vec,it_img3d); flag=0;
	sdum1=getCommOptString(ss,"kind","highpass"); // filter kind
	if (sdum1=="highpass"||sdum1=="bandpass"||sdum1=="lowpass") {
	  sdum2=getCommOptString(ss,"mode","relative"); // filter kind
	  flag=checkCommOptKey(ss,"binary",0);
	}
	if (sdum1=="highpass"||sdum1=="lowpass") { // single cutoff
	  getCommOptDouble(ss,"pxl_cut",avec,1,0);
	  it_img3d->threshold(sdum1,sdum2,avec,flag,id_vec);
	} // 	if (sdum1=="highpass"||sdum1=="lowpass") { // single cutoff
	else if (sdum1=="bandpass") { // bandpass: two cutoffs
	  getCommOptDouble(ss,"pxl_cut",avec,2,0);
	  it_img3d->threshold(sdum1,sdum2,avec,flag,id_vec);
	}
      	else if (sdum1=="gaussian") { // blur filter, gaussBlur2D function
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
	    getCommOptDouble(ss,"sy",&avec[1],1,1.5);
	  }
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
	  if (ivec[0]%2==0 || ivec[1]%2==0 ||  ivec[0]<3 || ivec[1]<3 )
	    error_dodri("img3d_cmd: wx, wy must be an odd integer 3 or above.");
	  getCommOptInt(ss,"pxl_cut",ivec1,1,0); // 0: include all pxls.
	  it_img3d->gaussBlur2D(avec,ivec,ivec1[0],id_vec);
	  checkCommand(ss); break;
	} //if (sdum1=="gaussian");
	else if (sdum1=="auto") { 
	  getCommOptDouble(ss,"pxl_cut",avec,1,0.75);
	  getCommOptInt(ss,"nbin",ivec,1,128);
	  if (avec[0]<0. || avec[0]>1.)
	    error_dodri("cutoff must be between 0.0 and 1.0.");
	  it_img3d->thresholdAuto2D(avec[0],ivec[0],id_vec);
	  checkCommand(ss); break;
	}
	else if (sdum1=="histo") {
	  getCommOptDouble(ss,"locut",avec,1,-1.);
	  getCommOptDouble(ss,"hicut",&avec[1],1,-1.);
	  getCommOptInt(ss,"nbin",ivec,1,32);
	  it_img3d->thresholdHisto2D(avec,ivec[0],id_vec);
	}	
	else 
	  error_dodri("Unknown filter kind: "+sdum1);
	checkCommand(ss); break;
      } // else if (s_action=="filter") {
      
      else if (s_action=="invert") {
	getImg3dIdRange(ss,id_vec,it_img3d);
	it_img3d->invert2D(id_vec);
	checkCommand(ss); break;
      } //      else if (s_action=="invert") {
      
      else if (s_action=="merge") {
	svec[0]=getCommOptString(ss,"from_tag","");
	svec[1]=getCommOptString(ss,"to_tag","");
	if (svec[0]=="" || svec[1]=="")
	  error_dodri("Destination image name missing.");
	sdum=getCommOptString(ss,"mode","add");
	getImg3dIdRange(ss,id_vec,it_img3d); 
	it_img3d->merge(svec,sdum,id_vec);
	checkCommand(ss); break;
      } //else if (s_action=="merge")

      else if (s_action=="remove_debris") { 
	getCommOptInt(ss,"size_cut",&ivec[0],1,1);
	if (ivec[0]<0) error_dodri("Provide a positive integer cutoff.");
	getCommOptUchar(ss,"pxl_bg",&pvec[0],1,0);
	getCommOptUchar(ss,"pxl_cut",&pvec[1],1,0);
	getCommOptInt(ss,"fill",ivec1,1,0);
	sdum2=getCommOptString(ss,"mode","rank"); // area or rank
	getImg3dIdRange(ss,id_vec,it_img3d);
	it_img3d->removeDebris2D(ivec[0],pvec,(unsigned char)ivec1[0],
				 sdum2,id_vec);
	checkCommand(ss); break;
      } //      else if (s_action=="remove_debris") {
      
      else if (s_action=="shift") {
	getCommOptInt(ss,"dx",ivec,1,0);
	getCommOptInt(ss,"dy",&ivec[1],1,0);
	getCommOptInt(ss,"fill",ivec1,1,0);
	if (ivec1[0]<0) ivec1[0]=0;
	else if (ivec1[0]>255) ivec1[0]=255; 
	getImg3dIdRange(ss,id_vec,it_img3d);
	it_img3d->shift(ivec,(unsigned char)ivec1[0],id_vec);
	checkCommand(ss); break;
      } //else if (s_action=="shift")
      
      else if (s_action=="write") {
	sdum1=getCommOptString(ss,"format",""); // output format
	ofname=getCommOptString(ss,"name",it_img3d->tag); // output filename
	getImg3dIdRange(ss,id_vec,it_img3d);
	if (checkCommOptKey(ss,"rgb",0)==1) {  // colo
	  getCommOptDouble(ss,"color_bg",avec,3,0.0);	  	  
	  getCommOptDouble(ss,"scale",avec1,3,1.0);
	  it_img3d->writeImage(ofname,sdum1,3,avec,avec1,id_vec);
	}
	else {
	  if (sdum1!="emap") { //grayscale
	    getCommOptDouble(ss,"color_bg",avec,1,0.0);	  	  	    
	    getCommOptDouble(ss,"scale",avec1,1,1.0);
	    it_img3d->writeImage(ofname,sdum1,1,avec,avec1,id_vec);
	  }
	  else { //emap
	    getCommOptDouble(ss,"dimension",avec1,3,1.);
	    sdum2=getCommOptString(ss,"mode","image");
	    //getImg3dIdRange(ss,id_vec,it_img3d); //use full range
	    it_img3d->writeEmap(avec1,ofname,sdum2,id_vec);
	  }
	}
	checkCommand(ss); break;
      } // else if (s_action=="write") {
      
      /************************/
      
      else { // no proper action keyword
	cout <<"WARNING: Unrecognized action keyword: "<<s_action<<endl;
	checkCommand(ss); break;
      }
      
    } // while(ss.good()) {
    ss.clear();
    cout<<endl; // blank line after each command
  } //while (bflag==1)
  return 0;
}

/************************************************************************/
void getImg3dIdRange(stringstream& ss, int *id_vec,
		     vector<img3d>::iterator it_img3d)
{/* Find index range and save to id_vec :
    No range given (default): 0-(ns-1) (full range in it_img3d)
    img3d_id: work with id range within it_img3d 
    img: work with img tag range
    img_id: work with img id range
 */
  int ns=(*it_img3d).nStack,i,j,k,i0,j0, flag0,flag1,flag2;
  string sdum,sdum0,sdum1;
  stringstream ss0;
  size_t found; char *cdum;
  
  id_vec[0]=0; id_vec[1]=ns-1; // default
  flag0=checkCommOptKey(ss,"img3d_id",1); // priority ordered
  flag1=checkCommOptKey(ss,"img",1);
  flag2=checkCommOptKey(ss,"img_id",1);

  if (flag0==1) { // img3d_id
    if (flag1==1)
      cout<<"  WARNING: Both img3d_id and img given. Ignoring img."<<endl;
    if (flag2==1)
      cout<<"  WARNING: Both img3d_id and img_id given. Ignoring img_id."<<endl;
    sdum=getCommOptString(ss,"img3d_id","");
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
    else error_dodri("img3d_id keyword must be followed by id_ini:id_fin.");
  }
  else { // img3d_id absent
    if (flag1==1) { // img
      if (flag2==1)
	cout<<"  Both img and img_id given. Ignoring img_id."<<endl;
      sdum=getCommOptString(ss,"img","");
      found=sdum.find(":");
      if (found!=string::npos) { 
	sdum0=sdum.substr(0,found);
	sdum1=sdum.substr(found+1);
	// only process provided tags
	if (sdum0!="") {
	  i=it_img3d->imgTag2Id(sdum0);
	  if (i==-1) error_dodri("img tag "+sdum1+" not found in img3d "
				+(*it_img3d).tag);
	  else id_vec[0]=i;
	}
	if (sdum1!="") {
	  i=it_img3d->imgTag2Id(sdum1);
	  if (i==-1) error_dodri("img tag "+sdum1+" not found in img3d "
				+(*it_img3d).tag);
	  else id_vec[1]=i;
	}
      }
      else error_dodri("img keyword must be followed by tag_ini:tag_fin.");
    } // if (flag1==1) { // img
    else { // img not given
      if (flag2==1) { // img_id
	sdum=getCommOptString(ss,"img_id","");
	found=sdum.find(":");
	if (found!=string::npos) { 
	  sdum0=sdum.substr(0,found);
	  sdum1=sdum.substr(found+1);
	  i=(sdum0=="")?0:(int)strtol(sdum0.c_str(),&cdum,10);
	  if (*cdum) error_dodri(sdum0+" is not a number.");
	  j=(sdum1=="")?(int)vimg.size()-1:(int)strtol(sdum1.c_str(),&cdum,10);
	  if (*cdum) error_dodri(sdum0+" is not a number.");
	  if (i>j) {
	    cout<<"  WARNING: id_ini larger than id_fin. Swapping."<<endl;
	    k=j; j=i; i=k;
	  }
	  k=(int)vimg.size()-1;
	  if ((i<0)||(i>k)) {
	    ss0<<"id_ini out of range. Must be between 0 and "<<k;
	    error_dodri(ss0.str());
	  }
	  if ((j<0)||(j>k)) {
	    ss0<<"id_fin out of range. Must be between 0 and "<<k;
	    error_dodri(ss0.str());
	  }
	  // find corresponding img tags
	  sdum0=vimg[i].tag;
	  sdum1=vimg[j].tag;
	  i0=it_img3d->imgTag2Id(sdum0);
	  j0=it_img3d->imgTag2Id(sdum1);
	  if (i0==-1) {
	    ss0<<"img_id "<<i<< " not found in img3d "+(*it_img3d).tag;
	    error_dodri(ss0.str());
	  }
	  else id_vec[0]=i0;
	  if (j0==-1) {
	    ss0<<"img_id "<<j<< " not found in img3d "+(*it_img3d).tag;
	    error_dodri(ss0.str());
	  }
	  else id_vec[1]=j0;
	} // if (found!=string::npos) {
	else error_dodri("img_id keyword must be followed by id_ini:id_fin.");
      } // if (flag2==1) { // img_id
    } // else { // img not given
  } // else { // img3d_id absent
  cout<<"  img3d_id selected: "<<id_vec[0]<<" - "<<id_vec[1]<<endl;
  return;
}

