#include "dodri.h"

/************************************************************************/
// Constructors
/************************************************************************/

img3d::img3d(string from_tag, string to_tag, string tag0)
{ /* Constructor using existing img: 
     In vimg, select from_tag to to_tag. */
  
  vector<img>::iterator it_img,it_img_ini, it_img_fin;
  
  tag=tag0; nStack=0; img_stack.clear();
  for (it_img=vimg.begin();it_img!=vimg.end();++it_img) {
    if (it_img->tag==from_tag) it_img_ini=it_img; 
    if (it_img->tag==to_tag) {it_img_fin=it_img; ++it_img_fin; break;}
  }
  nStack=(int)distance(it_img_ini,it_img_fin);
  
  for (it_img=it_img_ini;it_img!=it_img_fin;++it_img)
    img_stack.push_back(*it_img);
  
  cout <<"  img: "<<from_tag<< " - "<< to_tag
       <<" ( "<<nStack<<" images ) assigned to img3d: "<<tag<<endl;
}

/************************************************************************/
img3d::img3d(string ifname, string tag0)
{ /* Constructor from ifname containing list of image filenames. */

  int i= readImageStack(&vimg,ifname,"",tag0);
  string from_tag=tag0+"0", to_tag=tag0+itostr(i-1);
  
  vector<img>::iterator it_img,it_img_ini, it_img_fin;
  
  tag=tag0; nStack=0; img_stack.clear();
  for (it_img=vimg.begin();it_img!=vimg.end();++it_img) {
    if (it_img->tag==from_tag) it_img_ini=it_img; 
    if (it_img->tag==to_tag) {it_img_fin=it_img; ++it_img_fin; break;}
  }
  nStack=(int)distance(it_img_ini,it_img_fin);
  
  for (it_img=it_img_ini;it_img!=it_img_fin;++it_img)
    img_stack.push_back(*it_img);
  
  cout <<"  img: "<<from_tag<< " - "<< to_tag
       <<" ( "<<nStack<<" images ) assigned to img3d: "<<tag<<endl;
}

/************************************************************************/
img3d::img3d(img3d *img3d0, string tag0, int *id_vec)
{ /* Copy constructor */
  
  int i,j;
  string from_tag,to_tag,tag1;
  vector<img>::iterator it_img,it_img1, it_img_ini, it_img_fin;
  
  from_tag=img3d0->img_stack[id_vec[0]].tag;
  to_tag=img3d0->img_stack[id_vec[1]].tag;

  tag=tag0; nStack=0; img_stack.clear();
  for (it_img=vimg.begin();it_img!=vimg.end();++it_img) {
    if (it_img->tag==from_tag) it_img_ini=it_img; 
    if (it_img->tag==to_tag) {it_img_fin=it_img; ++it_img_fin; break;}
  }
  nStack=(int)distance(it_img_ini,it_img_fin);
  string tag_list[nStack]; 

  i=0;
  for (it_img=it_img_ini;it_img!=it_img_fin;++it_img) {
    tag_list[i]=it_img->tag;
    i++;
  }
  assert(i==nStack);

  for (i=0;i<nStack;i++) {
    j=i+id_vec[0]; // in case there is an ID offset
    tag1=tag0+to_string(j);
    it_img=findImg(tag_list[i]);
    if (it_img==vimg.end())
      error_dodri("Image "+tag_list[i]+" not found."); 
    copyImg((*it_img),it_img->pxl_gray,tag1);
    it_img1=findImg(tag1);
    img_stack.push_back(*it_img1);
  }
  assert(j==id_vec[1]);
  
  from_tag=tag0+"0";to_tag=tag0+itostr(i-1);
  cout<<"  img:"<<from_tag<<" - "<<to_tag
      <<" ( "<<nStack<<" images ) assigned to img3d: "<<tag<<endl;
}

/************************************************************************/
// Member functions
/************************************************************************/

void img3d::fillRegion2D(int *iniPxl, int *finPxl, int wSize,int foreground,
			 int background,string mode,int *id_vec)
{ /* Recognize regions encosing iniPxl and finPxl. 
     Calls img::fillRegion */
  
  int i,ndigit,zWidth,refPxl[2];
  double mRef[2],tempRefPxl[2];
  string sdum; stringstream ss;
  
  zWidth=id_vec[1]-id_vec[0];
  mRef[0]=(double)(finPxl[0]-iniPxl[0])/(zWidth);
  mRef[1]=(double)(finPxl[1]-iniPxl[1])/(zWidth);
  for (i=0;i<2;i++) {
    refPxl[i]=iniPxl[i];
    tempRefPxl[i]=(double)iniPxl[i];
  }
  
  // count number of digits to use
  ss << id_vec[1]; sdum=ss.str(); ndigit=sdum.length(); 
  for (i=id_vec[0];i<=id_vec[1];i++) {
    if ((wSize> (int)img_stack[i].rows -1 ) ||(wSize> (int)img_stack[i].columns -1))
      error_dodri("wsize must be smaller than image size.");
    cout << "  ID= " <<setw(ndigit)<<i <<" tag= " << img_stack[i].tag;
    refPxl[0]=(int)round(tempRefPxl[0]); refPxl[1]=(int)round(tempRefPxl[1]);
    img_stack[i].fillRegion(refPxl,wSize,foreground,background,mode,"");
    tempRefPxl[0]+=mRef[0];tempRefPxl[1]+=mRef[1];
    cout << endl;
  }
}

/************************************************************************/
void img3d::gaussBlur2D(double *s,int *w,  unsigned char pxlCut, int *id_vec)
{ /* 3D gaussian blur filter for thresholding. Calls 2D gaussBlur */

  int istk,ndigit;
  string sdum;
  stringstream ss;
  ss<< id_vec[1]; sdum=ss.str(); ndigit=sdum.length(); 

  for (istk=id_vec[0];istk<=id_vec[1];istk++) {
    cout << "  ID= " <<setw(ndigit)<<istk <<" tag= " << img_stack[istk].tag;
    img_stack[istk].gaussBlur(s,w,pxlCut,"");
    cout<<endl; 
  }
}

/************************************************************************/
int img3d::imgTag2Id(string sdum)
{ /* Find id number for the img with tag sdum.
     Returns: id number or -1 if tag not found 
  */
  
  int i,j;
  vector<img>::iterator it_img;
  j=-1;
  for (i=0;i<nStack;i++) {
    if (img_stack[i].tag==sdum) {j=i; break;}
  }
  return j;
}

/************************************************************************/
void img3d::invert2D(int *id_vec)
{  /* Inverts pixel intensity of image. */ 
  
  stringstream ss;  string sdum;
  int i,ndigit;
  
  // count number of digits to use
  ss << id_vec[1]; sdum=ss.str(); ndigit=sdum.length(); 
  
  for (i=id_vec[0];i<=id_vec[1];i++) {
    cout << "  ID= " <<setw(ndigit)<<i ; //<<" tag= " << img_stack[i].tag;
    img_stack[i].invert("");
    cout<<endl; 
  }
}

/************************************************************************/
void img3d::merge(string *merge_tag,string mode,int *id_vec)
{ /* Calls image merge. Merges image in merge_tag into image 
     specified by vimg3d *id_vec. 
     Resulting image saved in merge_tag. */

  int i,ndigit;
  vector<img>::iterator it_img1,it_img1_ini, it_img1_fin; //merged tag, saved
  vector<img>::iterator it_img; //this tag
  stringstream ss;
  string sdum;
  ss<< id_vec[1]; sdum=ss.str(); ndigit=sdum.length(); 
  
  it_img1_ini=findImg(merge_tag[0]); //img from another stack
  it_img1_fin=findImg(merge_tag[1]);
  it_img1=it_img1_ini;    

  for (i=id_vec[0];i<=id_vec[1];i++) {
    cout << "  ID= " <<setw(ndigit)<<i <<" tag= " << img_stack[i].tag;
    it_img=findImg(img_stack[i].tag); //current stack
    mergeImg(&(*it_img),&(*it_img1),mode,"");
    ++it_img1;     if (it_img1>it_img1_fin) break;
  } 
  return;  
}

/************************************************************************/
void img3d::removeDebris2D(int sizeCut, unsigned char *pxlCut,
			   unsigned char fillPxl,string mode,int *id_vec)
{ /* Image noise removal. Calls img::removeDebris. */
  
  string sdum; stringstream ss;
  int i,ndigit;
  
  // count number of digits to use
  ss << id_vec[1]; sdum=ss.str(); ndigit=sdum.length(); 
  for (i=id_vec[0];i<=id_vec[1];i++) {
    cout << "  ID= " <<setw(ndigit)<<i <<" tag= " << img_stack[i].tag;
    img_stack[i].removeDebris(sizeCut,pxlCut,fillPxl,mode,"");
    cout<<endl; 
  }
}

/************************************************************************/
void img3d::shift(int *dr, unsigned char marginPxl,int *id_vec)
{ /* Perform shifting. pxl_gray is overwritten. */
  
  stringstream ss;  string sdum;
  int i,ndigit;

  // count number of digits to use
  ss << id_vec[1]; sdum=ss.str(); ndigit=sdum.length(); 
  
  for (i=id_vec[0];i<=id_vec[1];i++) {
    cout << "  ID= " <<setw(ndigit)<<i <<" tag= " << img_stack[i].tag;
    img_stack[i].shift(dr,marginPxl,"","");
    cout<<endl; 
  }
}

/************************************************************************/
void img3d::threshold(string kind,string mode, double *cutoff, int flag,
		      int *id_vec)
{ /* Perform thresholding. pxl_gray is overwritten.*/
  
  string sdum; stringstream ss;
  int i,ndigit;
  
  // count number of digits to use
  ss << id_vec[1]; sdum=ss.str(); ndigit=sdum.length(); 
  for (i=id_vec[0];i<=id_vec[1];i++) {
    cout << "  ID= " <<setw(ndigit)<<i <<" tag= " << img_stack[i].tag;
    img_stack[i].threshold(kind,mode,cutoff,flag,"");
    cout<<endl; 
  }
}

/************************************************************************/
void img3d::thresholdAuto2D(double cutoff,int nbin,int *id_vec)
{ /* Auto filter. */
  
  string sdum; stringstream ss;
  int i,ndigit;
  
  // Count number of digits to use
  ss << id_vec[1]; sdum=ss.str(); ndigit=sdum.length(); 
  for (i=id_vec[0];i<=id_vec[1];i++) {
    cout << "  ID= " <<setw(ndigit)<<i <<" tag= " << img_stack[i].tag;
    img_stack[i].thresholdAuto(cutoff,nbin,"");
  }
}

/************************************************************************/
void img3d::thresholdHisto2D(double *cutoff, int nbin, int *id_vec)
{ /* Histogram filter. */
  
  string sdum; stringstream ss;
  int i,ndigit;

  // Count number of digits to use
  ss << id_vec[1]; sdum=ss.str(); ndigit=sdum.length(); 
  for (i=id_vec[0];i<=id_vec[1];i++) {
    cout << "  ID= " <<setw(ndigit)<<i <<" tag= " << img_stack[i].tag;
    img_stack[i].thresholdHisto(cutoff,nbin,"");
  }
}

/************************************************************************/
void img3d::writeEmap(double *dr,string fname,string outMode,int *id_vec)
{ /* write pxl_gray into EMDB file.
     dr[3]: pixel size   */
  
  // use nubmer of rows/columns for the first image
  int rows0=img_stack[0].rows, columns0=img_stack[0].columns;
  int i,j,istack,npixel, idum,i0,j0;
  char map[4]="MAP", machst[4]="DA"; // DA=0x44,0x41 (little endian)
  string sdum="#DODRI IMG3d TAG: "+tag;
  char label[80]=""; //="::::DODRI-OUTPUT:::: ::::";
  float adum, fstat[4]; // min/max/avg/std of pixel density
  double rdum0,rdum1;
  int nSlice=id_vec[1]-id_vec[0];
  
  string ofname=addFileNameExt(fname,".mrc"); // add filename extension
  ofstream ofs(ofname.c_str(), std::ios::out | std::ios::binary);
  
  strcpy(label,sdum.c_str());
  cout<<"  Pixel dimensions used for EMAP: ";
  for (i=0;i<3;i++) cout<<" "<<dr[i];
  cout<<endl;
  
  // find min/max/avg/std of pixel density
  fstat[0]=9999.; fstat[1]=-9999.; rdum0=rdum1=0.;
  npixel=0;
  for (istack=id_vec[0];istack<id_vec[1];istack++) {
    if (rows0!=(int)img_stack[istack].rows) 
      error_dodri("Unequal number of rows in img3d. Cannot write emap.");
    if (columns0!=(int)img_stack[istack].columns) 
      error_dodri("Unequal number of columns in img3d. Cannot write emap.");
    for (i=0;i<rows0;i++) {
      for (j=0;j<columns0;j++) {
	if (img_stack[istack].pxl_gray[i][j]>0) { // non-zero pixel
	  adum=(float)(img_stack[istack].pxl_gray[i][j]);
	  if (adum<fstat[0]) fstat[0]=adum;
	  if (adum>fstat[1]) fstat[1]=adum;
	  rdum0+=(double)adum; rdum1+=((double)adum*(double)adum); 
	  ++npixel;
	}
      }
    }
  } // for (istack=0;istack<nStack;istack++) {
  rdum0/=(double)npixel; rdum1/=(double)npixel;
  rdum1-=(rdum0*rdum0);
  if (rdum1<TOL) {    stringstream ss; ss<< rdum1;  }
  fstat[2]=(float)rdum0; fstat[3]=(float)(sqrt(rdum1));
  
  //  // find the largest number of pixels in one dir
  //  npixel=(rows0>columns0)?rows0:columns0;
  //  npixel=(npixel>nstack)?npixel:nstack; 
  // Write header
  // NC,NR,NS (1-3): all should be equal by convention
  ofs.write((const char*) &columns0, sizeof(columns0));
  ofs.write((const char*) &rows0, sizeof(rows0));
  ofs.write((const char*) &nSlice, sizeof(nSlice));
  // mode (4)
  idum=2; ofs.write((const char*) &idum, sizeof(idum));
  
  //// {NC,NR,NS}START (5-7)
  //idum=-1*columns0/2; ofs.write((const char*) &idum, sizeof(idum));
  //idum=-1*rows0/2; ofs.write((const char*) &idum, sizeof(idum));
  //idum=-1*nSlice/2; ofs.write((const char*) &idum, sizeof(idum));
  idum=0;
  ofs.write((const char*) &idum, sizeof(idum));
  ofs.write((const char*) &idum, sizeof(idum));
  ofs.write((const char*) &idum, sizeof(idum));

  // NX,NY,NZ (8-10)
  ofs.write((const char*) &columns0, sizeof(columns0));
  ofs.write((const char*) &rows0, sizeof(rows0));
  ofs.write((const char*) &nSlice, sizeof(nSlice));
  // {X,Y,Z}_LENGTH (11-13)
  adum=(float)dr[0]*(float)rows0;
  ofs.write((const char*) &adum, sizeof(adum));
  adum=(float)dr[1]*(float)columns0;
  ofs.write((const char*) &adum, sizeof(adum));
  adum=fabs(dr[2])*(float)nSlice;
  ofs.write((const char*) &adum, sizeof(adum));
  // ALPHA,BETA,GAMMA (14-16)
  adum=90.0; 
  ofs.write((const char*) &adum, sizeof(adum));
  ofs.write((const char*) &adum, sizeof(adum));
  ofs.write((const char*) &adum, sizeof(adum));
  // MAP{C,R,S} (17-19)
  idum=1; ofs.write((const char*) &idum, sizeof(idum));
  idum=2; ofs.write((const char*) &idum, sizeof(idum));
  idum=3; ofs.write((const char*) &idum, sizeof(idum));
  // AMIN,AMAX,AMEAN (20-22)
  ofs.write((const char*) &fstat[0], sizeof(fstat[0]));
  ofs.write((const char*) &fstat[1], sizeof(fstat[1]));
  ofs.write((const char*) &fstat[2], sizeof(fstat[2]));
  // ISPG (23)
  idum=1; ofs.write((const char*) &idum, sizeof(idum));
  // NSYMBT, LSKFLG (24,25)
  idum=0; ofs.write((const char*) &idum, sizeof(idum));
  idum=0; ofs.write((const char*) &idum, sizeof(idum));
  // SKWMAT,SKWTRN (word 26-37)
  adum=0.;
  for (i=0;i<12;i++) ofs.write((const char*) &adum, sizeof(adum));
  // EXTRA (word 38-52): Later put comments
  for (i=0;i<15;i++) ofs.write((const char*) &adum, sizeof(adum));
  // MAP,MACHST (53,54)
  ofs.write((const char*) &map, sizeof(map));
  ofs.write((const char*) &machst, sizeof(machst));
  // RMS (55)
  ofs.write((const char*) &fstat[3], sizeof(fstat[3]));
  // NLABEL (put only first label)
  idum=1; ofs.write((const char*) &idum, sizeof(idum));
  // LABEL: (57-256), 200 words
  // LABEL_1 (first label; each label is 80 chars (20 words)
  for (i=0;i<80;i++) ofs.write((const char*) &label[i], sizeof(label[i]));
  adum=0.;
  for (i=0;i<180;i++) ofs.write((const char*) &adum, sizeof(adum));
  
  // Write voxel data
  for (istack=id_vec[0];istack<=id_vec[1];istack++) {
    for (i=0;i<rows0;i++){
      idum=(dr[2]>0)?istack:(nStack-1-istack); // determine z-direction
      for (j=0;j<columns0;j++){
	i0=(outMode=="verbose")?i:(rows0-1-i);
	j0=j; //(outMode=="verbose")?j:(columns0-1-j);
	adum=(float)(img_stack[istack].pxl_gray[i0][j0]); 
	ofs.write((const char*) &adum, sizeof(adum));
      }
    }
  }
  ofs.close();
  cout <<"  Writing EMAP to "<<ofname<<endl;
  cout<<"  Min/Max/Avg/Std of pixel intensities: ";
  for (i=0;i<4;i++) cout<<setw(7)<<setprecision(4)<<fstat[i]<<" ";
  cout<<endl;

  return;
}

/************************************************************************/
void img3d::writeImage(string fname,string fmt,int nColor,
		       double *colorBg,double *colorScale,int *id_vec)
{ /* Write img3d into a set of images. */ 
  
  string ofname,sdum; stringstream ss;
  int i,ndigit;
  
  // count number of digits to use
  ss << id_vec[1]; sdum=ss.str(); ndigit=sdum.length(); 
  
  for (i=id_vec[0];i<=id_vec[1];i++) {
    ss.str(""); ss.clear();
    ss<<fname<<"_"<<setfill('0')<<setw(ndigit)<<i; //<<"."<<fmt;
    ofname=ss.str();
    cout << "  ID= " <<setw(ndigit)<<i <<"  tag= " << img_stack[i].tag;
    img_stack[i].writeImage(ofname,fmt,nColor,colorBg,colorScale);
    cout<<endl; 
  }
}

/************************************************************************/
// Non-member function
/************************************************************************/

void checkDuplicateImg3dTag(string tag)
{ // In vimg3d, find iterator for the element with tag
  vector<img3d>::iterator it_img3d;
  for (it_img3d=vimg3d.begin();it_img3d!=vimg3d.end();++it_img3d) {
    if ( it_img3d->tag==tag) error_dodri("Duplicate img3d tag: " + tag);
  }
  return;
}

/************************************************************************/
vector<img3d>::iterator findImg3d(string tag)
{ /* In vimg3d, find iterator for the element with tag */
  
  vector<img3d>::iterator it_img3d;
  for (it_img3d=vimg3d.begin();it_img3d!=vimg3d.end();++it_img3d) {
    if ( (*it_img3d).tag==tag) break;
  }
  return it_img3d;
}
