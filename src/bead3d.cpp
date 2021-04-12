# include "dodri.h"

/************************************************************************/
// Constructors
/************************************************************************/

bead3d::bead3d(img3d *img3d0, string tag0, int iDiam, string tag1,int dir)
{ /* Local pixel area-based constructor */
  int i,ndigit;
  vector<img>::iterator it_img;
  vector<bead>::iterator it_bead;
  string img_tag,bead_tag,sdum,save_img_tag="",start_tag,end_tag;
  stringstream ss;
  
  nStack=img3d0->nStack; img3d_tag=img3d0->tag;
  tag=tag0;
  
  ndigit=nDigit(nStack);

  for (i=0;i<nStack;i++) {
    ss.str(""); ss.clear();
    ss<< tag<<i; bead_tag=ss.str();
    if (tag1 != "") {
      ss.str(""); ss.clear();
      ss<< tag1<<i; save_img_tag=ss.str();
      if (i == 0) start_tag = save_img_tag;
      if (i == nStack-1) end_tag = save_img_tag;
    }

    img_tag=img3d0->img_stack[i].tag;
    it_img=findImg(img_tag);
    cout<<"  ID= "<<setw(ndigit)<<std::right<<i<<" tag= "<<bead_tag;
    vbead.push_back( bead(&(*it_img),bead_tag,iDiam,save_img_tag,dir)); 
    it_bead=findBead(&vbead,bead_tag);
    bead_stack.push_back(*it_bead);
  }
  cout<<"  "<<nStack<<" images used for assigning beads."<<endl;
}

/************************************************************************/
bead3d::bead3d(img3d *img3d0,string tag0,int iDiam,string tag1,
	       int verbosity,int *id_vec)
{ /* Directional area-based bead assignment
     Generate temporary beads based on img0 for scanning 4 directions,
     combine them into a new bead, then erase temporary beads */

  int i,ndigit;
  vector<img>::iterator it_img;
  vector<bead>::iterator it_bead;
  string img_tag,bead_tag,sdum,save_img_tag="",start_tag,end_tag;
  stringstream ss;

  nStack=id_vec[1]-id_vec[0]+1;
  img3d_tag=img3d0->tag;
  tag=tag0;  
  ndigit=nDigit(nStack);

  for (i=id_vec[0];i<=id_vec[1];i++) {
    ss.str(""); ss.clear();
    ss<< tag<<i; bead_tag=ss.str();
    if (tag1 != "") {
      ss.str(""); ss.clear();
      ss<< tag1<<i; save_img_tag=ss.str();
      if (i == id_vec[0]) start_tag = save_img_tag;
      if (i == id_vec[1]) end_tag = save_img_tag;
    }

    img_tag=img3d0->img_stack[i].tag;
    it_img=findImg(img_tag);
    cout<<"  ID= "<<setw(ndigit)<<std::right<<i<<" tag= "<<bead_tag;

    vbead.push_back(bead(&(*it_img),bead_tag,iDiam,save_img_tag,
			 -1.,-1.,-1,-1,-1,-1.,verbosity,"area"));
    it_bead=findBead(&vbead,bead_tag);
    bead_stack.push_back(*it_bead);
  }
  cout<<"  "<<nStack<<" images used for assigning beads."<<endl;
}

/************************************************************************/
bead3d::bead3d(string tag, string ifname, string image3d_tag0)
{ /* Constructs a bead from the data in a .dat file with the name ifname */
  
  this->tag = tag;
  readData(ifname,image3d_tag0);
}

/************************************************************************/
// Member functions
/************************************************************************/

int bead3d::beadTag2Id(string sdum)
{ /* Find id number for the bead with tag sdum.
     Returns: id number or -1 if tag not found */

  int i,j;
  vector<bead>::iterator it_bead;
  j=-1;
  for (i=0;i<nStack;i++) {
    if (bead_stack[i].tag==sdum) {j=i; break;}
  }
  return j;
}

/************************************************************************/
void bead3d::orient(double tol, string ax0,int iorie0, int *id_vec)
{ /* Orient BOC based on offset from slice of longest axis with x-axis. 
     Modifies every 2D bead based on ^ offset. 
     
     If beads are moved outside of starting bounds, 
     BOC will be repositioned.    */
  
  int i,j,istk,maxBeadCount,xstk=-1,ystk=-1,iorie=-1;
  double max_stk,max_dist,dist,sig,offset;
  double *bm=new double[2],**pos;
  
  // for position adjustment within quadrant I
  int xmove,ymove;
  double xmin,ymin,xmax,ymax,dx=0,dy=0,cm_temp[2];
  
  cout<<"  orient: tolerance="<<tol<<" axis="<<ax0<<endl; 
  
  // get largest maxBead for pos sizing
  maxBeadCount=TOL; 
  for (istk=id_vec[0];istk<=id_vec[1];istk++) 
    maxBeadCount=((int)bead_stack[istk].maxBead>maxBeadCount)?
      ((int)bead_stack[istk].maxBead):(maxBeadCount);

  pos=new double*[2];
  for (j=0;j<2;j++) pos[j]=new double[maxBeadCount]; 
  
  // Find istk with longest axis. Save to iorie. 
  max_stk=TOL; // max dist overall
  for (istk=id_vec[0];istk<=id_vec[1];istk++) {
    if (!bead_stack[istk].nBead) continue; // skip if no beads
    
    max_dist=TOL; // max distance between any point and bead centroid
    for (i=0;i<bead_stack[istk].nBead;i++) {
      for (j=0;j<2;j++) pos[j][i]=bead_stack[istk].cm[i][j];
    }
    for (j=0;j<2;j++)
      getavgw(&bm[j],&sig,pos[j],bead_stack[istk].beadMass,bead_stack[istk].nBead); 

    for (i=0;i<bead_stack[istk].nBead;i++) {
      dist=getDistPos(bead_stack[istk].cm[i],bm,0);
      if (dist>max_dist) max_dist=dist;
    }
    if (max_dist>max_stk) { // istk w/ longest axis
      iorie=istk;
      max_stk=max_dist;
    }
  }
  
  iorie=(iorie0==-1)?(iorie):(iorie0); 
  cout<<"          orienting with slice:"<<iorie<<endl; 
  
  do {
    offset=bead_stack[iorie].getAngleOffset(); // get offset of iorie
    for (istk=id_vec[0];istk<=id_vec[1];istk++) 
      bead_stack[istk].rotate(offset*rad2degree,bm); // rotate all by offset
  } while (abs(offset)>tol); 
  
  if (ax0=="y") {
    for (istk=id_vec[0];istk<=id_vec[1];istk++) 
      bead_stack[istk].rotate((PI/2.)*rad2degree,bm); // orient to y-axis 
  }

  // If beads were moved outside of QI, shift beads to place in QI
  xmove=0; ymove=0;
  xmin=ymin=1./TOL; xmax=ymax=TOL; 
  for (istk=id_vec[0];istk<=id_vec[1];istk++) {
    for (i=0;i<bead_stack[istk].nBead;i++) {
      if (bead_stack[istk].cm[i][0]<0) {
	xmove=1; 
	xmin=(bead_stack[istk].cm[i][0]<xmin)?(bead_stack[istk].cm[i][0]):(xmin); 
      }
      if (bead_stack[istk].cm[i][0]>bead_stack[istk].columns) {
	xmove=1; xstk=istk; 
	xmax=(bead_stack[istk].cm[i][0]>xmax)?(bead_stack[istk].cm[i][0]):(xmax); 	
      }
      if (bead_stack[istk].cm[i][1]<0) {
	ymove=1;
	ymin=(bead_stack[istk].cm[i][1]<ymin)?(bead_stack[istk].cm[i][1]):(ymin); 
      }
      if (bead_stack[istk].cm[i][1]>bead_stack[istk].rows) {
	ymove=1; ystk=istk; 
	ymax=(bead_stack[istk].cm[i][1]>ymax)?(bead_stack[istk].cm[i][1]):(ymax); 
      }
    }
  } //  for (istk=id_vec[0];istk<=id_vec[1];istk++) {

  if (xmove) { // translate BOC in x-dir by xmin
    cout<<" x-dir move activated with: xmin "<<xmin<<" and xmax:"<< xmax<<endl; 
    if (xmin!=1./TOL && xmax!=TOL)
      error_dodri("bead3d::orient: error during rotation.");
    if (xmin!=1./TOL) dx=-1.*xmin;
    else if (xmax!=TOL)  dx=-1.*(xmax-(double)bead_stack[xstk].columns);
  }
  if (ymove) { // translate BOC in y-dir by ymin 
    cout<<" y-dir move activated"<<endl;    
    if (ymin!=1./TOL && ymax!=TOL)
      error_dodri("bead3d::orient: error during rotation.");
    if (ymin!=1./TOL)  dy=-1.*ymin;
    else if (ymax!=TOL)  dy=-1.*(ymax-(double)bead_stack[ystk].rows);
  }

  // use same dx,dy shift for al slices
  for (istk=id_vec[0];istk<=id_vec[1];istk++) {  
    for (i=0;i<bead_stack[istk].nBead;i++) {
      cm_temp[0]=bead_stack[istk].cm[i][0]+dx;
      cm_temp[1]=bead_stack[istk].cm[i][1]+dy;       
      bead_stack[istk].moveBeadPos(i,cm_temp);
    }
  }
  
  // clear memory 
  for (j=0;j<2;j++) delete [] pos[j];
  delete [] pos; delete [] bm; 

  return; 
}

/************************************************************************/
void bead3d::readData(string fname,string image3d_tag0)
{ /* Read a .dat file made from bead::writeData to expedite beads. */
  
  stringstream ss;
  string image3d_tag, image_tag,line, ntag,ofname; 
  int i,ivec[1], istk = 0, ndigit;
  ifstream fin;
  vector<bead>::iterator it_bead;
  bool i3d_tag, nstk, flag;
  vector<img3d>::iterator it_img3d; 
  flag = i3d_tag = nstk = false;

  ofname=addFileNameExt(fname,".dat"); 
  cout<<"  read:: Read data from: " << ofname<<endl;
  if (img3d_tag!="") {
    cout <<"  "<<std::right<<
      "WARNING: bead3d " << tag << " is being overwritten.\n";
    for (i=0; i < nStack; i++) {
      delete[] bead_stack[i].beadMass;
      for (i = 0; i < (int)bead_stack[i].maxBead; i++)
	delete[] bead_stack[i].cm[i];
      delete[] bead_stack[i].cm;
    
      bead_stack[i].bead_cell.clear();
      bead_stack[i].cell_bead.clear();
    }
    bead_stack.clear();
  }

  nStack = 0;
  fin.open(ofname.c_str(), ios::in);
  if (fin.fail()) error_dodri("dodri: Invalid bead3d file");

  while (!flag && fin.good()) {
    flag = i3d_tag && nstk;
    if (flag) break;
    
    if (line.find("Stack: ")!=std::string::npos)
      error_dodri("READ:: Incomplete bead3d header");
    getline(fin, line);
    if (line == "") continue;
    ss.str(line);
    
    if (!i3d_tag && line.find("img3d_tag=")!=string::npos) {
      if (image3d_tag0=="")
	image3d_tag=getCommOptString(ss,"img3d_tag=","");
      else image3d_tag=image3d_tag0;
      if (image3d_tag!="") {
	img3d_tag = image3d_tag;
	i3d_tag = true;
	continue;
      }
      else error_dodri("READ:: missing image3d tag");
    }

    if (!nstk && line.find("nStack=")!=string::npos) {
      getCommOptInt(ss, "nStack=", ivec, 1, -1);
      if (ivec[0] != -1) {
	nStack = ivec[0];
	nstk = true;
	continue;
      }
      else error_dodri("READ:: invalid or missing nStack");
    }
  }
  
  ndigit = nDigit(nStack);
  it_img3d=findImg3d(img3d_tag);
  if (it_img3d==vimg3d.end()) error_dodri("READ:: img3d tag "+img3d_tag+" not found");
  
  while (fin.good()) {
    getline(fin, line);
    if (line == "") continue;
    ss.str(line);
    if (line.find("Stack:")!=string::npos) {      
      ntag = tag.c_str() + to_string(istk);
      image_tag=(*it_img3d).img_stack[istk].tag;
      cout<<"  ID= "<<setw(ndigit)<<std::right<<istk<<" tag= "<<ntag;      
      vbead.push_back(bead(ntag, ofname,image_tag, &fin));
      it_bead=findBead(&vbead,ntag);
      bead_stack.push_back(*it_bead);
      istk++;
    }
  }
  return; 
}

/************************************************************************/
void bead3d::rotate(double theta,int *id_vec)
{ /* Rotate by theta (in deg). Use bcm of slice with longest axis for 
     turning point.    */ 

  int i,j,istk,maxBeadCount,ndigit;
  int xstk,ystk=-1,xmove,ymove;  
  double *bm=new double[2],**pos;
  double xmin,ymin,xmax,ymax,dx=0,dy=0,cm_temp[2];

  // get largest maxBead for pos sizing
  maxBeadCount=TOL; 
  for (istk=id_vec[0];istk<=id_vec[1];istk++) 
    maxBeadCount=((int)bead_stack[istk].maxBead>maxBeadCount)?
      ((int)bead_stack[istk].maxBead):(maxBeadCount);

  pos=new double*[2];
  for (j=0;j<2;j++) pos[j]=new double[maxBeadCount]; 
    
  ndigit=nDigit(nStack);
  for (istk=id_vec[0];istk<=id_vec[1];istk++) {
    cout << "  ID= " <<setw(ndigit)<<istk;
    bead_stack[istk].rotate(theta,bm);
  }

  // If beads were moved outside of QI, shift beads to place in QI
  // for position adjustment within quadrant I  
  xmove=0; ymove=0;
  xmin=ymin=1./TOL; xmax=ymax=TOL; 
  for (istk=id_vec[0];istk<=id_vec[1];istk++) {
    for (i=0;i<bead_stack[istk].nBead;i++) {
      if (bead_stack[istk].cm[i][0]<0) {
	xmove=1; 
	xmin=(bead_stack[istk].cm[i][0]<xmin)?(bead_stack[istk].cm[i][0]):(xmin); 
      }
      if (bead_stack[istk].cm[i][0]>bead_stack[istk].columns) {
	xmove=1; xstk=istk; 
	xmax=(bead_stack[istk].cm[i][0]>xmax)?(bead_stack[istk].cm[i][0]):(xmax); 	
      }
      if (bead_stack[istk].cm[i][1]<0) {
	ymove=1;
	ymin=(bead_stack[istk].cm[i][1]<ymin)?(bead_stack[istk].cm[i][1]):(ymin); 
      }
      if (bead_stack[istk].cm[i][1]>bead_stack[istk].rows) {
	ymove=1; ystk=istk; 
	ymax=(bead_stack[istk].cm[i][1]>ymax)?(bead_stack[istk].cm[i][1]):(ymax); 
      }
    }
  } //  for (istk=id_vec[0];istk<=id_vec[1];istk++) {

  if (xmove) { // translate BOC in x-dir by xmin
    cout<<" xmove activated with: xmin "<<xmin<<" and xmax:"<< xmax<<endl; 
    if (xmin!=1./TOL && xmax!=TOL)
      error_dodri("bead3d::orient: error during rotation.");
    if (xmin!=1./TOL) dx=-1.*xmin;
    else if (xmax!=TOL)  dx=-1.*(xmax-(double)bead_stack[xstk].columns);
  }
  if (ymove) { // translate BOC in y-dir by ymin 
    cout<<" ymove activated"<<endl;    
    if (ymin!=1./TOL && ymax!=TOL)
      error_dodri("bead3d::orient: error during rotation.");
    if (ymin!=1./TOL)  dy=-1.*ymin;
    else if (ymax!=TOL)  dy=-1.*(ymax-(double)bead_stack[ystk].rows);
  }

  // use same dx,dy shift for al slices
  for (istk=id_vec[0];istk<=id_vec[1];istk++) {  
    for (i=0;i<bead_stack[istk].nBead;i++) {
      cm_temp[0]=bead_stack[istk].cm[i][0]+dx;
      cm_temp[1]=bead_stack[istk].cm[i][1]+dy;       
      bead_stack[istk].moveBeadPos(i,cm_temp);
    }
  }
  return; 
}

/************************************************************************/
void bead3d::setz(double dz, string origin)
{ /* origin: zCoord=0 for 0th bead ("ini"); (nStack-1)-th bead ("fin");
     at (nStack-1)/2) ("mid")
     dz: adds to previous bead in the stack. So, for "fin" and dz>0,
     0-th bead has negative zCoord.  */
  
  int i,i0;
  if (origin=="ini") {
    for (i=0;i<nStack;i++) bead_stack[i].zCoord=dz*(double)i;
  }
  else if (origin=="fin") {
    for (i=(nStack-1);i>=0;i--) {
      i0=(i-nStack+1);
      bead_stack[i].zCoord=dz*(double)i0;
    }
  }
  else if (origin=="mid") {
    for (i=0;i<nStack;i++) 
      bead_stack[i].zCoord=dz*((double)i-0.5*(double)(nStack-1));
  }
  else 
    error_dodri("  Unrecognized origin keyword: "+origin);
  cout<<"  setz: dz= "<<dz<<" origin= "<<origin<<
    " zCoord range: "<<bead_stack[0].zCoord<<" to "
      <<bead_stack[nStack-1].zCoord<<endl;
}

/************************************************************************/
void bead3d::writeCor(string fname, string outMode, int *id_vec)
{ /* Write bead coordinate file. */
  
  FILE *ff;
  string ofname,segname,sdum;
  stringstream ss;
  
  double x1, y1, z1, massMin, massMax, occupancy;
  double massMin0,massMax0;
  int i,j,imin, imax,iBead=0,nBead0=0,ndigit;
  multimap<int,int>::iterator it;
  
  // determine filename extension
  ofname = addFileNameExt(fname,".cor");
  cout<<"  Write cor for bead3d IDs "<<id_vec[0]<<" - "<<
    id_vec[1]<<" to "<<ofname<<endl;
  ff=fopen(ofname.c_str(),"w");
  
  // total number of beads in selection
  for (i=id_vec[0];i<=id_vec[1];i++) nBead0+=bead_stack[i].nBead;
  
  // find min/max mass: normalized mass entered as occupancy (col 55-60)
  massMin=1e308; massMax=-1;
  for (i=id_vec[0];i<=id_vec[1];i++) {
    MinMax(bead_stack[i].nBead, bead_stack[i].beadMass, 
	   &imin, &imax, &massMin0, &massMax0);
    massMin=(massMin0<massMin)?massMin0:massMin;
    massMax=(massMax0>massMax)?massMax0:massMax;
  }
  
  // use extended format for more than 99999 atoms
  if (nBead0>99999) fprintf(ff,"%10d%s\n",nBead0,"  EXT");
  else fprintf(ff,"%d\n",nBead0);
  
  ndigit=nDigit(nStack);
  
  for (i=id_vec[0];i<=id_vec[1];i++) {
    ss.str(""); ss.clear();
    if (ndigit<5) ss<<i;
    else ss<<hex<<i;
    segname=ss.str();
    for (j=0;j<bead_stack[i].nBead;j++) {
      if (outMode=="verbose") {
	x1=bead_stack[i].cm[j][0]; y1=bead_stack[i].cm[j][1]; 
	z1=bead_stack[i].zCoord;
      }
      else { // outMode=="image"
	y1=((double)(bead_stack[i].rows)-bead_stack[i].cm[j][0]); 
	x1=bead_stack[i].cm[j][1]; z1=bead_stack[i].zCoord;
      } 
      occupancy=bead_stack[i].beadMass[j]/massMax;
      
      /* From coorio.src:
	 fm2='(2I5,1X,A4,1X,A4,3F10.5,1X,A4,1X,A4,F10.5)'
	 I,IRES,RES(IRES),ATYPE(I),Xyz(1,I),xYz(2,I),xyZ(3,I),SID,RID,WMAIN(I)
      */
      if (nBead0>99999)  // use extended format
	fprintf (ff,"%10d%10d  %-8s  %-8s%20.10f%20.10f%20.10f  %-8s  %-8d%20.10f\n",
		 ++iBead,j,"BD3D","O",x1,y1,z1,segname.c_str(),j,occupancy);
      else 
	//    fprintf (ff,"%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %4s %-4d%10.5f\n",
	fprintf (ff,"%5d%5d %-4s %-4s%9.4f %9.4f %9.4f  %4s %-4d%10.5f\n",
		 ++iBead,j,"BD3D","O",x1,y1,z1,segname.c_str(),j,occupancy);
    } // for (j=0;j<bead_stack[i].nBead;j++) {
  } // for (i=id_vec[0];i<=id_vec[1];i++) {
  fclose(ff); return; 
}

/************************************************************************/
void bead3d::writeData(string fname, int *id_vec)
{ /* Writes .dat file of bead coordinates to open in DODRI. */

  FILE *dat;
  string datname;
  int i,ndigit;
  datname = addFileNameExt(fname,".dat");
  dat = fopen(datname.c_str(), "w");
  cout<<"  Write bead data to "<<datname<<endl;
  
  fprintf(dat, "nStack= %d\nimg3d_tag= %s\n", id_vec[1]-id_vec[0]+1,
	  img3d_tag.c_str());
  
  ndigit=nDigit(nStack);
  for (i = id_vec[0]; i <= id_vec[1]; i++) {
    cout << "  ID= " <<setw(ndigit)<<i<<" tag= " << bead_stack[i].tag<<endl; 
    fprintf(dat, "\nStack: %d\n\n", i-id_vec[0]);
    bead_stack[i].writeData(fname, dat);
  }
  fclose(dat); return; 
}

/************************************************************************/
void bead3d::writeImage(int beadL, int nColor, double *col_bg,double *col_bead,
			string fname, string fmt,int *id_vec)
{ /* Write bead3d into a set of images. */ 

  string ofname,sdum; stringstream ss;
  int i,ndigit,k;

  // count number of digits to use
  ss << id_vec[1]; sdum=ss.str(); ndigit=sdum.length(); 

  for (i=id_vec[0];i<=id_vec[1];i++) {
    ss.str(""); ss.clear();
    ss<<fname<<"_"<<setfill('0')<<setw(ndigit)<<i; 
    ofname=ss.str();
    cout << "  ID= " <<setw(ndigit)<<i <<"  tag= " << bead_stack[i].tag;
    // when beadL=0, use beadDiam for each bead_stack
    k=(beadL==0)?(bead_stack[i].beadDiam):beadL; 
    bead_stack[i].writeImage(k,nColor,col_bg,col_bead,ofname,fmt);
  }
  return;
}

/************************************************************************/
void bead3d::writePsf(string fname, int *id_vec)
{ /* Writes psf file.
     id_vec[2]: begin/end stack numbers
     Use stack number as SEGID. If number of stacks is greater than 1e4,
     use hexadecimal expression to fit it into 4-letter SEGID.*/
  
  int i,j,ndigit,iBead=0,nBead0=0;
  double rdum;
  string ofname,sdum;  stringstream ss;
  string segname,psfname;
  FILE *psf;

  ofname = addFileNameExt(fname,".psf"); 
  cout<<"  Write psf for bead3d IDs "<<id_vec[0]<<" - "<<
    id_vec[1]<<" to "<<ofname<<endl;

  psf=fopen(ofname.c_str(),"w");
  fprintf(psf,"PSF \n\n !NTITLE\n");
  fprintf(psf,"* PSF for %s \n* \n\n",fname.c_str());

  // total number of beads in selection
  for (i=id_vec[0];i<=id_vec[1];i++) nBead0+=bead_stack[i].nBead;
  fprintf(psf,"%6d !NATOM\n",nBead0);

  ndigit=nDigit(nStack);

  for (i=id_vec[0];i<=id_vec[1];i++) {
    ss.str(""); ss.clear();
    if (ndigit<5) ss<<i;
    else ss<<hex<<i;
    segname=ss.str();
    for (j=0;j<bead_stack[i].nBead;j++) {
      rdum=bead_stack[i].beadMass[j];
      fprintf(psf,"%8d %4s %04d BD3D O       0 %14.6f %14.6f %d\n",
	      ++iBead,segname.c_str(),j,0.0,rdum,0);
      /* from psfres.src: fmt01=(I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8)
	 I(8),SEGID(ISEG)(4),RESID(IRES)(4),RES(IRES)(4),
	 ATYPE(I)(4),IAC(I)(4; atomID),CG(I)(14.6;charge),
	 AMASS(I)(14.6;mass),IMOVE(I)(8)
      */
    } // for (j=0;j<bead_stack[i].nBead;j++) {
  } //   for (i=id_vec[0];i<=id_vec[1];i++) {
  fprintf(psf,"\n%8d !NBOND\n",0);
  fclose(psf);  return;
}

/************************************************************************/
// Non-member functions
/************************************************************************/
void checkDuplicateBead3dTag(vector<bead3d> * vbead3d,string tag)
{ // In vbead3d, find iterator for the element with tag
  vector<bead3d>::iterator it_bead3d;
  for (it_bead3d=(*vbead3d).begin();it_bead3d!=(*vbead3d).end();++it_bead3d) {
    if ( (*it_bead3d).tag==tag) error_dodri("Duplicate bead3d tag: " + tag);
  }
  return;
}

/************************************/
vector<bead3d>::iterator findBead3d(vector<bead3d> *vbead3d,string tag)
{ // In vbead3d, find iterator for the element with tag
  vector<bead3d>::iterator it_bead3d;
  for (it_bead3d=(*vbead3d).begin();it_bead3d!=(*vbead3d).end();++it_bead3d) {
    if ( (*it_bead3d).tag==tag) break;
  }
  return it_bead3d;
}
