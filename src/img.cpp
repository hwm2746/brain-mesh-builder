#include "dodri.h"

/************************************************************************/
// Constructors
/************************************************************************/

img::img(string ifname, string tag0) { 
  /* Constructor: read an image file & set global vars. 
     This constructore uses Magick++ API.    */ 
  
  ifstream infile(ifname);
  if (!infile.good()) error_dodri("Input image file does not exist.");
  
  int i, irow, icol, channel0;
  double rdum, *luminosity;
  unsigned int Type; 
  image_name=ifname; tag=tag0;
  
  Image img0(ifname);
  
  //Read Image Attributes
  rows=img0.rows();
  columns=img0.columns();
  BitsPerSample=img0.modulusDepth(); // no. bits req'd to support rgb components
  SamplesPerPixel=img0.colorSpace(); 
  TotalColors = img0.totalColors();  //>255 = multi-channel, 2=binary
  Type = img0.type(); 
  range=pow(2,BitsPerSample); 
  ncolors = SamplesPerPixel;
  if (TotalColors>255 && Type==8) ncolors=SamplesPerPixel*4;   //rgb
  else if (TotalColors>255 && Type>4) ncolors=SamplesPerPixel*3;   //cymk
  
  printf ("  Processing %s\n",ifname.c_str());
  printf ("  Tag = %s\n", tag.c_str());
  printf ("  nCell = %d\n", nCell);
  printf ("  rows,columns = %d %d\n", rows,columns);
  printf ("  BitsPerSample = %d\n", BitsPerSample);
  printf ("  TotalColors = %d\n", TotalColors);
  printf ("  Image Type = %d\n", Type);
  printf ("  Range = %d\n", range);
  printf ("  Number of colors = %d\n", ncolors);
  
  ///////////////////////////////////////////////////////
  // Read pixel values row by row 
  PixelPacket *pixels = img0.getPixels(0, 0, columns, rows); 
  irow=0;icol=0;
  Color color=pixels[columns * irow +icol]; //initialize object
  assert(color.isValid());

  pxl_raw=new unsigned char*[rows];
  pxl_gray=new unsigned char*[rows];
  col_tot=columns*ncolors;
  for (irow=0;irow<(int)rows;irow++) {
    pxl_raw[irow]=new unsigned char[col_tot];
    pxl_gray[irow]=new unsigned char[columns];
  }

  channel0=-1;
  if (ncolors==1)  { // binary or gray
    // choose the first channel with nonzero intensity
    for (irow = 0; irow < (int)rows; irow++) {
      if (channel0>-1) break;
      for (i=0;i<(int)col_tot;i++)  { 
	color = pixels[columns * irow + i];
	if (color.redQuantum()) { channel0=0; break;}
	if (color.greenQuantum()) { channel0=1; break;}
	if (color.blueQuantum()) { channel0=2; break;}
      }
    } // for (irow = 0; irow < (int)rows; irow++) {
    if (channel0==0) cout <<"  Red channel used."<<endl;
    else if (channel0==1) cout <<"  Green channel used."<<endl;
    else if (channel0==2) cout <<"  Blue channel used."<<endl;
    else cout<<"    MESSAGE: Image file has no content. Setting pixel values to 0."<<endl;  
  } // if (ncolors==1)  { // binary or gray
  
  for (irow = 0; irow < (int)rows; irow++) {
    if (ncolors==1)  { // binary or gray
      for (i=0;i<(int)col_tot;i++)  {
	color = pixels[columns * irow + i];
	if (channel0==0)
	  pxl_raw[irow][i]=(unsigned char)round(color.redQuantum() / range);
	else if (channel0==1)
	  pxl_raw[irow][i]=(unsigned char)round(color.greenQuantum() / range);
	else if (channel0==2)
	  pxl_raw[irow][i]=(unsigned char)round(color.blueQuantum() / range);
	else
	  pxl_raw[irow][i]=(unsigned char)0; 
      }
    }
    else if (ncolors<4) { //rgb
      for (i=0;i<(int)col_tot;i++)  {
	color = pixels[columns * irow + i/3];
	pxl_raw[irow][i]=(color.redQuantum() / range) ; 
	i++;	color = pixels[columns * irow + i/3];
	pxl_raw[irow][i]=(color.greenQuantum() / range) ;
	i++;	color = pixels[columns * irow + i/3];
	pxl_raw[irow][i]=(color.blueQuantum() / range) ; 
      }
    }
  } 
  
  // Convert to grayscale
  luminosity=new double[ncolors]; // luminosity values
  if (ncolors==1) luminosity[0]=1.; // single-color
  else { // count first 3 colors (http://en.wikipedia.org/wiki/Luminance_(relative))
    luminosity[0]=0.2126; luminosity[1]=0.7152; luminosity[2]=0.0722;
  }
  
  for (irow = 0; irow < (int)rows; irow++) {
    for (icol=0;icol<(int)columns;icol++) {
      rdum=0.;
      for (i=0;i<ncolors;i++) 
        rdum+=(double)pxl_raw[irow][icol*ncolors+i]*luminosity[i];
      pxl_gray[irow][icol]=(unsigned char)rdum;
    }
  }
  delete [] luminosity; // free up memory
  
  dRow=(double)rows/(double)nCell; // cell size in the row dir
  dCol=(double)columns/(double)nCell; // cell size in the col dir
}

/************************************************************************/
img::img(const img &img0, unsigned char **pxl, string tag0)
{
  /* Copy constructor: creates an image object that is identical to the
     input image except for a new pxl_gray and tag. pxl[rows][columns] must
     have the same dimesion as img0->pxl_gray[rows][columns] (no check
     made)   */
  
  int i,j;
  image_name = img0.image_name;
  tag=tag0;
  rows = img0.rows;
  columns = img0.columns;
  col_tot = img0.col_tot;
  BitsPerSample = img0.BitsPerSample;
  SamplesPerPixel = img0.SamplesPerPixel;
  
  dRow = img0.dRow;
  dCol = img0.dCol;
  ncolors = img0.ncolors;
  range = img0.range;
  // cannot use: pxl_gray = pxl since this doesn't allow copying
  // to different images

  pxl_raw = new unsigned char*[rows];
  for (i=0;i<(int)rows;i++) pxl_raw[i]=new unsigned char[col_tot];
  for (i=0;i<(int)rows;i++) {
    for (j=0;j<(int)col_tot;j++) pxl_raw[i][j]=img0.pxl_raw[i][j];
  }
  pxl_gray=new unsigned char*[rows];
  for (i=0;i<(int)rows;i++) pxl_gray[i]=new unsigned char[columns];
  for (i=0;i<(int)rows;i++) {
    for (j=0;j<(int)columns;j++) pxl_gray[i][j]=pxl[i][j];
  }

  cout<<"    Created "<<tag0<<" based on image "<< img0.tag<<endl;  
}

/************************************************************************/
// Member functions
/************************************************************************/

void img::fillRegion(int *refPxl,int wSize,int foreground,
		     int background,string mode,string tag1)
{ /* Fill region contiguous in color intensity to reference pixel refPxl.
     If any of the pixels in the mask window of size wSize differs in
     intensity from refPxl, that window is marked 'sealed.' In the end,
     refPxl will be in a region surrounded by an edge consisting of
     sealed windows.

     foreground: pixel intensity for filling areas
     mode=region: Fill the region with foreground.  
     mode=edge: Fill the edge with foreground.
     tag1: If non-empty, save output to new img. 
  */
  
  int i,j, *refGrid=new int [2]; // refGrid: grid coordinate for refPxl
  int ref_intensity; 
  unsigned char **pxl_tmp=new unsigned char*[rows];
  char **grid, cflag; // grid: grid of masks.
  int wx,wy,gx,gy; // wx,wy: size of grid
  // grid=0: unchecked, 1: checked, 2: sealed
  int x0,x1,y0,y1; // pixel range 
  vector<pair<int,int>>  edge, edge1,inner;
  vector<pair<int,int>>::iterator it;
  string sdum=(tag1=="")?tag:tag1;
  vector<img>::iterator it_img; 
  
  // for multi-dir search
  
  int idir,flag; 
  unsigned char ***pxl_dir=new unsigned char**[4];
  for (idir=0;idir<4;idir++) pxl_dir[idir]=new unsigned char*[rows];
  for (idir=0;idir<4;idir++) {
    for (i=0;i<(int)rows;i++) pxl_dir[idir][i]=new unsigned char[columns];
  }
  
  cout<<"  fill_region: reference pixel (x,y)= ("<<
    refPxl[0]<<","<<refPxl[1]<<")  wSize= "<<wSize<<endl;
  cout<<"         foreground= "<<foreground<<"  background= "
      <<background<<"  mode: "<<mode<<"  save= "<<sdum<<endl;
  
  wx=rows/wSize; wy=columns/wSize;
  if (rows%wSize!=0) ++wx; // take care of remaining pixels
  if (columns%wSize!=0) ++wy; // take care of remaining pixels
  grid=new char*[wx];
  for (i=0;i<wx;i++) grid[i]=new char[wy];
  for (i=0;i<(int)rows;i++) pxl_tmp[i]=new unsigned char[columns];
  
  for (idir=0;idir<4;idir++) { // 4 dirs
    edge.clear(); 
    
    // Initialize
    for (i=0;i<wx;i++) for (j=0;j<wy;j++) grid[i][j]=0; 
    for (i=0;i<(int)rows;i++) {
      for (j=0;j<(int)columns;j++) {
	if (background==-1)  pxl_dir[idir][i][j]=pxl_gray[i][j];
	else pxl_dir[idir][i][j]=(unsigned char)background;
      }
    }
    ref_intensity=pxl_gray[refPxl[1]][refPxl[0]]; // flip x/y coord to match image
    refGrid[0]=refPxl[1]/wSize;  refGrid[1]=refPxl[0]/wSize; // note x/y flipped
    edge.push_back(pair<int,int>(refGrid[0],refGrid[1]));
    cflag=0; // cflag of whether filling is done
    
    while (cflag==0) {
      for (it=edge.begin();it!=edge.end(); it++) {
	gx=(*it).first; gy=(*it).second;
	x0=gx*wSize; y0=gy*wSize; 
	
	/* search direction      */
	if (idir==0) {
	  x1=x0+wSize; x1=(x1<(int)rows)?x1:rows;
	  y1=y0+wSize; y1=(y1<(int)columns)?y1:columns;
	  for (i=x0;i<x1;i++) {
	    for (j=y0;j<y1;j++) {
	      if (pxl_gray[i][j]!=ref_intensity) {
		grid[gx][gy]=2; // seal
		goto fR0; 
	      }
	    }
	  }
	}
	else if (idir==1) {
	  x1=x0-wSize; x1=(x1>0)?x1:0;
	  y1=y0-wSize; y1=(y1>0)?y1:0;
	  for (i=x0;i>x1;i--) {
	    for (j=y0;j>y1;j--) {
	      if (pxl_gray[i][j]!=ref_intensity) {
		grid[gx][gy]=2; // seal
		goto fR0; 
	      }
	    }
	  }
	}
	else if (idir==2) {
	  x1=x0+wSize; x1=(x1<(int)rows)?x1:rows;
	  y1=y0-wSize; y1=(y1>0)?y1:0;
	  for (i=x0;i<x1;i++) {
	    for (j=y0;j>y1;j--) {
	      if (pxl_gray[i][j]!=ref_intensity) {
		grid[gx][gy]=2; // seal
		goto fR0; 
	      }
	    }
	  }
	}
	else if (idir==3) {
	  x1=x0-wSize; x1=(x1>0)?x1:0;
	  y1=y0+wSize; y1=(y1<(int)columns)?y1:columns;
	  for (i=x0;i>x1;i--) {
	    for (j=y0;j<y1;j++) {
	      if (pxl_gray[i][j]!=ref_intensity) {
		grid[gx][gy]=2; // seal
		goto fR0;
	      }
	    }
	  }
	}	
      fR0: continue; //if (cflag==0) grid[gx][gy]=1; // checked
      } //  for (it=edge.begin();it!=edge.end(); it++) {
      
      for (it=edge.begin();it!=edge.end(); it++) { // mark old edge
	gx=(*it).first; gy=(*it).second;
	if (grid[gx][gy]!=2) {
	  assert(grid[gx][gy]==0);
	  grid[gx][gy]=1; // mark non-sealed edge
	}
      }
      
      // find new edge around non-sealed edge
      edge1.clear();
      for (it=edge.begin();it!=edge.end(); it++) {
	gx=(*it).first; gy=(*it).second;
	if (grid[gx][gy]==2) continue; // sealed edge
	i=gx; // middle
	if (gy>0) { // up: not the first row
	  j=gy-1; if (grid[i][j]==0) edge1.push_back(pair<int,int>(i,j));
	}
	if (gy<(wy-1)) { // below: not the last row
	  j=gy+1; if (grid[i][j]==0) edge1.push_back(pair<int,int>(i,j));
	}
	if (gx>0) { // left: not the first col
	  i=gx-1; j=gy; if (grid[i][j]==0) edge1.push_back(pair<int,int>(i,j));
	  if (gy>0) { // up: not the first row
	    j=gy-1; if (grid[i][j]==0) edge1.push_back(pair<int,int>(i,j));
	  }
	  if (gy<(wy-1)) { // below: not the last row
	    j=gy+1; if (grid[i][j]==0) edge1.push_back(pair<int,int>(i,j));
	  }
	}
	if (gx<(wx-1)) { // right: not the last col
	  i=gx+1; j=gy; if (grid[i][j]==0) edge1.push_back(pair<int,int>(i,j));
	  if (gy>0) { // up: not the first row
	    j=gy-1; if (grid[i][j]==0) edge1.push_back(pair<int,int>(i,j));
	  }
	  if (gy<(wy-1)) { // below: not the last row
	    j=gy+1; if (grid[i][j]==0) edge1.push_back(pair<int,int>(i,j));
	  }
	}
      } // for (it=edge.begin();it!=edge.end(); it++) {
      sort(edge1.begin(),edge1.end()); // remove duplicate entries
      edge1.erase(unique(edge1.begin(),edge1.end()),edge1.end());
      if ((int)(edge1.size())== 0) {cflag=1; break;}
      else { // copy edge1 to edge
	edge.clear(); edge=edge1;
      }
    } //   while (cflag==0) {
    
    
    // fill pxl_dir
    //if (mode=="region") { // Fill the region with foreground.
    for (gx=0;gx<wx;gx++) {
      for (gy=0;gy<wy;gy++) {
	if (((grid[gx][gy]==1)&&(mode=="region"))||
	    ((grid[gx][gy]==2)&&(mode=="edge"))) {
	  
	  x0=gx*wSize; y0=gy*wSize; 
	  
	  if (idir==0 ) {
	    x1=x0+wSize; x1=(x1<(int)rows)?x1:rows;
	    y1=y0+wSize; y1=(y1<(int)columns)?y1:columns;
	    for (i=x0;i<x1;i++) {
	      for (j=y0;j<y1;j++) {
		if (foreground==-1) pxl_dir[idir][i][j]=pxl_gray[i][j];
		else pxl_dir[idir][i][j]=(unsigned char)foreground;
	      }
	    }
	  }
	  else if (idir==1) {
	    x1=x0-wSize; x1=(x1>0)?x1:0;
	    y1=y0-wSize; y1=(y1>0)?y1:0;
	    for (i=x0;i>x1;i--) {
	      for (j=y0;j>y1;j--) {
		if (foreground==-1) pxl_dir[idir][i][j]=pxl_gray[i][j];
		else pxl_dir[idir][i][j]=(unsigned char)foreground;
	      }
	    }
	  }
	  else if (idir==2) {
	    x1=x0+wSize; x1=(x1<(int)rows)?x1:rows;
	    y1=y0-wSize; y1=(y1>0)?y1:0;
	    for (i=x0;i<x1;i++) {
	      for (j=y0;j>y1;j--) {
		if (foreground==-1) pxl_dir[idir][i][j]=pxl_gray[i][j];
		else pxl_dir[idir][i][j]=(unsigned char)foreground;
	      }
	    }
	  }
	  else if (idir==3) {
	    x1=x0-wSize; x1=(x1>0)?x1:0;
	    y1=y0+wSize; y1=(y1<(int)columns)?y1:columns;
	    for (i=x0;i>x1;i--) {
	      for (j=y0;j<y1;j++) {
		if (foreground==-1) pxl_dir[idir][i][j]=pxl_gray[i][j];
		else pxl_dir[idir][i][j]=(unsigned char)foreground;
	      }
	    }
	  }
	}
      } // for (gy=0;gy<wx;gy++) {      
    } // for (gx=0;gx<wx;gx++) {
  } //   for (idir=0;idir<4;idir++) { // 4 dirs
  
  
  // combine region from all four directions
  for (i=0;i<(int)rows;i++) {
    for (j=0;j<(int)columns;j++)
      pxl_tmp[i][j]=(unsigned char)background;  // clear all
  }
  
  for (idir=0;idir<4;idir++) {
    for (i=0;i<(int)rows;i++) {
      for (j=0;j<(int)columns;j++) {
	if (pxl_tmp[i][j]==background) pxl_tmp[i][j]+=pxl_dir[idir][i][j];
      }
    }
  } //   for (idir=0;idir<4;idir++) {
  
  // to get "mode edge" only print pixels on the edges of pxl_tmp
  // include diagonals
  if (mode=="edge") { // only keep pixels on the region edge
    inner.clear();
    for (i=0;i<(int)rows;i++) {
      flag=0;
      for (j=0;j<(int)columns;j++) {
	flag=0;       	
	if (pxl_tmp[i][j]==background) continue; 
	x1=i; y1=(j+wSize<(int)columns)?(j+wSize):(j);  	// top
	if (pxl_tmp[x1][y1]!=pxl_tmp[i][j])  flag=1; 	
	x1=i; y1=((j-wSize)>=0)?(j-wSize):(0);  	// bottom
	if (pxl_tmp[x1][y1]!=pxl_tmp[i][j]) flag=1; 
	x1=((i-wSize)>=0)?(i-wSize):(0); y1=j;  	// left
	if (pxl_tmp[x1][y1]!=pxl_tmp[i][j]) flag=1; 
	x1=((i+wSize)<(int)rows)?(i+wSize):(i); y1=j;  	// right
	if (pxl_tmp[x1][y1]!=pxl_tmp[i][j]) flag=1;
	if (!flag) inner.push_back(pair<int,int>(i,j));
	else continue;
	
	// diagonal directions:
	// if all four diagonals wsize away from i,j have same intensity --> edge
	flag=1; 
	x1=((i-wSize)>=0)?(i-wSize):(0);  y1=((j-wSize)>=0)?(j-wSize):(0);
	if (pxl_tmp[x1][y1]!=pxl_tmp[i][j]) flag=0; 	
	x1=((i-wSize)>=0)?(i-wSize):(0);  y1=(j+wSize<(int)columns)?(j+1):(j);  
	if (pxl_tmp[x1][y1]!=pxl_tmp[i][j]) flag=0; 		
	x1=((i+wSize)<(int)rows)?(i+wSize):(i); y1=((j-wSize)>=0)?(j-wSize):(0);
	if (pxl_tmp[x1][y1]!=pxl_tmp[i][j]) flag=0; 			
	x1=((i+wSize)<(int)rows)?(i+wSize):(i); y1=(j+wSize<(int)columns)?(j+wSize):(j);  
	if (pxl_tmp[x1][y1]!=pxl_tmp[i][j]) flag=0;
	if (!flag) {inner.pop_back(); }
      }
    }
    sort(inner.begin(),inner.end()); // remove duplicate entries
    inner.erase(unique(inner.begin(),inner.end()),inner.end());
    for (it=inner.begin();it!=inner.end(); it++)  {
      x0=(*it).first; y0=(*it).second;
      if (background==-1) pxl_tmp[x0][y0]=pxl_gray[x0][y0];
      else pxl_tmp[x0][y0]=(unsigned char)background;
    }
  }    
  
  if (tag1=="") {
    for (i=0;i<(int)rows;i++) {
      for (j=0;j<(int)columns;j++) pxl_gray[i][j]=pxl_tmp[i][j];
    }
    cout<<"    Img "<<tag<<" modified."<<endl;
  }
  else { // save to other existing tag
    it_img=findImg(tag1);
    if (it_img==vimg.end()) error_dodri(" image tag not found");
    for (i=0;i<(int)rows;i++) {
      for (j=0;j<(int)columns;j++) {
	it_img->pxl_gray[i][j]=pxl_tmp[i][j];
      }
    }
    cout<<"    Img "<<tag<<" modified."<<endl;
  }
    
  for (i=0;i<(int)rows;i++) delete[] pxl_tmp[i]; 
  for (i=0;i<wx;i++) delete [] grid[i];  
  delete [] pxl_tmp; delete [] grid;  delete [] refGrid;
  return;
}

/************************************************************************/
int img::findCluster(unsigned char pxlCut,
		     vector<vector<pair<int,int>>> &cluster)
{ /* Find clusters of pixels with intensity > pxlCut
     If two clusters are connected only diagonally (no horiz/vert connection),
     treat them separate.     
     Returns: number of clusters found
  */

  int **pxl_clstr=new int*[rows]; // stores cluster number for each pixel
  int i,j,idum,idum1,idum2,ncluster=0;

  vector<pair<int,int>> subCluster;
  vector<pair<int,int>>::iterator it;

  cout<<"         findCluster: " <<
    "  pxlCut= "<<(int)pxlCut<<endl; 
  
  // Initialize
  for (i=0;i<(int)rows;i++) {
    pxl_clstr[i]=new int[columns];
    for (j=0;j<(int)columns;j++) pxl_clstr[i][j]=0;
  }
  
  cluster.clear(); 

  if (pxl_gray[0][0]>pxlCut) {   // first row, first col
    assert(pxl_clstr[0][0]==0);
    ++ncluster;  pxl_clstr[0][0]=ncluster;
    subCluster.clear(); subCluster.push_back(pair<int,int>(0,0));
    cluster.push_back(subCluster);
  }
  
  for (j=1;j<(int)columns;j++) { // first row
    if (pxl_gray[0][j]>pxlCut) {
      assert (pxl_clstr[0][j]==0);
      idum=pxl_clstr[0][j-1];
      if (idum>0) { // add to existing cluster
	pxl_clstr[0][j]=idum;
	cluster[idum-1].push_back(pair<int,int>(0,j));
      }
      else { // new cluster
	++ncluster;  pxl_clstr[0][j]=ncluster;
	subCluster.clear(); subCluster.push_back(pair<int,int>(0,j));
	cluster.push_back(subCluster);
      } // else { // new cluster
    } // if (pxl_gray[0][j]>pxlCut) {
  } // for (j=1;j<columns;j++) { // first row
  
  for (i=1;i<(int)rows;i++) { // first column
    if (pxl_gray[i][0]>pxlCut) {
      assert (pxl_clstr[i][0]==0);
      idum=pxl_clstr[i-1][0];
      if (idum>0) { // add to existing cluster
	pxl_clstr[i][0]=idum;
	cluster[idum-1].push_back(pair<int,int>(i,0));
      }
      else { // new cluster
	++ncluster;  pxl_clstr[i][0]=ncluster;
	subCluster.clear(); subCluster.push_back(pair<int,int>(i,0));
	cluster.push_back(subCluster);
      } // else { // new cluster
    } // if (pxl_gray[i][0]>pxlCut) {
  } // for (i=1;i<rows;i++) { // first column

  for (i=1;i<(int)rows;i++) {
    for (j=1;j<(int)columns;j++) {
      if (pxl_gray[i][j]>pxlCut) {
	assert (pxl_clstr[i][j]==0);
	idum=pxl_clstr[i-1][j]; idum1=pxl_clstr[i][j-1];
	if (idum==0 && idum1==0) { // new cluster
	  ++ncluster;  pxl_clstr[i][j]=ncluster;
	  subCluster.clear(); subCluster.push_back(pair<int,int>(i,j));
	  cluster.push_back(subCluster);
	}
	else if (idum!=0 && idum1==0) { // add to idum
	  pxl_clstr[i][j]=idum;
	  cluster[idum-1].push_back(pair<int,int>(i,j));
	}
	else if (idum==0 && idum1!=0) {  // add to idum1
	  pxl_clstr[i][j]=idum1;
	  cluster[idum1-1].push_back(pair<int,int>(i,j));
	}
	else { // two clusters meet.
	  if (idum==idum1) {  // two clusters are the same
	    pxl_clstr[i][j]=idum;
	  }
	  else { // merge to smaller cluster number
	    idum2=(idum<idum1)?idum:idum1; // smaller
	    idum1=(idum<idum1)?idum1:idum; idum=idum2;
	    for (it=cluster[idum1-1].begin();it!=cluster[idum1-1].end();it++) {
	      cluster[idum-1].push_back(pair<int,int>((*it).first,(*it).second));
	      pxl_clstr[(*it).first][(*it).second]=idum;
	    }
	    idum2=ncluster-1;
	    if (idum1<ncluster) {
	      cluster[idum1-1].clear(); // move ncluster to idum1
	      for (it=cluster[idum2].begin();it!=cluster[idum2].end();it++) {
		cluster[idum1-1].push_back(pair<int,int>
					   ((*it).first,(*it).second));
		pxl_clstr[(*it).first][(*it).second]=idum1;
	      }
	    }
	    cluster.pop_back(); // remove the last element of cluster[]
	    --ncluster;
	  } // else { // merge to smaller cluster number
	  cluster[idum-1].push_back(pair<int,int>(i,j));
	} // else { // two clusters meet.
      } // if (pxl_gray[i][j]>pxlCut) {
    } //  for (j=0;j<columns;j++) {
  } //  for (i=0;i<rows;i++) {

  for (i=0;i<(int)rows;i++) delete[] pxl_clstr[i];
  delete [] pxl_clstr;

  cout<<"         Number of clusters found: "<<ncluster<<endl;
  
  return ncluster;
}

/************************************************************************/
void img::gaussBlur(double *s,int *w, unsigned char pxlCut, string tag1)
{ /* 2D gaussian blur filter for thresholding. Used to smoothen image.
     Only include pixels with intensity > pxlCut. This is to avoid low
     intensity pixels to dim images further, which can be problematic when
     features are dim and contain pixels w/ intensity comparable to the
     background.
  */
  
  int i,j,cnt;
  double sx2=2.0*s[0]*s[0], sy2=2.0*s[1]*s[1]; // 2*var of gaussian kernel
  int wx=w[0], wy=w[1]; //kernel width, can be asymetric, must be int
  double kernel[wx][wy]; //setup kernel
  int kx=(wx-1)/2, ky=(wy-1)/2; // center of window
  int ix,iy,ix0,ix1,iy0,iy1; // kernel related indices
  double sum=0., rdum;
  int irow,icol;
  unsigned char **pxl_blur=new unsigned char*[rows], pxl0;
  string sdum=(tag1=="")?tag:tag1;
  
  cout<<"  Gaussian:  wx= "<<wx<<"  wy= "<<wy<<"  sx= "<<s[0]
      <<"  sy= "<<s[1]<<" pxlCut= "<<(int)pxlCut<<" save= "<<sdum<<endl;
  
  for (irow=0;irow<(int)rows;irow++) pxl_blur[irow]=new unsigned char[columns];
  for (irow=0;irow<(int)rows;irow++) {
    for (icol=0;icol<(int)columns;icol++) pxl_blur[irow][icol]=0;
  }
  
  // Create Gaussian kernel
  for (i=0;i<wx;i++) {
    for (j=0;j<wy;j++) {
      // need to use kx, ky as ref
      rdum= pow((double)(i-kx),2)/sx2 + pow((double)(j-ky),2)/sy2;
      kernel[i][j] = exp(-1.*rdum ); 
      sum+=kernel[i][j];
    }
  }  
  for (i=0;i<wx;i++) { //normalize
    for (j=0;j<wy;j++)  kernel[i][j]/=sum;
  }

  // Apply the kernel
  for (irow=0;irow<(int)rows;irow++) {
    for (icol=0;icol<(int)columns;icol++) {
      // set the bounaries of window to apply kernel
      ix0=(irow>=kx)?0:(kx-irow);
      iy0=(icol>=ky)?0:(ky-icol);
      ix1=(((int)rows-irow)>kx)?wx:(kx+rows-irow);
      iy1=(((int)columns-icol)>ky)?wy:(ky+columns-icol);
      
      sum=rdum=0.; cnt=0;
      for (ix=ix0;ix<ix1;ix++) { 
	for (iy=iy0;iy<iy1;iy++) {
	  pxl0=pxl_gray[irow-kx+ix][icol-ky+iy];
	  if (pxl0>pxlCut) {
	    rdum+=kernel[ix][iy]; ++cnt; // count included pixels
	    sum+=(double)pxl0*kernel[ix][iy];
	  }
	}
      }
      if (cnt==0) pxl_blur[irow][icol]=pxl_gray[irow][icol]; // do nothing
      else pxl_blur[irow][icol]=(unsigned char)(sum/rdum);
    } // for (icol=0;icol<columns;icol++) {
  } // for (irow=0;irow<rows;irow++) { //skipping edges for now
  
  if (tag1=="") { //save blur
    for (irow = 0; irow < (int)rows; irow++) {
      for (icol = 0; icol < (int)columns; icol++)
	pxl_gray[irow][icol]=pxl_blur[irow][icol];
    }
    cout<<"    Img "<<tag<<" modified." <<endl;
  }
  else copyImg(*this,pxl_blur,tag1);

  for (irow=0;irow<(int)rows;irow++) delete[] pxl_blur[irow];
  delete[] pxl_blur;
  return;
}

/************************************************************************/
void img::getPxlAreaFraction(unsigned char *pvec, double locut,
			     double hicut)
{ /* Find pixel intensity ranges pvec[2] where the fraction of the number
     of pixels is between locut and hicut. */
  
  if ((locut<0.)||(locut>1.)||(hicut<0.)||(hicut>1.))
    error_dodri(" Fractional cutoff must be between 0 and 1.");
  int i,j,*pxl_cnt=new int[256], idum;
  int const npixel=(int)rows*(int)columns;
  double *pxl_frac=new double[256], rdum;
  unsigned char cdum;
  
  for (i=0;i<256;i++) {pxl_cnt[i]=0; pxl_frac[i]=0.;}
  for (i=0;i<(int)rows;i++) {
    for (j=0;j<(int)columns;j++) {
      cdum=pxl_gray[i][j];
      assert((cdum>=0)&&(cdum<256));
      ++pxl_cnt[(int)cdum];
    }
  }
  idum=0; rdum=1./(double)npixel;
  for (i=0;i<256;i++) {
    idum+=pxl_cnt[i];
    pxl_frac[i]=rdum*(double)idum;
  }
  assert(idum==npixel);
  pvec[0]=0; pvec[1]=255; // default
  for (i=0;i<255;i++) {
    if ((pxl_frac[i]<locut)&&(pxl_frac[i+1]>=locut))
      pvec[0]=(unsigned char)i;
    if ((pxl_frac[i]<hicut)&&(pxl_frac[i+1]>=hicut))
      pvec[1]=(unsigned char)(i+1);
  }  

  delete [] pxl_frac; delete[] pxl_cnt;
  return;
}

/************************************************************************/
void img::getStatistics(int loCut, int hiCut, double *avg, double *sd)
{ /* Get avg and std of intensity for pxl_gray with intensity values in the
     range loCut and hiCut */
  
  int pxl_cnt=0, irow,icol,idum, pmin,pmax;
  double avg0=0.,sd0=0., rdum;
  
  pmin=hiCut; pmax=0;
  for (irow = 0; irow < (int)rows; irow++) { 
    for (icol=0;icol<(int)columns;icol++) {
      idum=(int)pxl_gray[irow][icol];
      if ((idum>=loCut)&&(idum<hiCut)) {
	if ((idum<pmin)&&(idum>0)) pmin=idum;
	if (idum>pmax) pmax=idum;
	rdum=(double)idum; avg0+=rdum; sd0+=rdum*rdum;
	++pxl_cnt;
      }
    } // for (icol=0;icol<(int)columns;icol++) {
  } // for (irow = 0; irow < (int)rows; irow++) {
  
  cout << "  Number of pixels in intensity range " <<loCut<<" - "<<hiCut-1
       << " : "<<pxl_cnt<< endl;
  cout << "  nonzero minimum and max pixel intensities in the range: "
       << pmin << "  "<<pmax<<endl;
  if (pxl_cnt==0) {  *avg=0.; *sd=0.; }
  else {
    avg0/=(double)pxl_cnt;  sd0/=(double)pxl_cnt; sd0-=(avg0*avg0); 
    assert(sd0>=0.); sd0=sqrt(sd0);
    *avg=avg0; *sd=sd0;
  }
  cout<<"  Avg/std of grayscale pixel intensity: "<<setw(8)<<setprecision(3)
      << *avg <<"  "<< *sd <<endl;
  return;
}

/************************************************************************/
void img::invert(string tag1)
{ /* Invert pixel intensity of image. Subtract all pxl_gray values from 
     max intensity.   */ 
  
  int irow,icol,idum,pmax,umax;
  unsigned char **pxl_tmp=new unsigned char*[rows];
  for (irow=0;irow<(int)rows;irow++) pxl_tmp[irow]=new unsigned char[columns];    
  string sdum=(tag1=="")?tag:tag1;
  
  cout<<"  invert: save= "<<sdum<<endl;
  
  pmax=0;
  for (irow = 0; irow < (int)rows; irow++) { 
    for (icol=0;icol<(int)columns;icol++) {
      idum=pxl_gray[irow][icol];
      pmax=(idum>pmax)?idum:pmax; 
    } // for (icol=0;icol<(int)columns;icol++) {
  } // for (irow = 0; irow < (int)rows; irow++) {


  umax=(unsigned char)pmax; 
  for (irow = 0; irow < (int)rows; irow++) { 
    for (icol=0;icol<(int)columns;icol++) 
      pxl_tmp[irow][icol]=umax-pxl_gray[irow][icol];
  } 
  
  if (tag1=="") {
    for (irow = 0; irow < (int)rows; irow++) { 
      for (icol=0;icol<(int)columns;icol++)
	pxl_gray[irow][icol]=pxl_tmp[irow][icol];
    }
    cout<<"    Image "<<tag<<" modified."<<endl;
  }
  else copyImg(*this,pxl_tmp,tag1);

  //clear memory
  for (irow=0;irow<(int)rows;irow++) delete[] pxl_tmp[irow];
  delete [] pxl_tmp;

  return; 
}

/************************************************************************/
void img::removeDebris(int sizeCut, unsigned char *pxlCut,
		       unsigned char fillPxl,string mode,string tag1)
{ /* Find isolated pixel clusters with pixel intensity greater than
     pxlCut[0]. Then delete pixels of the cluster under the following
     conditions:
     
     pxlCut[1]=0: Only perform removal based on cluster size.
     pxlCut[1]>0: Remove clusters that satisfy both the size cutoff 
                  and max intensity of the cluster less than pxlCut[1].
		  This allows removing small dim areas while dim pixels 
		  in a cluster containing bright ones are kept.

		  If sizeCut=0, remove all pixels where max intensity is
		  less than pxlCut[1]. In this case, mode is ignored.

     fillPxl: pixel intensity for filling deleted areas.
     mode=area: delete clusters whose area (pixel count) is less than sizeCut
          rank: delete clusters whose size is less than sizeCut-th among
	        identified clusters (largest cluster is rank 1).
     tag1: If non-empty, save output to new img. 
  */
  
  unsigned char **pxl_tmp=new unsigned char*[rows], pmax;
  int i,j,idum,idum1,idum2,ncluster=0,pxl_cnt;
  vector<vector<pair<int,int>>> cluster;
  vector<pair<int,int>>::iterator it;
  multimap<int,int> cluster_size;
  multimap<int,int>::iterator it1;
  multimap<int,int>::reverse_iterator rit1;
  string sdum=(tag1=="")?tag:tag1;

  cout<<"  remove_debris: mode= "<<mode<<"  size_cut= "<<sizeCut
      <<"  pxl_bg= "<< (int)pxlCut[0] <<endl;
  cout<<"         pxl_cut= "<< (int)pxlCut[1]
      <<"  fill= "<<(int)fillPxl <<"  save= "<<sdum<<endl;
  
  // Initialize
  for (i=0;i<(int)rows;i++) pxl_tmp[i]=new unsigned char[columns];
  for (i=0;i<(int)rows;i++) {
    for (j=0;j<(int)columns;j++)  pxl_tmp[i][j]=pxl_gray[i][j]; 
  }
  
  ncluster=findCluster(pxlCut[0], cluster);
  if (ncluster==0) {
    cout<<"  removeDebris: No cluster found."<<endl;
    return;
  }

  // rank cluster sizes
  for (i=0;i<ncluster;i++) {
    idum=cluster[i].size(); // sort according to cluster size
    cluster_size.insert(pair<int,int>(idum,i+1)); 
  }
  pxl_cnt=0;
  
  if (sizeCut>0) {
    if (mode=="area") { // remove small clusters
      for (it1=cluster_size.begin();it1!=cluster_size.end();it1++) {
	if (it1->first < sizeCut) {
	  idum1=it1->second -1 ; // cluster number
	  pmax=0;
	  if (pxlCut[1]>0) { // find max pixel intensity
	    for (it=cluster[idum1].begin();it!=cluster[idum1].end();it++) {
	      i=it->first; j=it->second;
	      if (pmax<pxl_tmp[i][j]) pmax=pxl_tmp[i][j];
	    }
	  }
	  if ((pxlCut[1]==0)||(pmax<pxlCut[1])) {
	    for (it=cluster[idum1].begin();it!=cluster[idum1].end();it++) {
	      i=(*it).first; j=(*it).second; pxl_tmp[i][j]=fillPxl;
	      ++pxl_cnt;
	    }
	  }
	} // if (it1->first < sizeCut) {
      } // for (it1=cluster_size.begin();it1!=cluster_size.end();it1++) {
    }  //  if (mode=="area") { // remove small clusters
    else if (mode=="rank") { // keep large clusters
      idum2=0; // cluster count
      for (rit1=cluster_size.rbegin();rit1!=cluster_size.rend();rit1++) {
	++idum2;
	if (idum2>sizeCut) {
	  idum1=rit1->second - 1; // cluster number
	  pmax=0;
	  if (pxlCut[1]>0) { // find max pixel intensity
	    for (it=cluster[idum1].begin();it!=cluster[idum1].end();it++) {
	      i=it->first; j=it->second;
	      if (pmax<pxl_tmp[i][j]) pmax=pxl_tmp[i][j];
	    }
	  }
	  if ((pxlCut[1]==0)||(pmax<pxlCut[1])) {
	    for (it=cluster[idum1].begin();it!=cluster[idum1].end();it++) {
	      i=(*it).first; j=(*it).second; pxl_tmp[i][j]=fillPxl;
	      ++pxl_cnt;
	    }
	  }
	} // if (idum2>sizeCut) {
      } //for (rit1=cluster_size.rbegin();rit1!=cluster_size.rend();rit1++) {
    } //  else if (mod=="rank") { // keep large clusters
    else error_dodri("Unrecognized mode: "+mode);
  } // if (sizeCut>0) {
  else { // use only pxlCut[1]
    if (pxlCut[1]<1) 
      cout <<"  WARNING: Both size_cut and pxl_cut zero. Nothing done."<<endl;
    else {
      for (it1=cluster_size.begin();it1!=cluster_size.end();it1++) {
	idum1=it1->second -1 ; // cluster number
	pmax=0;
	if (pxlCut[1]>0) { // find max pixel intensity
	  for (it=cluster[idum1].begin();it!=cluster[idum1].end();it++) {
	    i=it->first; j=it->second;
	    if (pmax<pxl_tmp[i][j]) pmax=pxl_tmp[i][j];
	  }
	}
	if ((pxlCut[1]==0)||(pmax<pxlCut[1])) {
	  for (it=cluster[idum1].begin();it!=cluster[idum1].end();it++) {
	    i=(*it).first; j=(*it).second; pxl_tmp[i][j]=fillPxl;
	    ++pxl_cnt;
	  }
	}
      } // for (it1=cluster_size.begin();it1!=cluster_size.end();it1++) {
    }
  } //  else { // use only pxlCut[1]
  
  cout<<"         Number of pixels removed: "<<pxl_cnt<<endl;
  if (tag1=="") {
    for (i=0;i<(int)rows;i++) {
      for (j=0;j<(int)columns;j++) {
	pxl_gray[i][j]=pxl_tmp[i][j];
      }
    }
    cout<<"    Image "<<tag<<" modified."<<endl;
  }
  else copyImg(*this,pxl_tmp,tag1);

  for (i=0;i<(int)rows;i++) delete[] pxl_tmp[i]; 
  delete [] pxl_tmp; 

  return;
}

/************************************************************************/
void img::shift(int *dr, unsigned char marginPxl, string tag1,string tag_ref) 
{/* Shift img by (dr[0],dr[1]).
    marginPxl: Intensity of pixels left out by shift
    
    Output written to new img tag1 and is added to vimg. 
    If tag1 is empty, the original pxl_gray is overwritten.  

    If tag_ref!="", when shifting, instead of using marginPxl to fill, 
    use respective pixels from tag_ref. 
 */
  int i,irow, icol,a,b;
  unsigned char **pxl_shift=new unsigned char*[rows] ;
  string sdum=(tag1=="")?tag:tag1;
  vector<img>::iterator it_img; // for tag_ref 

  cout<<"  shift: dx= "<<dr[0]<<"  dy= "<<dr[1]<<"  fill= "
      <<(int)marginPxl<<"  save= "<<sdum<<endl;

  if (tag_ref!="") cout<<"    Use image="<<tag_ref<<" to fill shift."<<endl; 
  
  for (irow=0;irow<(int)rows;irow++) pxl_shift[irow]=new unsigned char[columns];

  if (tag_ref=="") { // no reference tag for filling, use marginPxl
    for (irow=0;irow<(int)rows;irow++) {
      for (icol=0;icol<(int)columns;icol++) pxl_shift[irow][icol]=marginPxl;
    }
  }
  else { // use intensity of pixels from reference tag for filling
    it_img=findImg(tag_ref);
    if (it_img==vimg.end())  error_dodri("Image "+tag_ref+" not found."); 
    for (irow=0;irow<(int)rows;irow++) {
      for (icol=0;icol<(int)columns;icol++)
	pxl_shift[irow][icol]=it_img->pxl_gray[irow][icol];
    }
  }
    
  for (irow=0;irow<(int)rows;irow++){
    for (icol=0;icol<(int)columns;icol++) {
      a = irow-dr[1];  b = icol+dr[0];
      if ((a<0) || (a>((int)rows-1))) continue; //checking boundaries of image
      if ((b<0) || (b>((int)columns-1))) continue;
      pxl_shift[a][b] = pxl_gray[irow][icol];
    }
  }

  if (tag1=="") {
    for (irow = 0; irow < (int)rows; irow++) {
      for (icol = 0; icol < (int)columns; icol++) {
	pxl_gray[irow][icol]=pxl_shift[irow][icol];
      }
    }
    cout<<"    Image "<<tag<<" modified."<<endl;
  }
  else copyImg(*this,pxl_shift,tag1);
  
  for (i=0;i<(int)rows;i++) delete[] pxl_shift[i];
  delete[] pxl_shift;

  return;
}

/************************************************************************/
void img::threshold(string kind,string mode,double *cutoff,int flag, string tag1)  
{ /* Perform thresholding (filtering)
     kind=highpass/lowpass/bandpass
     mode=absolute (pixel intensity value in the 0--255 range)
          relative (avg+cutoff*sd)
	  fraction: keep given percent number of pixels
	  auto: works only for bipass as of 170818
	  cutoff: two values for bandpass, single value otherwise.
     If the flag for binary mode is true, then the thresholding is binary.

     Output written to new img tag1 and is added to vimg. If tag1 is empty,
     the original pxl_gray is overwritten.  
  */
  
  int i,irow, icol;
  double rdum,avg,sd;
  unsigned char pvec[2];
  unsigned char **pxl_thresh = new unsigned char*[rows];
  string sdum=(tag1=="")?tag:tag1;
  
  cout<<"  Threshold: kind= "<<kind<<"  mode= "<<mode
      <<"  cutoff= "<<cutoff[0]<<endl; 
  if (mode=="bandpass") cout<<"          "<<cutoff[1];
  cout <<"         binary_flag= "<<flag<<"  save= "<<sdum  << endl;
  
  if (mode=="relative") getStatistics(10,range,&avg,&sd);
  
  if (kind!="bandpass") { // not bandpass
    if (mode=="relative") {
      rdum=avg+sd*cutoff[0];
      pvec[0]=(rdum<0.)?0:(unsigned char)rdum; 
      pvec[0]=(rdum>255)?255:(unsigned char)rdum;
    }
    else if (mode=="absolute") {
      pvec[0]=(cutoff[0]<0)?0:(unsigned char)cutoff[0];
      pvec[0]=(cutoff[0]>255)?255:(unsigned char)cutoff[0];
    }
    else { // fraction
      getPxlAreaFraction(pvec,cutoff[0],1.);
    }
    pvec[1]=255; // highpass by default
    if (kind=="lowpass") {
      pvec[1]=pvec[0]; pvec[0]=0;
    }
  } // if (kind!="bandpass") { // single cutoff
  else { // bandpass: two cutoffs
    if (mode=="relative") {
      rdum=avg+sd*cutoff[0];
      pvec[0]=(rdum<0.)?0:(unsigned char)rdum; 
      pvec[0]=(rdum>255)?255:(unsigned char)rdum;
      rdum=avg+sd*cutoff[1];
      pvec[1]=(rdum<0.)?0:(unsigned char)rdum; 
      pvec[1]=(rdum>255)?255:(unsigned char)rdum;
    }
    else if (mode=="absolute") {
      pvec[0]=((unsigned char)cutoff[0]<0)?0:(unsigned char)cutoff[0];
      pvec[0]=((unsigned char)cutoff[0]>255)?255:(unsigned char)cutoff[0];
      pvec[1]=((unsigned char)cutoff[1]<0)?0:(unsigned char)cutoff[1];
      pvec[1]=((unsigned char)cutoff[1]>255)?255:(unsigned char)cutoff[1];
    }
    else { // fraction
      getPxlAreaFraction(pvec,cutoff[0],cutoff[1]);
    }
  }
  
  if (pvec[0]>=255.) {
    cout<<"  WARNING: locut greater than 255. Nothing to do"<<endl;
    return;
  }
  if (pvec[1]<1.) {
    cout<<"  WARNING: hicut less than 1. Nothing to do"<<endl;
    return;
  }
  
  for (irow = 0; irow < (int)rows; irow++) 
    pxl_thresh[irow] = new unsigned char[columns];

  for (irow = 0; irow < (int)rows; irow++) {
    for (icol = 0; icol < (int)columns; icol++) {
      if (pxl_gray[irow][icol] >= pvec[0] && pxl_gray[irow][icol] <= pvec[1]) {
	if (flag) pxl_thresh[irow][icol] = 255; // binary
	else pxl_thresh[irow][icol] = pxl_gray[irow][icol];
      }
      else pxl_thresh[irow][icol] = 0;
    }
  }
  cout<<"         Thresholding with low/high pixel intensity: "
      << (int)pvec[0] << " " << (int)pvec[1] <<endl;
  
  if (tag1=="") {
    for (irow = 0; irow < (int)rows; irow++) {
      for (icol = 0; icol < (int)columns; icol++) 
	pxl_gray[irow][icol]=pxl_thresh[irow][icol];
    }
    cout<<"    Image "<<tag<<" modified."<<endl;
  }
  else copyImg(*this,pxl_thresh,tag1);
  
  for (i=0;i<(int)rows;i++) delete[] pxl_thresh[i];
  delete[] pxl_thresh;
  
  return;
}

/************************************************************************/
void img::thresholdAuto(double cutoff, int nbin, string tag1)
{ /* Builds histogram of pixel intensities, then at low end, cut at the
   first minimum of the histogram after the low-intensity background peak
   (call the intensity as pxl0).  The remaining intensity interval is
   [pxl0,255]. set high-intensity cutoff as pxl1=[255-pxl0]*cutoff+pxl0. Then
   squeeze [pxl0,255] into [pxl0,pxl1].
   nbin: number of bins for the histogram

   NOTE: This function assumes that the peak for the background is the highest
   in the image.
  */
  
  int npixel=rows*columns;
  double *pxl_1d=new double[npixel], dBin, rdum0,rdum1,rdum2;
  double **pxl_histo, pxl0,pxl1,avg,sd;
  unsigned char **pxl_tmp=new unsigned char*[rows];
  int i,j,i0=0,i1=0,i2;
  vector<int> local_min, local_max; // local histo min/max
  vector<int>::iterator it;
  
  string sdum=(tag1=="")?tag:tag1;

  if (cutoff>1.) cutoff=1.;
  for (i=0;i<(int)rows;i++) pxl_tmp[i]=new unsigned char[columns];
  if (npixel<(2*nbin)) error_dodri("Image too small for auto threshold.");
  pxl_histo=new double*[2];
  for (i=0;i<2;i++) pxl_histo[i]=new double[nbin];
  rdum0=1.; rdum1=0.;

  for (i=0;i<(int)rows;i++) {
    for (j=0;j<(int)columns;j++) pxl_1d[i*rows+j]=(double)pxl_gray[i][j];
  }
  histogram(pxl_histo,pxl_1d,npixel,nbin,&dBin,&rdum0,&rdum1,0);

  // set highest pixel intensity to work with as avg+sd.
  getStatistics(10,range,&avg,&sd);
  i2=nbin-1;
  for (i=0;i<(nbin-1);i++) {
    if (pxl_histo[0][i]>(avg+sd)) {i2=i; break;}
  }

  // find local min/max
  for (i=1;i<i2;i++) {
    rdum0=pxl_histo[1][i-1]; rdum1=pxl_histo[1][i]; rdum2=pxl_histo[1][i+1];
    
    if ((rdum1<rdum0)&&(rdum1<rdum2)) local_min.push_back(i);
    if ((rdum1>rdum0)&&(rdum1>rdum2)) local_max.push_back(i);
  }
  
  rdum0=1.0e5; i0=0; // find minimum point
  for (i=0;i<(int)local_min.size();++i) {
    j=local_min[i];
    rdum2=pxl_histo[1][j];
    if ((rdum2>0.)&&(rdum2<rdum0)) { rdum0=rdum2; i0=j; }
  }
     
  i1=(int)(cutoff*(double)(nbin-i0))+i0;
  if (i1>=nbin) i1=nbin-1;
  // i0, i1 are indices. Get corresponding pixel intensities
  pxl0=pxl_histo[0][i0]+dBin*0.5; // choose half point of each bin
  pxl1=pxl_histo[0][i1]+dBin*0.5; // choose half point of each bin
  // stretch intensity range [pxl0,pxl1] to [0,255]
  rdum2=pxl1-pxl0; assert(rdum2>0.);
  cout<<"  auto: Pixel intensity range to keep: "<<setw(3)<<(int)pxl0
      <<"  "<<setw(3)<<(int)pxl1<<"  cutoff= "<<setw(3)<<cutoff
      <<"  nbin= "<<nbin<<"  save= "<<sdum<<endl;
  
  for (i=0;i<(int)rows;i++) {
    for (j=0;j<(int)columns;j++) {
      rdum0=(double)pxl_gray[i][j];
      if (rdum0<pxl0) pxl_tmp[i][j]=0;
      else if (rdum0>=pxl0 && rdum0<pxl1) {
	rdum1=(rdum0-pxl0)*255./rdum2;
	if (rdum1>255.) pxl_tmp[i][j]=255;
	else if (rdum1<0.) pxl_tmp[i][j]=0;
	else pxl_tmp[i][j]=(unsigned char)rdum1;
      }
      else pxl_tmp[i][j]=255;
    }
  }

  if (tag1=="") {
    for (i = 0; i < (int)rows; i++) {
      for (j = 0; j < (int)columns; j++) pxl_gray[i][j]=pxl_tmp[i][j];
    }
    cout<<"    Image "<<tag<<" modified."<<endl;
  }
  else copyImg(*this,pxl_tmp,tag1);
  
  for (i=0;i<(int)rows;i++) delete [] pxl_tmp[i];
  for (i=0;i<2;i++) delete [] pxl_histo[i];
  delete [] pxl_tmp; delete [] pxl_histo;   delete [] pxl_1d;
  
  return;
}

/************************************************************************/
void img::thresholdHisto(double *cutoff, int nbin, string tag1)
{ /* Apply histogram-based cutoffs and stretch the pixel intensity.
     
   cutoff[0]<0: cut at the first minimum of the histogram after the
   low-intensity background peak.  This
   assumes that a minimum exists after the peak in histogram.
   cutoff[0]>0: Cut at cutoff[0]*[peak intensity]
   -> Call the intensity as pxl0.

   The remaining intensity interval is [pxl0,255] (pxl0<255). set
   high-intensity cutoff pxl1 as:
   0<cutoff[1]<1.: pxl1=[255-pxl0]*cutoff+pxl0. 
   cutoff[1]<0.||cutoff[1]>1.: pxl1=[highest pixel intensity)

   Then expand [pxl0,pxl1] into [pxl0,255]. If a pixel intensity > pxl1, it
   is set to 255.

   cutoff[0]<0, 0<cutoff[1]<1 corresponds to the old thresholdAuto() ftn
    but more properly handles pixel intensity cutoffs.

   nbin: number of bins for the histogram    */
  
  int npixel=rows*columns;
  double *pxl_1d=new double[npixel], dBin, rdum0,rdum1,rdum2;
  double **pxl_histo, pxl0,pxl1,avg,sd;
  unsigned char **pxl_tmp=new unsigned char*[rows];
  int i,j,i0,i2,imax;
  vector<int> local_min; 
  vector<int>::iterator it;
  
  string sdum=(tag1=="")?tag:tag1;

  for (i=0;i<(int)rows;i++) pxl_tmp[i]=new unsigned char[columns];
  if (npixel<(2*nbin)) error_dodri("Image too small for auto threshold.");
  pxl_histo=new double*[2];
  for (i=0;i<2;i++) pxl_histo[i]=new double[nbin];
  rdum0=1.; rdum1=0.;

  for (i=0;i<(int)rows;i++) {
    for (j=0;j<(int)columns;j++) pxl_1d[i*rows+j]=(double)pxl_gray[i][j];
  }
  histogram(pxl_histo,pxl_1d,npixel,nbin,&dBin,&rdum0,&rdum1,1);
  
  getStatistics(10,range,&avg,&sd);
  // set highest pixel intensity to work with as avg+sd.
  i2=nbin-1;
  for (i=0;i<(nbin-1);i++) {
    if (pxl_histo[0][i]>(avg+sd)) {i2=i; break;}
  }
  // find local min/max
  for (i=1;i<i2;i++) {
    rdum0=pxl_histo[1][i-1]; rdum1=pxl_histo[1][i]; rdum2=pxl_histo[1][i+1];
    if ((rdum1<rdum0)&&(rdum1<rdum2)) local_min.push_back(i);
  }
  
  // set pxl0
  if (cutoff[0]<0.) {
    rdum0=1.0e5; i0=0; // find minimum point
    for (i=0;i<(int)local_min.size();++i) {
      j=local_min[i];
      rdum2=pxl_histo[1][j];
      if ((rdum2>0.)&&(rdum2<rdum0)) { rdum0=rdum2; i0=j; }
    }
    pxl0=pxl_histo[0][i0]+dBin*0.5; // choose half point of each bin
  }
  else {
    rdum0=0.;imax=0;
    for (i=0;i<nbin;i++) {  // get peak of pxl_histo
      if (rdum0<pxl_histo[1][i]) {imax=i; rdum0=pxl_histo[1][i];}
    }
    pxl0=cutoff[0]*pxl_histo[0][imax]+dBin*0.5;
  }
  
  // set pxl1
  if ((cutoff[1]>0.)&&(cutoff[1]<1.)) {
    rdum0=pxl0+dBin*0.5;
    pxl1=cutoff[1]*(pxl_histo[0][nbin-1]-rdum0)+rdum0;
  }
  else { // pxl_histo[*][nbin-1] stores the highest pixel intensity 
    pxl1=pxl_histo[0][nbin-1]+dBin*0.5;
  }
  cout<<"  Using pxl0:"<<pxl0<<" pxl1:"<<pxl1<<endl; 
  if (pxl0>=pxl1) error_dodri("Problem with locut or hicut.");    
    
  // stretch intensity range [pxl0,pxl1] to [0,255]
  rdum2=pxl1-pxl0; assert(rdum2>0.);
  cout<<"  Pixel intensity range to keep: "<<setw(3)<<(int)pxl0
      <<"  "<<setw(3)<<(int)pxl1<<endl;
  cout<<"  locut/hicut= "<<setw(4)
      <<cutoff[0]<<" "<<cutoff[1]
      <<"  nbin= "<<nbin<<"  save= "<<sdum<<endl;
  
  for (i=0;i<(int)rows;i++) {
    for (j=0;j<(int)columns;j++) {
      rdum0=(double)pxl_gray[i][j];
      if (rdum0<pxl0) pxl_tmp[i][j]=0;
      else if (rdum0>=pxl0 && rdum0<pxl1) {
	rdum1=(rdum0-pxl0)*255./rdum2;
	if (rdum1>255.) pxl_tmp[i][j]=255;
	else if (rdum1<0.) pxl_tmp[i][j]=0;
	else pxl_tmp[i][j]=(unsigned char)rdum1;
      }
      else pxl_tmp[i][j]=255;
    }
  }

  if (tag1=="") {
    for (i = 0; i < (int)rows; i++) {
      for (j = 0; j < (int)columns; j++) pxl_gray[i][j]=pxl_tmp[i][j];
    }
    cout<<"    Image "<<tag<<" modified."<<endl;
  }
  else copyImg(*this,pxl_tmp,tag1);
  
  for (i=0;i<(int)rows;i++) delete [] pxl_tmp[i];
  for (i=0;i<2;i++) delete [] pxl_histo[i];
  delete [] pxl_tmp; delete [] pxl_histo;   delete [] pxl_1d;
  
  return;
}

/************************************************************************/
void img::writeImage(string fname, string fmt, int nColor,
		     double *colorBg,double *colorScale)
{ /* Write pxl_gray to an image fname. Each color channel is
     multiplied by  colorScale (between 0 and 1). */
  
  int i,irow, icol,flag;
  string ofname;
  
  if (fmt=="") { // if no format given, check if fname already has an extension
    flag=checkImageNameExt(fname);
    ofname=(flag==0)?(addFileNameExt(fname,"tif")):fname;
  }
  else ofname=addFileNameExt(fname,fmt); //set image file name
  
  for (i=0;i<nColor;i++) { // limit color scale
    colorScale[i]=(colorScale[i]<0.)?0.:colorScale[i];
    colorScale[i]=(colorScale[i]>1.)?1.:colorScale[i];
    
    colorBg[i]=(colorBg[i]<0.)?0.:colorBg[i];
    colorBg[i]=(colorBg[i]>1.)?1.:colorBg[i];    
  }
  
  cout<<"  write:  name= "<<ofname<<"  color scale= ";
  for (i=0;i<nColor;i++) cout <<setw(5)<<setprecision(2)<<colorScale[i]<<" ";
  cout<<endl;
  cout<<"          color_bg= ";
  for (i=0;i<nColor;i++) cout <<setw(5)<<setprecision(2)<<colorBg[i]<<" ";
  cout<<endl;
  
  /* Use Magick++ functions to write image name of extension fmt. 
     Reference: https://www.imagemagick.org/Magick++/Pixels.html */
  double rdum[nColor];
  Image img0(Geometry(0,0),"black");
  
  if (nColor==3) img0.extent(Geometry(columns,rows),
   			     ColorRGB(0.,0.,0.));
  else img0.extent(Geometry(columns,rows),ColorGray(0.));
  
  //image attributes
  img0.type(TrueColorType);
  img0.modulusDepth(BitsPerSample);
  img0.modifyImage(); 
  
  PixelPacket *pixel_cache=img0.getPixels(0,0,columns,rows);
  for (irow=0;irow<(int)rows;irow++) { 
    for (icol=0;icol<(int)columns;icol++) {
      PixelPacket *pixel=pixel_cache+irow*columns+icol;
      for (i=0;i<nColor;i++) {
  	rdum[i]=(colorScale[i]*(double)pxl_gray[irow][icol])/range;
      }
      if (nColor==3) *pixel=ColorRGB(rdum[0],rdum[1],rdum[2]);
      else {
	if (rdum[0]==0) *pixel=ColorGray(*colorBg);
	else *pixel=ColorGray(rdum[0]);
      }
    }
  }
  img0.syncPixels();   //save changes to image
  img0.write(ofname); //wirte image to file   

  return; 
}

/************************************************************************/
// Non-member function
/************************************************************************/

void checkDuplicateImgTag(vector<img> * vimg,string tag)
{ /* In vimg, check if user inputted an already defined image tag */ 
  vector<img>::iterator it_img;
  for (it_img=(*vimg).begin();it_img!=(*vimg).end();++it_img) {
    if ( (*it_img).tag==tag) error_dodri("Duplicate image tag: " + tag);
  }
  return;
}

/************************************************************************/
void copyImg(const img &img0, unsigned char **pxl, string tag1)
{ /* Copy img0 to a new image with tag1, with pixel data pxl. If tag1
     exists, that image is overwritten. If img with tag1 exists,
     pxl[rows][columns] must have the same dimension as tag1's
     pxl_gray[rows][columns] (no check made).
     
     NOTE: Cannot call copyImg more than once within an img member function
     since vimg changes.     
  */
  
  vector<img>::iterator it_img= findImg(tag1);
  int i,j;
  cout<< "  Image "<<img0.tag<<" being copied to "<<tag1<<endl;

  if ( it_img != vimg.end()) {
    cout<<"  WARNING: Image " + tag1 + " is overwritten."<<endl;
    // The following operations are the same as constructor
    // img::img(img *img, unsigned char **pxl, string tag0)
    it_img->image_name = img0.image_name;
    if ((it_img->rows!=img0.rows)||(it_img->columns!=img0.columns))
      error_dodri("Dimensions of the two images do not match.");
    it_img->col_tot = img0.col_tot;
    it_img->BitsPerSample = img0.BitsPerSample;
    it_img->SamplesPerPixel = img0.SamplesPerPixel;
    it_img->pxl_raw = img0.pxl_raw;
    //it_img->pxl_gray = pxl; CANNOT DO THIS WAY - the originally allocated
    // memory for it_img->pxl_gray will be leaked. Do element-wise copy
    for (i=0;i<(int)it_img->rows;i++) {
      for (j=0;j<(int)it_img->columns;j++) it_img->pxl_gray[i][j]=pxl[i][j];
    }
    it_img->pxl_flag = img0.pxl_flag;
    it_img->dRow = img0.dRow;
    it_img->dCol = img0.dCol;
    it_img->ncolors = img0.ncolors;
    it_img->range = img0.range;
  } // if ( it_img != vimg.end()) {
  else { // new image
    vimg.push_back(img(img0, pxl, tag1));
  }
  return;
}

/************************************************************************/
vector<img>::iterator findImg(string tag)
{ /* In vimg, find iterator for the element with tag */ 
  vector<img>::iterator it_img;
  for (it_img=vimg.begin();it_img!=vimg.end();++it_img) {
    if ( (*it_img).tag==tag) break;
  }
  return it_img;
}

/************************************************************************/
void getImgTagRange(stringstream& ss, string *svec)
{
  string sdum,sdum0,sdum1;
  size_t found;  char *cdum;
  int i,j,k;
  
  svec[0]=vimg.begin()->tag; svec[1]=vimg.rbegin()->tag; // default

  sdum=getCommOptString(ss,"img","");
  if (sdum!="") { // range provided
    found=sdum.find(":");
    if (found!=string::npos) { 
      svec[0]=sdum.substr(0,found);
      svec[1]=sdum.substr(found+1);
      if (svec[0]=="") svec[0]=vimg.begin()->tag; // max range for missing part
      if (svec[1]=="") svec[1]=vimg.rbegin()->tag;
      if (findImg(svec[0])==vimg.end())
	error_dodri("tag "+svec[0]+" doesn't exist.");
      if (findImg(svec[1])==vimg.end())
	error_dodri("tag "+svec[1]+" doesn't exist.");
    }
    else error_dodri("img keyword must be followed by tag_ini:tag_fin.");
  }
  sdum=getCommOptString(ss,"img_id","");
  if (sdum!="") { // range provided
    if (checkCommOptKey(ss,"img",0)==1) {
      cout<<"  Both img and img_id given, ignoring img_id."<<endl;
      return;
    }
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
      svec[0]=vimg[i].tag;
      svec[1]=vimg[j].tag;
      cout<<"  img id_range= "<<i<<" "<<j<<endl;
    }
    else error_dodri("img_id keyword must be followed by id_ini:id_fin.");
  }
  cout<< "  img tag_range= "<<svec[0]<<" "<<svec[1]<<endl;
  return;
}

/************************************************************************/
void mergeImg(img *img0, img *img1, string mode,string tag2)
{ /* Merge img0 and img1. If tag2="", save to img0. 
     If tag2 given create new img, img0 and img1 will be unmodified. */
  
  int i,j,k=0;
  double avgdum,sddum,rdum;
  unsigned char **pxl_tmp;
  vector<img>::iterator it_img,it_img1;
  string sdum=(tag2=="")?img0->tag:tag2; 

  cout<<"  merge: " + img0->tag + " with "+img1->tag +"  mode= " << mode
      << "  save= " + sdum << endl;
  
  if ((img1->rows!=img0->rows)||(img1->columns!=img0->columns))
    error_dodri("Size of two images do not match.");

  avgdum=sddum=0; 
  if (tag2=="") tag2=img0->tag; // merge into img1 by default
  it_img= findImg(tag2);
  if (it_img != vimg.end()) {
    if ((it_img->rows!=img1->rows)||(it_img->columns!=img1->columns))
      error_dodri("Destination image has different dimension.");
    // The following operations are the same as constructor
    // img::img(img *img, unsigned char **pxl, string tag0)
    it_img->image_name = img1->image_name;
    it_img->col_tot = img1->col_tot;
    it_img->BitsPerSample = img1->BitsPerSample;
    it_img->SamplesPerPixel = img1->SamplesPerPixel;
    it_img->TotalColors = img1->TotalColors; 

    it_img->pxl_raw = img1->pxl_raw;
    it_img->pxl_flag = img1->pxl_flag;
    it_img->dRow = img1->dRow;
    it_img->dCol = img1->dCol;
    it_img->ncolors = img1->ncolors;
    
    for (i=0;i<(int)it_img->rows;i++) {
      for (j=0;j<(int)it_img->columns;j++) {
	if (mode == "add") 
	  k=(int)img0->pxl_gray[i][j]+(int)img1->pxl_gray[i][j];
	else if (mode == "subtract") 
	  k=(int)img0->pxl_gray[i][j]-(int)img1->pxl_gray[i][j];
	if (k>255) k=255;
	if (k<0) k=abs(k);
	it_img->pxl_gray[i][j]=(unsigned char)k;
	rdum=(double)k; avgdum+=rdum; sddum+=rdum*rdum;
      }
    }
    cout<<"    Image "<<tag2<<" modified."<<endl;
  } // if ( it_img != vimg.end()) {
  else { // new image
    pxl_tmp=new unsigned char*[img1->rows];
    for (i=0;i< (int)img1->rows;i++)
      pxl_tmp[i]=new unsigned char[img1->columns];
    for (i=0;i< (int)img1->rows;i++) {
      for (j=0;j< (int)img1->columns;j++) {
	if (mode == "add")
	  k=(int)(img0->pxl_gray[i][j])+(int)(img1->pxl_gray[i][j]);
	if (mode == "subtract") 
	  k=(int)(img0->pxl_gray[i][j])-(int)(img1->pxl_gray[i][j]);
	if (k>255) k=255;
	if (k<0) k=abs(k);
	pxl_tmp[i][j]=(unsigned char)k;
      }
    }
    vimg.push_back(img(*img1, pxl_tmp, tag2));
    it_img1= findImg(tag2);

    cout<<"  New image "<<tag2<<" created."<<endl;
    for (i=0;i< (int)img1->rows;i++) delete [] pxl_tmp[i];
    delete [] pxl_tmp;
  }
  return;
}
