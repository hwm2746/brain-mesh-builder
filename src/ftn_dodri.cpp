#include "dodri.h"

/************************************************************************/
string addFileNameExt(string fname,string fmt)
{ /* Add fmt as a filename extension. 
     fmt may or may not have a dot in front. */
  
  string newfname;
  std::size_t pos;
  if (fmt.find(".")==string::npos) fmt="."+fmt; // add a dot 
  pos=fname.find(fmt); //fmt: dat, psf, tif, cor, etc
  if ( pos == string::npos) newfname=fname+fmt;
  else { // ensure extension is at the end of the filename
    newfname = (fname.length()- pos != fmt.length())?(fname+fmt):fname;
  }
  return newfname;
}

/************************************************************************/
int checkImageNameExt(string fname)
{ /* Check if fname contains a common image filename extension */

  vector<string>
    extName{"bmp","eps","gif","jpg","jpeg","png","pdf","pnm","ppm","ps","tga","tif"};
  vector<string>::iterator it;
  for (it=extName.begin();it!=extName.end();it++) {
    if (fname.find("."+(*it))!=string::npos) return 1; // extension found
  }
  return 0; // extension absent
}

/************************************************************************/
void error_dodri(string sdum)
{ 
  cout << "ERROR: "<<sdum << endl; exit(-1);
}

/**********************************************************************/
int histogram(double **pdf, double *a, int N, int nBin, double *dBin, 
		   double *amin, double *amax, char normalize)
{ /* Build histogram of size nBin for a[N]. dBin: size of each bin. 
     if amin<amax: histogram covers range of given amin/amax
     if amin<0 or amax<0: respective min/max range is found in a[].
     normalize=0: do not normalize. !=0: normalize.
     returns number of pixels considered.
     
     Ends in the range of histogram are
     min and max of a[] (no padding of dBin/2 at ends.
     histogram[0][0]=min of a[], histogram[0][nBin-1]={(max of a[])-dBin}.
     histogram[0][nBin]: ordinate, histogram[1][nBin]: count
  */
  
  double dr,rdum,rdum1;
  double amin0=1.0e308, amax0=-1.0e308;
  int i,j,cnt;
  
  for (i=0;i<2;i++) {
    for (j=0;j<nBin;j++) pdf[i][j]=0.;
  }
  if (*amin< 0.) { // find min
    for (i=0;i<N;i++) if (a[i]<amin0) amin0=a[i];
    *amin=amin0;
  }
  else {amin0=*amin;}
  if (*amin> *amax) { // find max
    for (i=0;i<N;i++) if (a[i]>amax0) amax0=a[i];
    *amax=amax0;
  }
  else {amax0=*amax;}
  
  dr=(amax0-amin0)/(double)nBin; *dBin=dr;
  for (i=0;i<nBin;i++) {pdf[0][i]=amin0+dr*(double)i;}
  cnt=0;
  for (i=0;i<N;i++) {
    for (j=0;j<nBin;j++) { // include amin0/amax0 in the range
      rdum=(j==0)?(pdf[0][j]-(fabs(amin0)*1.0e-64)):pdf[0][j]; 
      rdum1=(j==(nBin-1))?(amax0+(fabs(amax0)*1.0e-64)):pdf[0][j+1];
      if ((a[i]>=rdum)&&(a[i]<rdum1)) { pdf[1][j]+=1.0; ++cnt;}
    }
  }
  if (normalize!=0) {
    rdum=1.0/(double)cnt;
    for (i=0;i<nBin;i++) pdf[1][i]*=rdum;
  }
  return cnt;
}

/************************************************************************/
string itostr (int i)
{ // integer to string
  stringstream ss;
  ss << i;
  return ss.str();
}

/************************************************************************/
void MinMax(int N, double *a, int *imin, int *imax, double *min, double *max)
{ /* Find min and max values for array a of size N.
     imin/imax: indexes of min and max */
  
  double amin=1.0e308, amax=-1.0e308;
  int i;
  
  for (i=0;i<N;i++) {
    if (a[i]<amin) {amin=a[i]; *imin=i;}
    if (a[i]>amax) {amax=a[i]; *imax=i;}
  }
  *min=amin; *max=amax;
}

/************************************************************************/
int nDigit(int i)
{ /* Count number of digits of i */ 
  stringstream ss; string sdum;
  ss.str("");ss.clear(); ss << i; sdum=ss.str();
  return sdum.length();
}
