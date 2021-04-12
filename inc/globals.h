/* globals.h: Global parameters and variables */

#ifdef SET_EXT
  #define EXT_NAME 
#else
  #define EXT_NAME extern
#endif

//#include <algorithm> 
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
//#include <cstdio>
#include <cmath>
#include <iomanip>
//#include <iterator>
#include <string>
#include <cstring>
//#include <cstdlib>
//#include <ctime>
//#include <list>
#include <map>
#include <vector>
#include <numeric>

//magick++ related
#include <Magick++.h> 

using namespace Magick; 
using namespace std;

#ifndef DODRI
#define DODRI

// Global Parameters
#define nPIXEL 3 // 3 for rgb, 1 for grayscale output
#define maxVec 200 // max number of bonds to consider from one bead 
#define PI 3.14159265358979323846264338
#define degree2rad 0.017453292519943295 // pi/180
#define rad2degree 1/degree2rad // 180/pi
#define TOL 1.0e-5 // tolerance for rounding off a double number

/************************************************************************/
// Global variables

/* nCell: number of cells subdiving an image
   cellMark: for neighbor search
   nstack: number of image stacks. For single-image analysis, nstack=0 
*/
EXT_NAME int nCell, cellMark[9]; 
EXT_NAME clock_t timer; 

/************************************************************************/
// Global functions

extern string addFileNameExt(string fname,string fmt);
extern int checkBondCrossing(double *p1,double *p2,double *p3,double 
			     *p4,double *p0);
extern void crossProd(double a[], double b[], double *c); 
extern void distPoint2Line(double v[], double u[], double q[], int ndim,
			   double *p, char flag);
extern double dotprod(double a[], double b[], int N);
extern void error_dodri(string sdum);
extern void fit0(double x[], double y[], int ndata, double *a,
                double *b, double *siga, double *sigb, double *chi2);
extern void fit1(double x[], double y[], int ndata, double weight[], int mwt,
	  double *a, double *b, double *siga, double *sigb, double *chi2);
extern void fitPlane(int N, double arr[], double *pvec);
extern void fitQuad0(double *p1, double *p2, double *p3, double *vel, 
			  double *acc, int N);
extern int gaussj(double **a, int n, double **b, int m);
extern double getAngle(double *bond_1, double *bond_2, char flag);
extern double getavg0(double ia[], int N); 
extern void getavgw(double *avg, double *sig, double ia[], double weight[],
		    int N);
extern void getCentroid(double **pos, double *cm, int N, int ndim);
extern double getDistPos(double *pos1, double *pos2, char flag);
extern double getDistPosNdim(int dim,double pos1[], double pos2[], char flag);
extern void getNormal(double **d0, double *norm);
extern int histogram(double **pdf,double *a, int N, int nBin, double *dBin, 
		     double *amin, double *amax, char normalize);
extern string itostr (int i);
extern void line2vec(double a, double b, double *u, double *q,
		     signed char flag);
extern void MinMax(int N, double *a, int *imin, int *imax, double *min, 
		   double *max);
extern int nDigit(int i);
extern void rotate2dim(double *x, double *y, int N, double phi);
extern int verifyTriRHR(double v[3][3],double n[3]); 

/************************************************************************/
inline std::string& ltrim(std::string& s, const char* t = " \t\n\r\f\v")
{
  s.erase(0, s.find_first_not_of(t));
  return s;
}

/************************************************************************/
inline std::string& rtrim(std::string& s, const char* t = " \t\n\r\f\v")
{
  s.erase(s.find_last_not_of(t) + 1);
  return s;
}

/************************************************************************/
inline std::string& removeComment(std::string& s)
{ /* Remove comments from i/o line */
  
  s=s.substr(0, s.find("#"));  
  return s;
}

#endif // #ifndef DODRI
