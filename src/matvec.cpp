#include "dodri.h"
/* matvec.cpp: General matrix or vector-related functions. */ 

/************************************************************************/
int checkBondCrossing(double *p1,double *p2,double *p3,double *p4,double *p0)
{ /* check if line segments connecting points p1-p2 and p3-p4 cross.
     return value: 1: cross, 0: do not cross 
     p0: coords of crossing point saved. */

  double x1=p1[0],y1=p1[1],x2=p2[0],y2=p2[1];
  double x3=p3[0],y3=p3[1],x4=p4[0],y4=p4[1];
  double x0,y0,rdum1,rdum2;
  double *iVec=new double[2], *jVec=new double[2];

  /* x0,y0: crossing point. */
  x0=(y3-y1+(y2-y1)*x1/(x2-x1)-(y4-y3)*x3/(x4-x3))
    /((y2-y1)/(x2-x1)-(y4-y3)/(x4-x3));
  if (x2-x1==0) x0=x1; //vertical line

  if (fabs(x2-x1)>fabs(x4-x3)) // increase calc accuracy
    y0=(y2-y1)*(x0-x1)/(x2-x1)+y1;
  else 
    y0=(y4-y3)*(x0-x3)/(x4-x3)+y3;
  // If bonds cross, angles from x0 to the ends of bonds are both pi 
  iVec[0]=x1-x0; iVec[1]=y1-y0; jVec[0]=x2-x0; jVec[1]=y2-y0;

  rdum1=getAngle(iVec,jVec,0);
  iVec[0]=x3-x0; iVec[1]=y3-y0; jVec[0]=x4-x0; jVec[1]=y4-y0;
  rdum2=getAngle(iVec,jVec,0);
  delete [] jVec; delete [] iVec;
  if ((rdum1>(0.5*PI))&&(rdum2>(0.5*PI))) {
    p0[0]=x0; p0[1]=y0; return 1;
  }
  else return 0;
}

/************************************************************************/
void crossProd(double a[], double b[], double *c)
{ /* Calculates \vec a \times \vec b in 3D. No size check. */

  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];
  
  return; 
}

/************************************************************************/
void distPoint2Line(double v[], double u[], double q[], int ndim, double *p, 
		      char flag)
{ /* In ndim-dim, calculates distance from \vec v to \vec r=p[0]*\vec u+\vec q.
   p[1]:  flag==1: distance,   flag==2: distance^2 
          flag==0: Skip distance calculation and p[1]=0.

  NOTE: u must be a unit vector (not checked here)  */
  
  int i;
  double *vdum=new double[ndim], p0,rdum;

  for (i=0;i<ndim;i++)  vdum[i]=v[i]-q[i];
  
  p0=dotprod(vdum,u,ndim); p[0]=p0;
  if (flag==0) rdum=0.;
  else {
    rdum=dotprod(vdum,vdum,ndim);
    if ((flag==1)||(flag==2)) { // distance
      rdum=rdum-p0*p0;
      if (rdum<0.) rdum=0.;
      else { rdum=(flag%2==1)?rdum:sqrt(rdum);}
    }
    else error_dodri("distPoint2Line: flag "+to_string((int)flag)
		    +" out of range.");
  }
  p[1]=rdum;
  
  delete [] vdum;
  return;
}

/************************************************************************/
double dotprod(double a[], double b[], int N)
{ /* Calculates \vec a\cdot \vec b in N-dimension
     Note that this routine does not check if a & b are of the same size */
  
  double r=0.0;  
  for (int i=0;i<N;i++) r+=a[i]*b[i];
  return r;
}

/************************************************************************/
void fit0(double x[], double y[], int ndata, double *a,
                double *b, double *siga, double *sigb, double *chi2)
{
  /* Given a set of data points x[ndata],y[ndata] fit them to a straight
  line y=a+bx by minimizing chi^2. Returned are a,b and their respective
  probable uncertainties siga and sigb, and the chi-square chi2.  */

  int i;
  double t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;
  *b=0.0;
  for (i=0;i<ndata;i++) {  //  ...or without weights.
    sx += x[i];
    sy += y[i];
  }
  ss=ndata;
  sxoss=sx/ss;
  for (i=0;i<ndata;i++) {
    t=x[i]-sxoss;
    st2 += t*t;
    *b += t*y[i];
  }

  *b /= st2;   // Solve for a, b,  a, and  b.
  *a=(sy-sx*(*b))/ss;
  *siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
  *sigb=sqrt(1.0/st2);
  
  *chi2=0.0;  //  Calculate  chi^2.

  for (i=0;i<ndata;i++) *chi2 += pow((y[i]-(*a)-(*b)*x[i]),2.);
  sigdat=sqrt((*chi2)/(ndata-2));  // For unweighted data evaluate 
  *siga *= sigdat;                 // typical sig using chi2, and
  *sigb *= sigdat;                 // adjust the standard deviations.
}

/**********************************************************************/
void fit1(double x[], double y[], int ndata, double weight[], int mwt,
	 double *a, double *b, double *siga, double *sigb, double *chi2)
{
  /* Given a set of data points x[ndata],y[ndata] with individual
     weight[ndata], fit them to a straight line y=a+bx by minimizing
     chi^2. Returned are a,b and their respective probable uncertainties
     siga and sigb, and the chi-square chi2.
     
     mwt=1: use weight. mwt=0: don't use weight[].
     weight[] is like 1/sig[] in the original routine in NR.
     
     The following not implemented: the goodness-of-fit probability q (that
     the fit would have chi^2 this large or larger). If mwt=0 on input,
     then the std are assumed to be unavailable: q is returned as 1.0 and
     the normalization of chi^2 is to unit std on all points.  */

  int i;
  double wt, t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;
  double sw=0.; // sum of weights
  *b=0.0;
  if (mwt) { 
    ss=0.;
    for (i=0;i<ndata;i++) {
      wt=pow(weight[i],2.0);
      sw+=wt;
      ss+=wt; sx+=x[i]*wt; sy+=y[i]*wt;
    }
  } else {
    for (i=0;i<ndata;i++) {  //  ...or without weights.
      sx += x[i];
      sy += y[i];
    }
    ss=ndata;
  }
  sxoss=sx/ss;
  if (mwt) { 
    for (i=0;i<ndata;i++) {
      t=(x[i]-sxoss)*weight[i];
      st2+=t*t;
      *b += t*y[i]*weight[i];
    }
  } else {
    for (i=0;i<ndata;i++) {
      t=x[i]-sxoss;
      st2 += t*t;
      *b += t*y[i];
    }
  }
  *b /= st2;   // Solve for a, b,  a, and  b.
  *a=(sy-sx*(*b))/ss;
  *siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
  *sigb=sqrt(1.0/st2);
  
  *chi2=0.0;  //  Calculate  chi^2.

  if (mwt==0) { 
    for (i=0;i<ndata;i++) *chi2 += pow((y[i]-(*a)-(*b)*x[i]),2.);
    sigdat=sqrt((*chi2)/(ndata-2));  // For unweighted data evaluate 
    *siga *= sigdat;                 // typical sig using chi2, and
    *sigb *= sigdat;                 // adjust the standard deviations.
    //    sw=(double)ndata;
  } else {
    for (i=0;i<ndata;i++) *chi2 += pow((y[i]-(*a)-(*b)*x[i])*weight[i],2.);
  }
  *chi2 /= sw; // normalize chi^2
}

/************************************************************************/
void fitPlane(int N, double arr[], double *pvec)
{ /* Given a set of N/3 points (=npts) data points, fit them to a plane
     z=ax+by+c. Returned are a, b, and c in pvec.    */

  int i,j,npts=N/3; 
  double sum, **pt=new double*[3]; //3d coord (x,y,z)
  double **alpha=new double*[3];
  double **beta=new double*[3];
  
  for (j=0;j<3;j++) {pt[j]=new double[npts];}
  for (i=0;i<npts;i++) {
    for (j=0;j<3;j++) {
      if (arr[i*3+j]==-1)
	error_dodri("  Need N*3 reference points. " );
      pt[j][i]=arr[i*3+j];}
  }

  for (i=0;i<3;i++) {alpha[i]=new double[3];beta[i]=new double[1];}
  
  // fill matrix
  sum=0.; for (i=0;i<npts;i++) {sum+=pow(pt[0][i],2);}
  alpha[0][0]=sum;    //alpha[0][0]=x^2
  sum=0.;for (i=0;i<npts;i++) {sum+=(pt[0][i]*pt[1][i]);}
  alpha[1][0]=sum;alpha[0][1]=sum;   //alpha[1][0]=xy
  sum=0.;for (i=0;i<npts;i++) {sum+=(pt[0][i]);}
  alpha[2][0]=sum; alpha[0][2]=sum;   //alpha[2][0]=x
  sum=0.;  for (i=0;i<npts;i++) {sum+=pow(pt[1][i],2);}
  alpha[1][1]=sum;   //alpha[1][1]=y^2
  sum=0.;  for (i=0;i<npts;i++) {sum+=(pt[1][i]);}
  alpha[2][1]=sum; alpha[1][2]=sum;   //alpha[2][1]=y
  sum=0.;  for (i=0;i<npts;i++) {sum+=1;}
  alpha[2][2]=sum;   //alpha[2][2]=N
  
  sum=0.;  for (i=0;i<npts;i++) {sum+=(pt[0][i]*pt[2][i]);}
  beta[0][0]=sum;   //beta[0][0]=z*x
  sum=0.;  for (i=0;i<npts;i++) {sum+=(pt[1][i]*pt[2][i]);}
  beta[1][0]=sum;   //beta[1][0]=z*y
  sum=0.;  for (i=0;i<npts;i++) {sum+=(pt[2][i]);}
  beta[2][0]=sum;   //beta[2][0]=z

  gaussj(alpha,3,beta,1);   //solve matrix with Gauss-Jordan method
  for (i=0;i<3;i++) {pvec[i]=beta[i][0];}
  
  for (j=0;j<3;j++) { delete[] pt[j];delete[] alpha[j];delete[] beta[j];}
  delete[] pt; delete[] alpha; delete[] beta; 
  
  return;
}

/************************************************************************/
void fitQuad0(double *p1, double *p2, double *p3, double *vel, 
		   double *acc, int N)
{ /* For 3 vectors p1, p2, p3 (dimension N), get velocity vel and acceleration
    acc so that a general position pos= p1+vel*s+(1/2)s^2*acc, where 
    s is the contour length. To get vel & acc, approximate s to be 
    the straight distance s2=|p2-p1| for p2, and s3=|p3-p2|+|p2-p1|.
    Plug in s2 & s3 to the eq for pos, so that
    
    p2=p1+vel*s2+(1/2)s2^2*acc,    p3=p1+vel*s3+(1/2)s3^2*acc

    then solve for vel & acc.  */
  
  int i;
  double s2, s3, rdum;
  s2=s3=0;
  for (i=0;i<N;i++) {
    rdum=(p2[i]-p1[i]); s2+=(rdum*rdum);
    rdum=(p3[i]-p2[i]); s3+=(rdum*rdum);
  }
  s2=sqrt(s2); s3=sqrt(s3)+s2;
  for (i=0;i<N;i++) {
    acc[i]=(2/(s2*s3))*(p1[i]-p2[i]+s2*(p2[i]-p3[i])/(s2-s3));
    vel[i]=(p2[i]-p3[i])/(s2-s3) - 0.5*acc[i]*(s2+s3);
  }
  return;
}

/************************************************************************/
int gaussj(double **a, int n, double **b, int m)
{ /* Solve matrix a using Gauss-Jordan elimination. */
  
  int *indxc, *indxr, *ipiv;
  int i, icol=0, irow=0, j,k,l,l1;
  double big,dum,pivinv,temp;
  // indxc[i] : column of the i-th pivot element. 
  // indxr[i] : the row in which above pivot was originally located.
  
  indxc=new int[n];
  indxr=new int[n];
  ipiv =new int[n];
  for (i=0;i<n;i++) ipiv[i]=0;

  for (i=0;i<n;i++) { // Main loop over columns to be reduced
    big=0.;
    for (j=0;j<n;j++)
      if (ipiv[j] != 1)
        for (k=0;k<n;k++) {
          if (ipiv[k]==0) {  
            if ( fabs(a[j][k]) >=  big) {
              big=fabs(a[j][k]);
              irow=j; icol=k;
            }
          }
          else if (ipiv[k]>1)
            return -1; //cout<<"Singular Matrix! "<<endl;
	}
    ++(ipiv[icol]);
    // Now that we have the pivot element, interchange rows to put the
    // pivot element on the diagonal. The columns are not physically 
    // interchanged, only relabeled: indxc[i], the column of the i-th
    // pivot element is the i-th column that is reduced, while 
    // indxr[i] is the row in which that pivot element was originally
    // located. If indxr[i].ne.indxc[i], there is an implied column 
    // interchange. With this bookkeeping, the solution b's will end up in
    // the correct order, and the inverse matrix will be scrambled by
    // columns.
    if (irow!=icol) {  // Interchange two rows.    
      for (l=0;l<n;l++) { //taking care of SWAP(a[irow][l],a[icol][l]);
	temp=(a[irow][1]);
	(a[irow][1])=(a[icol][1]);
	(a[icol][1])=temp;
      }
      for (l=0;l<m;l++) { //SWAP(b[irow][l],b[icol][l]);
	temp=(b[irow][1]);
	(b[irow][1])=(b[icol][1]);
	(b[icol][1])=temp;
      }
    }
    indxr[i]=irow; 
    indxc[i]=icol;
    // Normalize pivot row
    if (a[icol][icol]==0.) 
      return -2; //cout<<"singular matrix-2"<<endl;
    pivinv=1./a[icol][icol];
    a[icol][icol]=1.;
    for (l=0;l<n;l++) a[icol][l] *=pivinv;
    for (l=0;l<m;l++) b[icol][l] *=pivinv;
    // Reduce the rows
    for (l1=0;l1<n;l1++)
      if (l1 != icol) {
        dum=a[l1][icol];
        a[l1][icol]=0.;
        for (l=0;l<n;l++) a[l1][l] -= a[icol][l]*dum;
        for (l=0;l<m;l++) b[l1][l] -= b[icol][l]*dum;
      }
  }  // End of Main loop over columns (i).
  // Now undo the column interchange
  for (l=(n-1);l>=0;l--) {
    if (indxr[l] != indxc[l])
      for (k=0;k<n;k++) {  //  SWAP(a[k][indxr[l]], a[k][indxc[l]]);
	temp=(a[k][indxr[1]]);
	(a[k][indxr[1]])=(a[k][indxc[1]]);
	(a[k][indxc[1]])=temp;
      }
  }
  delete [] ipiv; delete [] indxr; delete [] indxc;
  return 1; // Operation successful 
}

/************************************************************************/
double getAngle(double *bond_1, double *bond_2, char flag)
{ /* Find angle between two bonds bond_1 & bond_2. 
     flag={0,1}: 2-Dim angle. 
        flag=0: Range is (0,PI)
        flag=1: Range is (-PI,PI). The sign is defined in the
           counterclockwise manner from bond_1 to bond_2.  
     flag=3: 3-dim angle  */
  
  int i, ndim=-1;
  double r1, r2, r3, phi;
  double *bvec1, *bvec2;
  if ((flag==0)||(flag==1)) ndim=2;
  else if (flag==3) ndim=3;
  else {
    stringstream ss; ss << "getAngle: Unknown flag value: "<< flag;
    error_dodri(ss.str());
  }
  
  bvec1=new double[ndim];  bvec2=new double[ndim];  
  for (i=0;i<ndim;i++) { bvec1[i]=bond_1[i]; bvec2[i]=bond_2[i]; }
  r1=dotprod(bvec1,bvec1,ndim);  r2=dotprod(bvec2,bvec2,ndim);

  r1=sqrt(r1); r2=sqrt(r2);
  for (i=0;i<ndim;i++) { bvec1[i]/=r1; bvec2[i]/=r2; }
  r3=dotprod(bvec1,bvec2,ndim);
  if (r3>1.0) r3=1.0-1e-32; // take care of round off error
  if (r3<-1.0) r3=-1.0+1e-32; // take care of round off error
  phi=acos(r3);
  if (flag==1) { // get the sign 
    r3=bvec1[0]*bvec2[1]-bvec1[1]*bvec2[0]; // cross product
    if (r3<0) phi*=-1.0;
  }
  delete[] bvec1;  delete[] bvec2;
  return phi; // return the magnitude
}

/************************************************************************/
double getavg0(double ia[], int N)
/* Calculate average of an array ia[] of size N (no std calculated). */
{
  int i; double a;
  a=0;
  for (i=0;i<N;i++) {a+=ia[i];}
  a/=(double)N;
  return a;
}

/************************************************************************/
void getavgw(double *avg, double *sig, double ia[], double weight[], int N)
/* Calculate weighted average and s.d. of an array ia[] of size N. */
{
  int i; double a, b, wt=0., rdum;
  a=b=0;
  for (i=0;i<N;i++) { 
    wt+=weight[i];  rdum=ia[i];
    a+=rdum*weight[i]; b+=rdum*rdum*weight[i];
  }
  a/=wt; b/=wt;
  b=b-a*a; if (fabs(b)<TOL) b=0.; // take care of roundoff
  assert(b>=0.);
  *avg=a; *sig=sqrt(b);
  return;
}

/************************************************************************/
void getCentroid(double **pos, double *cm, int N, int ndim )
{ /* Calculate centroid of an array pos[ndim][N] and store in cm[ndim] */

  for (int i=0;i<ndim;i++) {
    cm[i]=getavg0(pos[i],N);
  }
  return ; 
}

/************************************************************************/
double getDistPos(double *pos1, double *pos2, char flag) 
{ /* Get distance between pos1 & pos2:
     returns distance (flag=0) distance squared (flag!=0) */
  double r, x1,y1,x2,y2;
  x1=pos1[0]; y1=pos1[1]; x2=pos2[0]; y2=pos2[1];
  r=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
  if (flag==0) return sqrt(r);  else return r;
}

/************************************************************************/
double getDistPosNdim(int dim,double pos1[], double pos2[], char flag) 
{ /* Get Eucledian distance between N-dimensional points pos1 & pos2:
     returns distance (flag=0) distance squared (flag!=0) */ 
  double r=0.;
  for (int i=0;i<dim;i++) r+=pow(pos1[i]-pos2[i],2);
  if (flag==0) return sqrt(r);  else return r;
}

/************************************************************************/
void getNormal(double **d0, double *norm)
{ /* Calc unit normal from any two 3D vectors. 
     Note: d0[x,y,z][vertex]. */
  
  double Nx, Ny, Nz, N,vec[2][3];
  
  for (int j=0;j<3;j++) {
    vec[0][j]=d0[j][1]-d0[j][0];
    vec[1][j]=d0[j][2]-d0[j][0];
  }
  
  Nx = vec[0][1]*vec[1][2]-vec[0][2]*vec[1][1];
  Ny = vec[0][2]*vec[1][0]-vec[0][0]*vec[1][2];
  Nz = vec[0][0]*vec[1][1]-vec[0][1]*vec[1][0];
  N = sqrt(Nx*Nx+Ny*Ny+Nz*Nz);
  norm[0]=Nx/N; norm[1]=Ny/N ; norm[2]=Nz/N;
}

/************************************************************************/
void line2vec(double a, double b, double *u, double *q, signed char flag)
{ /* Convert an equation for a line to the form \vec r=p\hat u +\vec q
   where p is the contour length parameter.
   flag>0: use y=a+bx
   flag<0: use x=a+by
   
   abs(flag)==1: \vec q=(0,a) (flag>0), (a,0) (flag<0)
   abs(flag)==2: \vec q \bot \vec u   */
  
  double bdum=1/(sqrt(1+b*b)), adum=a*bdum*bdum;
  if (flag>0) { // y=a+bx
    u[0]=bdum; u[1]=b*bdum;
    if (flag%2==1) {q[0]=0; q[1]=a;}
    else { q[0]=-b*adum; q[1]=adum; }
  }
  else if (flag<0) { // x=a+by
    u[0]=b*bdum; u[1]=bdum;
    if (flag%2==1) {q[0]=a; q[1]=0;}
    else { q[0]=adum; q[1]=-b*adum; }
  }
  else error_dodri("line2vec: Flag should be nonzero.");
}

/************************************************************************/
void rotate2dim(double *x, double *y, int N, double phi)
{ /* Rotate x[N],y[N] by angle phi */
  int i;
  double x1,y1;
  double cp=cos(phi),sp=sin(phi), sp1=-1.0*sp;
    for (i=0;i<N;i++) {
    x1=x[i]*cp+y[i]*sp; y1=x[i]*sp1+y[i]*cp;
    x[i]=x1; y[i]=y1;
  }
  return;
}

/************************************************************************/
int verifyTriRHR(double v[3][3],double n[3])
{ /* Verify triangle vertices satisfy RHR condition to normal. 
     If not, order needs to be adjusted (not done here). 

     RHR condition: triangle vertices listed CCW to triangle normal 
                    pointing "out of" surface

     Input: triangle vertices, surface normal oriented "out" of surface
     Output: vertex order sequence    */ 
  
  int k;
  double phi,e1[3],e2[3],ex[3]; 
  
  for (k=0;k<3;k++) { e1[k]=v[1][k]-v[0][k];  e2[k]=v[2][k]-v[0][k];  }  
  crossProd(e1,e2,ex);
  phi=getAngle(ex,n,3); // phi: CW==PI; CCW==0
  
  if (phi>PI/2) return 0; // v2 and v1 need to be switched (not done here)
  else return 1; // keep order
}
