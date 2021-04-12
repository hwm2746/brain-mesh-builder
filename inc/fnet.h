#ifndef FNET
#define FNET

#ifdef SET_EXT
  #define EXT_NAME 
#else
  #define EXT_NAME extern
#endif

#include "bead.h"

/************************************************************************/
class fnet : public bead {
 public:
  int nBond, nFilament; 
  string bead_tag, tag;
  int maxFilLength;

  multimap<int,int> bond,bond_del; 
  multimap<int,int> bead_fil, fil_bead;

  // Constructor
  fnet () {};// constructor for declaration
  fnet(bead *bd, string tag0); // from beads  
  fnet(img *img0,int iDiam,unsigned char pxlCut,string tag0);  // contour
  fnet(string tag, string ifname, string bead_tag0,ifstream *fin=NULL); //from data file

  // Member functions
  void addBead2Fn(double mass,double cmr,double cmc); 
  void assignFilSeqPos(int filSeq[], double x[], double y[], int N,
		       signed char dir);
  double contourLength(int *filSeq, int ipos, int jpos);
  void deleteBead(int iBead, char sort);
  char deleteMark(int p, int q);
  void equalizeFilamentBond(double rCut);
  int findCluster(multimap<int,int> *bead_cluster, 
		  multimap<int,int> *cluster_bead);
  void findFilament(int flag);
  int findFilamentSeq(int i, int *filSequence);
  int getBondVec(int iBead, int iSkip, int *bondVecList);
  double getDist(int i, int j, int flag);  
  void getPxlEdge(int ix0,int iy0,unsigned char pxlCut,
		  vector<pair<double,double>> &edge_pos,img *img0);
  int isBonded(int ibead,int jbead);  
  void listCorners(double x0,double y0,int dir, double dr,
		   vector<pair<double,double>> &corner );
  void joinCyclicFil(double rCut);
  void measureWidth(double arr[],int nslice,string ofname,
		    vector<array<double,2>> &axvec,
		    vector<double> &axwdth,
		    vector<array<double,2>> &epos0,
		    vector<array<double,2>> &epos1,
		    string ax0,double a0) ;
  void moveBeadPos(int p,double *newPos);
  void moveFilSeq(int *filSeq, int filLength, int i_fr, int i_to);
  void putBeadToFilament(int iBead, int iFil);
  void readData(string fname, string bead_tag0, ifstream *fin=NULL) ;
  int registerOneBond(int i, int j, int iCut, double rCutSq, img *img);
  void removeFloatingFil(int length_cut);   
  void removeOneBond(int i, int j);
  void smoothenFilAverage(int segLength,double rCut);
  void updateFilament(int iBead,int kBead); 
  void writeCor(string fname, string fmt, string occ, string outMode);
  void writeData(string fname, string bd_fname,FILE *dat=NULL,
		 FILE *bd_dat=NULL);  
  void writeImage(int beadL, int bondL, int nColor, double *colorBg,
		  double *colorBead, double *colorBond,
		  string fname, string fmt);
  void writePsf(string fname);
};

/************************************************************************/
// Non-member function

void checkDuplicateFnetTag(vector<fnet> * vfnet,string tag);
vector<fnet>::iterator findFnet(vector<fnet> * vfnet,string tag);
//void getFnetTagRange(stringstream& ss, string *svec);
//vector<fnet>::iterator beadToFnet(vector<fnet> *vfnet,string bead_tag);

int fnet_cmd(stringstream& ss, ifstream *fin0);
/************************************************************************/
// Global variable
EXT_NAME  vector<fnet> vfnet;

#endif // #define FNET
