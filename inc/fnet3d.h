#ifndef FNET3D
#define FNET3D

#ifdef SET_EXT
  #define EXT_NAME 
#else
  #define EXT_NAME extern
#endif

#include "fnet.h"

/************************************************************************/
class fnet3d { 
 public:
  int nStack;
  string tag, bead3d_tag,img3d_tag;
  vector<fnet> fnet_stack;
  // bond_a/b: bonds with imaging planes above/below
  multimap<int,int> *bond_a, *bond_b; 
  
  // Constructor
  fnet3d () {};
  fnet3d(bead3d *bd3d, int *id_vec, string tag0); // from bead 
  fnet3d(img3d *img3d0, int iDiam,unsigned char pxlCut,string tag0) ; // from image
  fnet3d(string tag, string ifname, string bead3d_tag0); //from input file  
  
  // Member functions
  void addMissingZBond(int istk,int jstk,
		       vector<int> iSeq, vector<int> jSeq);     
  void assignBondZ(double rCut, int flag, int *id_vec);
  void clearBondZ(int *id_vec) ;
  void equalizeFilamentBond(double rCut, int *id_vec);  
  void fillEnclosedHole(int *id_vec);
  void findFilament(int flag,int *id_vec);
  void findNextZ(int istk,int iBead, int iFil,int kstk,
		 vector<int> &z_bonded);
  void findSeqZ(int iBead0,int istk,multimap<int,int> &zlist,int *id_vec);  
  int fnetTag2Id(string sdum);  
  int getBondVecZ(int p, int istk, int iSkip, double **bondVec,
		  int *bondVecList, int dir);  
  void joinCyclicFil(double rCut, int *id_vec);
  int locateEnclosedHole(vector<array<int,2>> *hole,int *id_vec);
  void measureWidth3D(int N, double arr[], int nslice,string ofname,
		      int *id_vec,vector<array<double,2>> *axis_pos,
		      vector<double> *axis_width,
		      vector<array<double,2>> *end_pos0,
		      vector<array<double,2>> *end_pos1,
		      string ax0);
  void pruneBondZ(int *id_vec);
  int pruneBondSubZ(int ib, int istk,int dir);
  void readData(string fname, string bd3d_tag0); 
  void registerOneBondZ(int ib, int kb, double rCutSq, int istk);
  void removeFloatingFil(int length_cut,int *id_vec);
  void removeOneBondZ(int ib, int kb, int istk);
  void removeUnZbondedBeads(int *id_vec);  
  void setz(double dz, string origin);
  void smoothenFilAverage(int segLength,double rCut,int *id_vec);
  void smoothenZ(int segLength,int delZ, double rCut,int *id_vec);
  void trackFilSeq(int istk, int b0, int b1, vector<int> &list);
  void writeCor(string fname, string outMode, int *id_vec);
  void writeData(string fname, string bd_name,int *id_vec); 
  void writeImage(int beadL, int bondL, int nColor, double *col_bg,
		  double *col_bead, double *col_bond,
		  string fname, string fmt, int *id_vec);
  void writePsf(string fname, int *id_vec);  
};

/************************************************************************/
// Non-member function

void checkDuplicateFnet3dTag(vector<fnet3d> * vfnet3d, string tag);
vector<fnet3d>::iterator findFnet3d(vector<fnet3d> * vfnet3d, string tag);
void getFnet3dIdRange(stringstream& ss, int *ivec,
		      vector<fnet3d>::iterator it_fnet3d);

int fnet3d_cmd(stringstream& ss, ifstream *fin0);

/************************************************************************/
// Global variable
EXT_NAME  vector<fnet3d> vfnet3d;

#endif // #define FNET3D
