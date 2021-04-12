#ifndef BEAD
#define BEAD

#ifdef SET_EXT
  #define EXT_NAME 
#else
  #define EXT_NAME extern
#endif

#include "img.h"

/************************************************************************/
class bead { 
 public: 
  string img_tag, tag; // source image name and tag 
  int maxBead; 
  int nBead, beadDiam; 
  double **cm, *beadMass, zCoord;
  // cm: center of mass, beadMass: mass of each bead
  
  int rows,columns;
  double dRow,dCol;
  multimap<int,int> cell_bead; // bead indexes for a given cell number
  multimap<int,int> bead_cell; // cell indexes for a given bead number

  // Constructors
  bead () {}; // constructor for declaration
  bead(const bead &b0,double **cm1,double *beadMass1,string tag1); //copy constructor  
  bead(img *img0, string tag0, int iDiam,
       string tag1, int dir); // pixel area-based
  bead(img *img0, string tag0, int iDiam,string tag1, //for area
       double threshold, double loCut, int flag,  //for center
       int hiCut, int delCut, double wFrac, //for boundary
       int verbosity, string mode); //multi-dir constructor, all modes
  bead(string tag, string ifname, string image_tag0, ifstream *fin = NULL); // from data file
  
  // Member functions
  void createBead(double mass, double cmr, double cmc);
  int decimate(double rCut,int verbosity);
  double getAngleOffset();
  int getBeadDir(double rCut,multimap<int,array<double,2>> &cm_dir,
		 multimap<int,array<double,2>> &cm_buf);  
  void gradient(double rCut,int iRad,
		multimap<int,array<double,2>> &cm_dir);  
  void merge(bead *bd, double rCut, int verbosity);
  int moveBeadPos(int p, double *newPos);
  void readData(string fname, string image_tag0, ifstream *fin = NULL);
  void rotate(double theta, double *pt=NULL); 
  void setScanDir(int iDiam,int dir,int& rStart,int& cStart,int& rEnd,
		  int& cEnd,int& rDelta,int& cDelta);
  int withinRadius(int beadNum, double radius,
		   vector<pair<double,int> >& closeBeads);  
  void writeCor(string fname, string outMode);
  void writeData(string fname, FILE *dat = NULL);
  void writeImage(int beadL, int nColor, double *col_bg, double *col_bead,
		  string fname, string fmt);
  void writePsf(string fname);

};

/************************************************************************/
// Non-member function
void checkDuplicateBeadTag(vector<bead> * vbead,string tag);
void copyBead(const bead &b0, double **cm1, double *beadMass1, string tag1); 
vector<bead>::iterator findBead(vector<bead> * vbead,string tag);
int bead_cmd(stringstream& ss, ifstream *fin0);

/************************************************************************/
// Global variable

EXT_NAME  vector<bead> vbead;

#endif // #define BEAD
