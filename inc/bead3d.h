#ifndef BEAD3D
#define BEAD3D

#ifdef SET_EXT
  #define EXT_NAME 
#else
  #define EXT_NAME extern
#endif

#include "bead.h"
#include "img3d.h"

/************************************************************************/
class bead3d{
 private:
 public: 
  int nStack;
  string tag, img3d_tag;
  vector<bead> bead_stack;
  
  // Constructor
  bead3d() {};
  bead3d(img3d *img3d0, string tag0, int iDiam,
	 string tag1, int dir);// area-based
  bead3d(img3d *img3d0,string tag0,int iDiam,
	 string tag1,int verbosity,int *id_vec); //multi-dir area
  bead3d(string tag, string ifname, string image3d_tag0); // from data file  

  // Member functions
  int beadTag2Id(string sdum);
  void orient(double tol, string ax0,int iorie0, int *id_vec);
  void readData(string fname,string image3d_tag0);
  void rotate(double theta,int *id_vec);   
  void setz(double dz, string origin);
  void writeCor(string fname, string outMode, int *id_vec);
  void writeData(string fname, int *id_vec);
  void writeImage(int beadL,int nColor,double *col_bg,double *col_bead,
		  string fname, string fmt,int *id_vec);
  void writePsf(string fname, int *id_vec);  
};  

/************************************************************************/
// Non-member function
void checkDuplicateBead3dTag(vector<bead3d> * vbead3d, string tag);
vector<bead3d>::iterator findBead3d(vector<bead3d> *vbead3d,string tag);
void getBead3dIdRange(stringstream& ss, int *ivec,
		      vector<bead3d>::iterator it_bead3d);
int bead3d_cmd(stringstream& ss, ifstream *fin0);

/************************************************************************/
// Global variable
EXT_NAME  vector<bead3d> vbead3d;

#endif // #define BEAD3D
