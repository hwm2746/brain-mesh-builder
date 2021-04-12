#ifndef BRAIN3D
#define BRAIN3D

#ifdef SET_EXT
#define EXT_NAME
#else
#define EXT_NAME extern
#endif

/************************************************************************/
class brain3d {
 public:
  string tag;

  // Constructors
  brain3d() {};
  brain3d(string tag0);

  //Member functions
  void analyzeZfVent(fnet3d *fn3d,double arr[],int nslice,
		     int dz,int *id_vec,string ax0);
  void writeSTL(string surf_tag,string fname);   
};
/************************************************************************/
//Non-member functions

void checkDuplicateBrain3dTag(vector<brain3d> * vbrain3d, string tag);
vector<brain3d>::iterator findBrain3d(vector<brain3d> *vbrain3d,string tag);

int brain3d_cmd(stringstream&ss, ifstream *fin0);
/************************************************************************/
// Global variables
EXT_NAME vector<brain3d> vbrain3d;

#endif //#define BRAIN3D
  
