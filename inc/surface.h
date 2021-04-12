#ifndef SURFACE
#define SURFACE

#ifdef SET_EXT
  #define EXT_NAME 
#else
  #define EXT_NAME extern
#endif

#include "fnet3d.h"

/************************************************************************/
class surface: public fnet3d { 
 public:
  int nQuad,fn3d_id_vec[2];
  string tag,fnet3d_tag;
   multimap<int,int> *bead_quad; //quad indexes for a given bead
  vector<array<int,5>> *quad; //stores beads (1-4) of quad (0)

  //surface params
  multimap<int,array<double,3>> quad_cm;
  multimap<int,array<double,3>> quad_normal; 
  
  // Constructors
  surface() {};
  surface(fnet3d *fn3d, int *id_vec, string tag0) ;

  // Member functions
  void addOneQuad(int istk, array<int,5> quadbeads);
  void calcParam(int flag);   
  int checkDuplicateQuad(int istk,array<int,5> quadbeads);
  void drawSurf(string fname,string color,string outMode,int *id_vec);
  void drawVMD(FILE *ff, double d0[3], double d1[3], double a);
  void quadNormOrie(); 
  void registerHexagon();
  void registerMissingSurf();   
  void registerPent(); 
  void registerQuads();
};

/************************************************************************/
// Non-member function

void checkDuplicateSurfaceTag(vector<surface> * vsurface, string tag);
vector<surface>::iterator findSurface(vector<surface> *vsurface, string tag);

int surface_cmd(stringstream& ss, ifstream *fin0);

/************************************************************************/
// Global variable
EXT_NAME  vector<surface> vsurface;

#endif // #define SURFACE
