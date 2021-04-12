#ifndef IMG3D
#define IMG3D

#ifdef SET_EXT
  #define EXT_NAME 
#else
  #define EXT_NAME extern
#endif

#include "img.h"

/************************************************************************/
class img3d{ // Stores raw & grayscale pixel data
 private:
 public: 
  int nStack;
  string tag;
  vector<img> img_stack;

  // constructor
  img3d();
  img3d(string from_tag, string to_tag, string tag0);
  img3d(string ifname, string tag0);
  img3d(img3d *img3d0, string tag1, int *id_vec); //copy constructor

  // Member functions
  void fillRegion2D(int *iniPxl, int *finPxl,int wSize,int foreground,
		    int background,string mode,int *id_vec);  
  void gaussBlur2D(double *s,int *w, unsigned char pxlCut,int *id_vec);
  int imgTag2Id(string sdum);
  void invert2D(int *id_vec);
  void merge(string *merge_tag,string mode,int *id_vec);   
  void removeDebris2D(int sizeCut, unsigned char *pxlCut,
		      unsigned char fillPxl,string mode,int *id_vec);
  void shift(int *dr, unsigned char marginPxl, int *id_vec);
  void threshold(string kind,string mode, double *cutoff, int flag,
		 int *id_vec);  
  void thresholdAuto2D(double cutoff,int nbin,int *id_vec);
  void thresholdHisto2D(double *cutoff, int nbin, int *id_vec);    
  
  void writeEmap(double *dr,string ofname,string outMode,int *id_vec);
  void writeImage(string fname,string fmt,int nColor,
		  double *colorBg, double *colorScale,int *id_vec);
}; 

/************************************************************************/
// Non-member function

void checkDuplicateImg3dTag(string tag);
void getImg3dIdRange(stringstream& ss, int *id_vec, 
		     vector<img3d>::iterator it_img3d);

vector<img3d>::iterator findImg3d(string tag);
int img3d_cmd(stringstream& ss, ifstream *fin0);
/************************************************************************/
// Global variable
EXT_NAME  vector<img3d> vimg3d;

#endif // #define IMG3D
