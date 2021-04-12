#ifndef IMG
#define IMG

#ifdef SET_EXT
  #define EXT_NAME 
#else
  #define EXT_NAME extern
#endif

/************************************************************************/
class img { // Stores raw & grayscale pixel data
private:
public: 
  string image_name,tag; // source image
  
  // image i/o
  int rows, columns, col_tot; // num rows, columns, columns*ncolor
  unsigned int BitsPerSample;     // normally 8 for grayscale image
  unsigned int SamplesPerPixel;      // normally 1 for grayscale image
  unsigned int TotalColors; // <255 for grayscale 
  
  unsigned char **pxl_raw, **pxl_gray;
  char **pxl_flag; // for marking purpose 
  double dRow,dCol;
  // double avg,sd; 
  int ncolors,range; 
  
  // Constructors
  img () {}; // constructor for declaration
  img(string ifname, string tag0);
  img(const img &img0, unsigned char **pxl, string tag0); //copy constructor
  
  // Member functions
  void fillRegion(int *refPxl,int wSize,int foreground,
		  int background,string mode,string tag1);
  int findCluster(unsigned char pxlCut,
		  vector<vector<pair<int,int>>> &cluster);  
  void gaussBlur(double *s,int *w, unsigned char pxlCut,string tag1);
  void getPxlAreaFraction(unsigned char *pvec, double locut, double hicut);  
  void getStatistics(int loCut, int hiCut, double *avg, double *sd);
  void invert(string tag1);
  void removeDebris(int sizeCut, unsigned char *pxlCut,
		    unsigned char fillPxl,string mode,string tag1);   
  void shift(int *dr, unsigned char marginPxl, string tag1,string tag_ref);
  void threshold(string kind,string mode,double *cutoff,int flag, string tag1);
  void thresholdAuto(double cutoff,int nbin, string tag1);
  void thresholdHisto(double *cutoff, int nbin, string tag1);  
  void writeImage(string fname, string fmt, int nColor,
		  double *colorBg,double *colorScale);
};

/************************************************************************/
// Non-member function

void checkDuplicateImgTag(vector<img> * vimg,string tag);
void copyImg(const img &img0, unsigned char **pxl, string tag1);
vector<img>::iterator findImg(string tag);
void getImgTagRange(stringstream& ss, string *svec);
void mergeImg(img *img0, img *img1, string mode,string tag2);

int img_cmd(stringstream& ss, ifstream *fin0);

/************************************************************************/
// Global variable

EXT_NAME  vector<img> vimg;

#endif // #define IMG
