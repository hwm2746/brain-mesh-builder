#ifndef IO_DODRI
#define IO_DODRI

// I/O-related
extern void checkCommand(stringstream& command);
extern int checkCommOptKey(stringstream& command, string key, int iflag);
extern int checkImageNameExt(string fname);
extern string getCommNext(stringstream& command);
extern void getCommOptDouble(stringstream& command, string commOpt, 
			     double avec[], int nValue, double defaultVal);
extern void getRangeOptDouble(stringstream& command, string commOpt, 
			      double avec[], int nLayers, double defaultValue);
extern void getCommOptInt(stringstream& command, string commOpt, 
			  int ivec[], int nValue, int defaultVal);
extern void getRangeOptInt(stringstream& command, string commOpt, 
			   int ivec[], int nLayers, int defaultValue);
extern string getCommOptString(stringstream& command, string commOpt, 
			       string defaultVal);
extern void getCommOptUchar(stringstream& command, string commOpt, 
			    unsigned char pvec[], int nValue,
			    unsigned char defaultVal);
extern void processCommand(stringstream& ss, ifstream& fin0);
extern int readImageStack(vector<img>* idata, string fmt, string ifname, string tag);
extern void writeImageBead(list<Magick::Drawable> &drawBead, double **cm,
			   int nBead, int beadL, int nColor, double *colorBead);
extern void writeImageBond(list<Magick::Drawable> &drawBond, double **cm,
			   int nBead,int bondL, multimap<int,int> &bond,
			   int nColor, double *colorBond);

#endif //#ifndef IO_DODRI
