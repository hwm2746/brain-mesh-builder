#include "dodri.h"

/************************************************************************/
void checkCommand(stringstream& command)
{
  string sdum=command.str();
  ltrim(sdum);
  if (sdum!="") cout<<"  WARNING: extraneous keywords: "<< sdum<<endl;
  return;
}

/************************************************************************/
int checkCommOptKey(stringstream& command, string key, int iflag)
{ /* Check if command string contains a key.
     returns: 0 - key not found , 1 - found
     iflag=0: remove the key on return. !=0: do not remove
 */
  string sdum, sdum1;
  int flag=0;
  stringstream ss, ss1;
  
  getline(command,sdum1); // backup
  ss1<<sdum1;
  while (!ss1.eof()) {
    ss1 >> sdum;
    if (sdum==key) flag=1; // 
    else ss << sdum << " "; // store non-option params
    sdum="";
  }
  command.clear(); // restore
  if (iflag==0) command.str(ss.str()); // return
  else command.str(sdum1); // keep the original
  return flag;
}

/************************************************************************/
string getCommNext(stringstream& command)
{ /* Output the next keyword from command and delete it. */
  
  stringstream ss;
  string sdum,sdum1;

  command >> sdum1;
  while (!command.eof()) {
    command >> sdum;    ss << sdum << " "; // store non-option params
    sdum=""; // reset
  }
  command.clear(); // restore
  sdum=ss.str();
  command.str(sdum);
  return sdum1;
}
/************************************************************************/
void getCommOptInt(stringstream& command, string commOpt, 
		   int ivec[], int nValue, int defaultVal)
{ /* From command, read nValue of integer numbers after keyword commOpt.
     Returns: array of values read. */
  
  int i,flag=0;
  string sdum,sdum1, arg_s[nValue];
  stringstream ss;
  
  getline(command,sdum1);  command.clear(); command.str(sdum1);// backup
  while (!command.eof()) {
    command >> sdum;
    if (sdum==commOpt) { 
      for (i=0;i<nValue; i++) {
	if (!command.fail()) { 
	  command >> ivec[i];
	  if (command.fail()) { 
	    i--;
	    continue;
	  }
	}
	else ivec[i] = defaultVal;
      }
      ++flag;   //  break;
    }
    else ss << sdum << " ";
    sdum="";
  }
  
  if (flag==0) for (i=0;i<nValue; i++) ivec[i]=defaultVal; 
  command.clear();  // restore
  command.str(ss.str());
}

/************************************************************************/
void getCommOptDouble(stringstream& command, string commOpt, 
		      double avec[], int nValue, double defaultVal)
{ /* From command, read nValue of double numbers after keyword commOpt.
     Returns: array of values read.  */
  
  int i,flag=0;  
  string sdum,sdum1;
  stringstream ss;
  
  getline(command,sdum1);  command.clear(); command.str(sdum1);// backup
  while (!command.eof()) {
    command >> sdum;
    if (sdum==commOpt) { 
      for (i=0;i<nValue; i++) {
	if (!command.fail()) {
	  command >> avec[i];
	  if (command.fail()) { 
	    i--;
	    continue;
	  }
	}
	else avec[i] = defaultVal;
      }
      ++flag;   // break;
    }
    else ss << sdum << " ";
    sdum="";
  }
  if (flag==0) for (i=0;i<nValue; i++) avec[i]=defaultVal; 
  command.clear(); // restore
  command.str(ss.str());
}

/************************************************************************/
string getCommOptString(stringstream& command, string commOpt, 
			string defaultVal)
{ /* From command string read string command option value for commOpt. */
  string sdum, r, sdum1,sdum2; int flag=0;
  stringstream ss;
  
  getline(command,sdum1);  command.clear(); command.str(sdum1);// backup
  //command.clear(command.goodbit); 
  while (!command.eof()) {
    command >> sdum;
    if (sdum==commOpt) { command >> r; ++flag; } // break;}
    else ss << sdum << " ";
    sdum="";
  }
  sdum2=(flag==0)?defaultVal:r;
  command.clear();   // restore
  command.str(ss.str());  
  return sdum2;
}

/************************************************************************/
void getCommOptUchar(stringstream& command, string commOpt, 
		     unsigned char pvec[], int nValue, unsigned char defaultVal)
{ /* From command, read nValue of unsigned characters after keyword commOpt.
     Returns: array of values read.
  */
  
  int i,flag=0, idum;
  string sdum,sdum1, arg_s[nValue];
  stringstream ss;
  
  getline(command,sdum1);  command.clear(); command.str(sdum1);// backup
  while (!command.eof()) {
    command >> sdum;
    if (sdum==commOpt) { 
      for (i=0;i<nValue; i++) {
	if (!command.fail()) {
	  // read as integer and convert to unsigned char, since
	  // command >> pvec[i] only reads in a single character.
	  command >> idum; 
	  pvec[i]=(unsigned char)idum;
	  if (command.fail()) { 
	    i--;
	    continue;
	  }
	}
	else pvec[i] = defaultVal;
      }
      ++flag;   //  break;
    }
    else ss << sdum << " ";
    sdum="";
  }
  
  if (flag==0) for (i=0;i<nValue; i++) pvec[i]=defaultVal; 
  command.clear();  //command.str(sdum1); // restore
  command.str(ss.str());
}

/************************************************************************/
void processCommand(stringstream& ss, ifstream& fin0)
{  /*  Function to process commands from input files. */
  
  string sdum;
  string s_comm="";
  int flag = 0;
  do {
    getline(fin0,sdum);
    //  trim leading spaces & convert to lowercase
    ltrim(sdum); 
    std::transform(sdum.begin(), sdum.end(), sdum.begin(), ::tolower);
    if (flag && sdum.substr(0, sdum.find("#")) == "") continue;
    // remove comment & append to prior s_comm
    s_comm.append(sdum.substr(0, sdum.find("#")));
    rtrim(s_comm); // trim trailing spaces
    if (*(s_comm.end()-1) == '\\') { // check for append line operator
      flag = 1;
      s_comm=s_comm.substr(0, s_comm.size()-1); // remove the last '\'
      rtrim(s_comm);
      s_comm.append(" ");
    }
    else {
      flag = 0;
    }
  } while (flag);
  ss.str(""); ss.clear(); ss<< s_comm;
}

/************************************************************************/
int readImageStack(vector<img>* idata, string ifname, string fmt,
		   string tag)
{ /* Read image stack listed in ifname.
     fmt: unused for now.
     Returns: number of images read (nstack)  */
  
  int i;
  img idata0;
  ifstream fin0;
  stringstream ss;
  string sdum,tag0;

  fin0.open(ifname.c_str());
  if (!fin0.is_open()) error_dodri("readImage:: Failed to open file: "+ifname);

  i=0;
  while (!fin0.eof()) {
    getline(fin0,sdum);
    removeComment(sdum);
    ltrim(sdum); // trim leading spaces
    if (sdum!="")  cout<<endl; // add a blank line
    ss.str(""); ss<< sdum;
    while (ss.good()) {
      ss >> sdum; if (sdum.length()==0) continue;
      tag0=tag+itostr(i++);
      checkDuplicateImgTag(idata,tag0);
      idata0=img(sdum,tag0);  
      (*idata).push_back(idata0);
      break; 
    } // while (ss.good()) {
    ss.clear(); // clear state flag
  } // while (!fin0.eof()) {
  cout<<endl<<"  Read "<< i<<" image files from "
      <<ifname<<endl;
  fin0.close();
  return i;
}

/************************************************************************/
void writeImageBead(list<Magick::Drawable> &drawBead,  double **cm,
		    int nBead, int beadL, int nColor, double *colorBead)
{ /* Add bead rendering to the list drawBead */
  int iBead;
  double x0,y0, r0=0.5*(double)beadL;
  
  drawBead.push_back(DrawableStrokeAntialias(false));
  
  if (nColor==3) {
    drawBead.push_back(DrawableStrokeColor(ColorRGB(colorBead[0],colorBead[1],
						    colorBead[2])));
    drawBead.push_back(DrawableFillColor(ColorRGB(colorBead[0],colorBead[1],
						  colorBead[2])));
  }
  else {
    drawBead.push_back(DrawableFillColor(ColorGray(*colorBead)));
  }

  for (iBead=0;iBead<nBead;iBead++) {
    x0=cm[iBead][1];  y0=cm[iBead][0];    
    drawBead.push_back(DrawableCircle(x0,y0,x0+r0,y0));
  }
  return;
}
  
/************************************************************************/
void writeImageBond(list<Magick::Drawable> &drawBond, double **cm,
		    int nBead,int bondL, multimap<int,int> &bond,
		    int nColor, double *colorBond)
{ /* Add bond rendering to the list drawBond */

  int iBead,jBead;
  double x0,y0,x1,y1;

  multimap<int,int>::iterator it;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> it_range;
  
  drawBond.push_back(DrawableStrokeAntialias(false));
  if (nColor==3) {
    drawBond.push_back(DrawableStrokeColor(ColorRGB(colorBond[0],colorBond[1],
						    colorBond[2])));
    drawBond.push_back(DrawableFillColor(ColorRGB(colorBond[0],colorBond[1],
						  colorBond[2])));
    
  }
  else {
    drawBond.push_back(DrawableStrokeColor(ColorGray(*colorBond)));
    drawBond.push_back(DrawableFillColor(ColorGray(*colorBond)));
  }
  drawBond.push_back(DrawableStrokeWidth(bondL)); 

  for (iBead=0;iBead<nBead;iBead++) {
    x0=cm[iBead][1];  y0=cm[iBead][0];
    it_range=bond.equal_range(iBead);     
    for (it=it_range.first;it!=it_range.second;it++) {
      jBead=it->second;
      if (jBead>=iBead) continue; // avoid double rendering
      x1=cm[jBead][1];  y1=cm[jBead][0];
      drawBond.push_back(DrawableLine(x0,y0,x1,y1));
    }
  } // for (iBead=0;iBead<nBead;iBead++) {


  return;
}
