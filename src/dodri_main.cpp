/* dodri_main.cpp : Computer-Aided Feature Extraction

(c) Wonmuk Hwang 2013-2021, hwm@tamu.edu

Contributor: Ana Chang-Gonzalez

This program is copyrighted under the terms of the GNU General Public License
version 3. See COPYRIGHT.TXT that came along with this program.

The copyright includes all source files and documentation.

* HOW TO CITE:

Title: Building 3-Dimensional Model of Early-Stage Zebrafish Embryo Midbrain-Hindbrain Boundary

Authors: Ana C. Chang-Gonzalez, Holly C. Gibbs, Arne C. Lekven, Alvin T. Yeh, and Wonmuk Hwang

TO BE PUBLISHED (2021).


*/

#define SET_EXT 
#include "dodri.h"

int main(int argc, char *argv[])
{
  if (argc!=2) { 
    cout<<"Include argument in the command line: ./dodri input.dat" 
	<< endl; return 1;
  }

  // s_comm: for command line, s_action: command action
  string sdum,sdum1,sdum2,s_comm,s_action; stringstream ss;
  string ifname="",tag,tag1,outTag="";

  ifstream fin0(argv[1]); // Open input data file
  if (!fin0) error_dodri("img_cmd: Input file does not exist.");
  
  try {
    InitializeMagick(*argv);
  }
  catch (exception &error_) {
    cout << "Caught exception:" << error_.what()<<endl;
    return 1;
  }
  
  /*******************************************************
   //  Process input
   *******************************************************/
  while (!fin0.eof()) {
    
    sdum=""; //reset sdum
    processCommand(ss, fin0);
    if (ss.str().size() >0)  cout<<"CMD: "<< ss.str() <<endl;
    else continue;
    
    while (ss.good()) {
      sdum=getCommNext(ss); if (sdum.length()==0) continue;
      if (sdum=="stop") {
	cout<< endl<<"DODRI:: Program finished."<<endl; exit(-1);
      }
      
      /************ image-related ************/
      else if (sdum=="img") {
	img_cmd(ss,&fin0); break;
      }
      /************ img3d-related ************/
      else if (sdum=="img3d") {
      	img3d_cmd(ss,&fin0); break;
      }
      /************ bead-related ************/
      else if (sdum=="bead") {
      	bead_cmd(ss,&fin0); break;
      }
      /************ bead3d-related ************/
      else if (sdum=="bead3d") {
      	bead3d_cmd(ss,&fin0); break;
      }
      /************ fnet-related ************/
      else if (sdum=="fnet") {
      	fnet_cmd(ss,&fin0); break;
      }
      /************ fnet3d-related ************/
      else if (sdum=="fnet3d") {
      	fnet3d_cmd(ss,&fin0); break;
      }
      /************ surface-related ************/
      else if (sdum=="surface") {
      	surface_cmd(ss,&fin0); break;
      }
      ///************ brain analysis-related ************/
      //else if (sdum=="brain") {
      //	brain_cmd(ss,&fin0); break;
      //}
      /************ brain3d analysis-related ************/
      else if (sdum=="brain3d") {
      	brain3d_cmd(ss,&fin0); break;
      }
      
      /************ Output-related ************/
      else if (sdum=="timer_start")  {
	timer=clock();
	cout<<"Starting timer."<<endl; break;
      }
      else if (sdum=="timer_end") {
	//t_end=time(0);
	//difftime(t_end,t_start); //real time
	timer=clock()-timer; timer=(float)timer/CLOCKS_PER_SEC; 
	cout<<"Run time: "<<timer<<" seconds."<<endl;
	break; 
      }

      /******************************/
      else {
	sdum1="Unrecognized token: "+sdum;
	error_dodri(sdum1);
      }
    } //     while (ss.good()) {
    ss.clear(); // clear state flag
  }
  return 0;
}
