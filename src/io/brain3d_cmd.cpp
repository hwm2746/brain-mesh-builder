#include "dodri.h"

int brain3d_cmd(stringstream& ss, ifstream *fin0) {
  int bflag,ivec[2],ivec1[2],ivec2[2];
  double avec[10]; 
  string s_action,tag,tag1,ofname,sdum;
  
  vector<brain3d>::iterator it_brain3d;
  vector<fnet3d>::iterator it_fnet3d;    
  
  tag=getCommOptString(ss,"tag","");
  if (tag=="") error_dodri("Missing tag.");
  
  bflag=0; //check for do block
  if (checkCommOptKey(ss,"do",0)==0) goto BRAIN3D0;
  checkCommand(ss);
  
  bflag=1;
  while (bflag==1) {
    sdum="";
    processCommand(ss,*fin0);
    if (ss.str().size()>0)  cout<<"  "<<tag<<": "<< ss.str() <<endl;
    else continue;
    
  BRAIN3D0: it_brain3d=findBrain3d(&vbrain3d,tag); // find brain3d by input tag
    while (ss.good()) {
      s_action=getCommNext(ss); if (s_action.length()==0) continue;
      
      //check if tag exists
      if ((s_action!="stop")&&(s_action!="done")&&(s_action!="build")) {
	if (it_brain3d==vbrain3d.end())
	  error_dodri("Cannot find brain3d "+ tag);
      }
      
      if (s_action=="stop") {
	cout<< "dodri:: Program finished."<<endl; exit(-1);
      }
      else if (s_action=="done") {
	bflag=0; checkCommand(ss); break;
      }
      
      /************ Constructor ************/

      else if (s_action=="build") {
	checkDuplicateBrain3dTag(&vbrain3d,tag);
	sdum=getCommOptString(ss,"from","");
	tag1=getCommOptString(ss,"model_tag","");	
	if (sdum=="fnet3d") {
	  if (tag1=="") error_dodri("brain3d_cmd: No fnet3d specified.");
	  //it_fnet3d=findFnet3d(&vfnet3d,tag1); // find fnet3d with tag
	  //if (it_fnet3d==vfnet3d.end()) 
	  //  error_dodri("findFnet3d: cannot find fnet3d "+ tag1);
	  //getFnet3dIdRange(ss,ivec,it_fnet3d);
	  //vbrain3d.push_back(brain3d(&(*it_fnet3d),ivec,tag)); 
	} //	else if (sdum=="fnet3d_tag") {
	else vbrain3d.push_back(brain3d(tag)); //blank object
	checkCommand(ss); break;
      } //      else if (s_action=="build") {
      
      /************ Member functions ************/
      
      else if (s_action=="anal_zf_vent") {
	tag1 = getCommOptString(ss,"fnet3d_tag","");
	if (tag1=="") error_dodri("brain3d_cmd: no fnet3d tag specified"); 
	it_fnet3d=findFnet3d(&vfnet3d,tag1); 
	if (it_fnet3d==vfnet3d.end())
	  error_dodri("findFnet3d: cannot find fnet3d "+tag1);
	getFnet3dIdRange(ss,ivec,it_fnet3d);
	getCommOptDouble(ss,"pt",avec,9,-1); // 3 ref pts for plane, xyz
	getCommOptInt(ss, "nslice",ivec1, 1, 10);
	getCommOptInt(ss, "dz",ivec2,1,2);
	sdum=getCommOptString(ss,"axis","y"); // axis for auto-pt finding 
	it_brain3d->analyzeZfVent(&(*it_fnet3d),avec,ivec1[0],
				  ivec2[0],ivec,sdum); 
	checkCommand(ss); break;	
      }

      else if (s_action=="write") {
	ofname=getCommOptString(ss,"name",tag); // output filename
	sdum=getCommOptString(ss,"type",""); // output type
	if (sdum=="") error_dodri("No output type given");
	else if (sdum=="stl")	  {
	  tag1 = getCommOptString(ss,"ref_surf","");
	  if (tag1=="")
	    error_dodri("brain3d_cmd: no reference surface specified"); 
	  it_brain3d->writeSTL(tag1,ofname);
	}
	else error_dodri("Unrecognized output type: "+sdum);
	checkCommand(ss);break;
      } //      else if (s_action=="write") {
      
      /************************/
      
      else { // no proper action keyword
	cout <<"WARNING: Unrecognized action keyword: "<<s_action<<endl;
	checkCommand(ss); break;
      }
    } //while (ss.good())
    ss.clear();
    cout<<endl;
  } //while (bflag==1)
  return 0;
}
