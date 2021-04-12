#include "dodri.h"

int surface_cmd(stringstream& ss, ifstream *fin0) {
  int bflag,ivec[2]; 
  string s_action, tag,tag1,ofname,sdum1,sdum2; 
  vector<surface>::iterator it_surface;
  vector<fnet3d>::iterator it_fnet3d;  

  tag=getCommOptString(ss,"tag","");
  if (tag=="") error_dodri("Missing tag.");

  bflag=0; //check for do block
  if (checkCommOptKey(ss,"do",0)==0) goto SURFD0; //no do block
  checkCommand(ss);
  
  bflag=1;
  while (bflag==1) {
    processCommand(ss,*fin0);
    if (ss.str().size() >0) cout<<"  "<<tag<<": "<< ss.str() <<endl;
    else continue;
    
  SURFD0: it_surface=findSurface(&vsurface,tag);
    while(ss.good()) {
      s_action=getCommNext(ss); if (s_action.length()==0) continue;
      
      // check if tag exists
      if ((s_action!="stop")&&(s_action!="done")&&(s_action!="build")) {
	if (it_surface==vsurface.end())
	  error_dodri("Cannot find surface "+ tag);
      }
      if (s_action=="stop") {
	cout<< "dodri:: Program finished."<<endl; exit(-1);
      }
      else if (s_action=="done") {
	bflag=0; checkCommand(ss); break;
      }

      /************ Constructor ************/

      else if (s_action=="build") {
	checkDuplicateSurfaceTag(&vsurface,tag);
	tag1=getCommOptString(ss,"fnet3d",""); // tag of fnet3d to use
	it_fnet3d=findFnet3d(&vfnet3d,tag1); // find fnet3d with tag
	if (it_fnet3d==vfnet3d.end()) 
	  error_dodri("findFnet3d: cannot find fnet3d "+ tag1);
	getFnet3dIdRange(ss,ivec,it_fnet3d);
	vsurface.push_back(surface(&(*it_fnet3d),ivec,tag));
	checkCommand(ss);break;
      } //       else if (s_action=="build") {

      /************ Member functions ************/
      
      else if (s_action=="draw_surface") {
	ofname=getCommOptString(ss,"name","ffquad"); 
	sdum1=getCommOptString(ss,"color","green");
	sdum2=getCommOptString(ss,"mode","image"); 
	getCommOptInt(ss,"range",ivec,2,-1); 
	it_surface->drawSurf(ofname,sdum1,sdum2,ivec);
	checkCommand(ss);break;
      } //       else if (s_action=="draw_surface") {

      else if (s_action=="get_centroid") {
	it_surface->calcParam(0);
	checkCommand(ss);break;
      } //       else if (s_action=="get_centroid") {
      
      else if (s_action=="get_normal") {
	it_surface->calcParam(1); 
	checkCommand(ss);break;
      } //       else if (s_action=="get_normal") {

      else if (s_action=="register_hexagon") {
	it_surface->registerHexagon();
	checkCommand(ss);break;
      } //       else if (s_action=="register_hex") {
      
      else if (s_action=="register_missing") {
	it_surface->registerMissingSurf();
	checkCommand(ss);break;
      } //       else if (s_action=="register_missing")      

      else if (s_action=="register_pentagon") {
	it_surface->registerPent();
	checkCommand(ss);break;
      } //       else if (s_action=="register_pentagon")
      
      else if (s_action=="register_quads") {
	it_surface->registerQuads();
	checkCommand(ss);break;
      } //       else if (s_action=="register_quads") {

      /************************/
      
      else { // no proper action keyword
	cout <<"WARNING: Unrecognized action keyword: "<<s_action<<endl;
	checkCommand(ss);break;
      }      
    } //while(ss.good())  
    ss.clear();
    cout<<endl;
  } //while (bflag==1) 
  return 0;
}

