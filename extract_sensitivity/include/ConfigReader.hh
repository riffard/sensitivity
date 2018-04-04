#ifndef ConfigReader_hh
#define ConfigReader_hh 1

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

enum ScaningMode{Grid, Dichotomy};

template<typename T>
T GetVal(string line){
  string arg;
  T val;
  stringstream ss(line);
  ss >> arg >> val;
  return val;
  }

template<>
ScaningMode GetVal<ScaningMode> (string line) {

  if(line.find("Dichotomy") != string::npos)return Dichotomy;
  else if(line.find("Grid") != string::npos)return Grid;
  else{
    cout<<"Error, the command "<<line<<" is not recognized"<<endl;
    exit(1);
  }  
}
  
class ConfigReader{


public:
  ConfigReader(string filename){
    
    ifstream configFile(filename.c_str());

    string line;
    while(getline(configFile, line)){

      if(!FormatLine(line)) continue;
      
      if(line.find("TotalExpousre") != string::npos){ TotalExpousre =  GetVal<double>(line); }
      else if(line.find("UseNuisanceParam") != string::npos){ UseNuisanceParam = GetVal<bool>(line); }
      else if(line.find("N_pseudo_exp") != string::npos){ N_pseudo_exp = GetVal<int>(line); }
      else if(line.find("vervose_level") != string::npos){  vervose_level= GetVal<int>(line); }
      else if(line.find("p_value_target") != string::npos){  p_value_target= GetVal<double>(line); }
      else if(line.find("p_value_tolerence") != string::npos){ p_value_tolerence = GetVal<double>(line); }
      else if(line.find("scaningMode") != string::npos){  scaningMode = GetVal<ScaningMode>(line); }
      else if(line.find("grid_cross_section_end") != string::npos){  grid_cross_section_end= GetVal<double>(line); }
      else if(line.find("grid_Nsteps") != string::npos){  grid_Nsteps = GetVal<int>(line); }
      else if(line.find("dico_Nsteps_max") != string::npos){  dico_Nsteps_max = GetVal<int>(line); }
      
    }
    
  }
  ~ConfigReader(){}
  


  double Get_TotalExpousre(){return TotalExpousre;}
  bool Get_UseNuisanceParam(){return UseNuisanceParam;}
  int Get_N_pseudo_exp(){return N_pseudo_exp;}
  int Get_vervose_level(){return vervose_level;}
  double Get_p_value_target(){return p_value_target;}
  double Get_p_value_tolerence(){return p_value_tolerence;}
  ScaningMode Get_scaningMode(){return scaningMode; }
  double Get_grid_cross_section_end(){return grid_cross_section_end;}
  int Get_grid_Nsteps(){return grid_Nsteps;}
  int Get_dico_Nsteps_max(){return dico_Nsteps_max;}
  
  

  
private:

  static bool FormatLine(string &line){

    while(line.front() == ' ') line.erase(line.begin(), line.begin()+1);

    if(line.size() == 0) return false;
    if(line.front() == '#') return false;
    
    if(line.find("#") != string::npos)line.erase( line.begin()+line.find("#"), line.end());
    EraseChar(line, "=");
    EraseChar(line, ":");
    EraseChar(line, ",");
    EraseChar(line, ":");
    
    return true;
  }

  static void EraseChar(string &line, string to_erase){
    
    size_t found = line.find(to_erase);

    while(found != string::npos){
      line.erase(line.begin()+found, line.begin()+found+1);
      found = line.find(to_erase);
    }
    
  }

  
  double TotalExpousre;
  bool UseNuisanceParam;
  int N_pseudo_exp;
  int vervose_level;
  double p_value_target; 
  double p_value_tolerence;
  ScaningMode scaningMode;
  double grid_cross_section_end;
  int grid_Nsteps;
  int dico_Nsteps_max ;
  
  
 


};


#endif
