#ifndef _FILEREADING_H_
#define _FILEREADING_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <map>
#include "Vector.h"
#include <boost/algorithm/string.hpp>
using namespace std;


bool ReadDouble(std::ifstream &Filestream, double &Out)
{
    std::string DataLine;
    if(getline(Filestream, DataLine, '\n'))
    {
        stringstream linestream(DataLine);
        linestream >> Out;
        return true;
    }
    else
    {
        return false;
    }
}

bool ReadVector(std::ifstream &Filestream, Vector &Out)
{
    std::string DataLine;
    if(getline(Filestream, DataLine, '\n'))
    {
        stringstream linestream(DataLine);
        linestream >> Out.x >> Out.y >> Out.z;
        return true;
    }
    else
    {
        return false;
    }
}


bool ReadString(std::ifstream &Filestream, std::string &Out)
{
    std::string DataLine;
    if(getline(Filestream, DataLine, '\n'))
    {
        stringstream linestream(DataLine);
        linestream >> Out;
        return true;
    }
    else
    {
        return false;
    }
}


void AddToMapFromFile(std::ifstream &Filestream, std::map<std::string,std::string> &MapOut)
{
    std::string DataLine;
    while(getline(Filestream,DataLine,'\n'))
    {
        stringstream linestream(DataLine);
        std::string Key,Val;
        linestream >> Key >> Val;
        MapOut.insert(std::pair<std::string,std::string>(Key,Val));
    }
}

bool VectorFromString(std::string S, Vector &V)
{
    std::vector<std::string> strs;
    boost::split(strs, S, boost::is_any_of(std::string("\t ,")));
    if(strs.size() != 3)
    {
      return false;
    }
    else
    {
      V.x = atof(strs[0].c_str());
      V.y = atof(strs[1].c_str());
      V.z = atof(strs[2].c_str());
      return true;
    }
}

void DoubleFromString(std::string S, double &x)
{
    x = atof(S.c_str());
}

void IntFromString(std::string S, int &x)
{
    x = atoi(S.c_str());
}


bool VectorFromMap(std::string key, std::map<std::string,std::string> &Map, Vector &V)
{
    if( Map.count(key) == 0)
    {
        cout << "Warning: Unable to find key (" << key << ") in Map" << endl;
    }    
    
    std::string Val = Map.find(key)->second;
    return VectorFromString(Val, V);
}

void DoubleFromMap(std::string key, std::map<std::string,std::string> &Map, double &x)
{
    if( Map.count(key) == 0)
    {
        cout << "Warning: Unable to find key (" << key << ") in Map" << endl;
    }
    
    std::string Val = Map.find(key)->second;
    DoubleFromString(Val, x);
}

void StringFromMap(std::string key, std::map<std::string,std::string> &Map, std::string &S)
{
    if( Map.count(key) == 0)
    {
        cout << "Warning: Unable to find key (" << key << ") in Map" << endl;
    }
    
    S = Map.find(key)->second;
}

void IntFromMap(std::string key, std::map<std::string,std::string> &Map, int &i)
{
    if( Map.count(key) == 0)
    {
        cout << "Warning: Unable to find key (" << key << ") in Map" << endl;
    }
    
    std::string Val = Map.find(key)->second;
    IntFromString(Val,i);
}




#endif
