#ifndef _ABSORBCOEFFDATA_H_
#define _ABSORBCOEFFDATA_H_


#include "DataLoader.h"

class AbsorbCoeffData : public DataLoader
{
public:
    
    AbsorbCoeffData( float MinWavelength, float MaxWavelength, int NumDataColumns, int NumDataPoints, const char* Filename):
    DataLoader( MinWavelength, MaxWavelength, NumDataColumns, NumDataPoints, Filename){}
    
    double GetAbsorbCoeffDataPoint( double Wavelength )
    {
        return GetDataPoint(Wavelength).ColumnValues[0];
    }
    
};


#endif



/*
 //This Class loads SORTED (smallest to largest wavelength) data from a file and then creates a fixed number of points within that data between the requested intervals.
 
 
 #include <iostream>
 #include <fstream>
 #include <cstdlib>
 #include <vector>
 #include <string>
 #include <sstream>
 #include <cmath>
 #include <algorithm>
 
 using namespace std;
 
 struct AbsorbCoeffDataPoint
 {
 float Wavelength;
 float Mu;
 };
 
 bool AbsorbCoeffLssThnComp( AbsorbCoeffDataPoint A, AbsorbCoeffDataPoint B)
 {
 if(A.Wavelength < B.Wavelength)
 {
 return true;
 }
 else
 {
 return false;
 }
 }
 
 class AbsorbCoeffData
 {
 private:
 float Min, Max;
 int NumDataPoints;
 float Delta;
 
 std::vector<float> AbsorbCoeffDataPoints;
 
 private:
 
 AbsorbCoeffDataPoint FindAbsorbCoeffDataPointFromFile( float Wavelength, std::vector<AbsorbCoeffDataPoint> *AbCoeffData)
 {
 AbsorbCoeffDataPoint DummySearch;
 DummySearch.Wavelength = Wavelength;
 
 vector<AbsorbCoeffDataPoint>::iterator C = lower_bound( AbCoeffData->begin(),
 AbCoeffData->end(),
 DummySearch,
 AbsorbCoeffLssThnComp );
 if( C == AbCoeffData->begin() )
 {
 //then we're done
 DummySearch = *C;
 //cout << DummySearch.Energy << endl; //can only get first element if it's a match
 }
 else
 {
 float weight = (DummySearch.Wavelength - (*(C-1)).Wavelength)/((*C).Wavelength - (*(C-1)).Wavelength); //requested energy is always >= C-1.energy
 
 //interpolate between values
 DummySearch.Mu = lerp( (*(C-1)).Mu, (*C).Mu, weight);
 //could also interpolate the last one too
 }
 
 return DummySearch;
 }
 
 public:
 
 
 AbsorbCoeffData(float Min, float Max, int NumDataPoints)
 {
 this->Min = Min;
 this->Max = Max;
 this->NumDataPoints = NumDataPoints;
 }
 
 
 public:
 void LoadData( const char* Filename )
 {
 ifstream datafile(Filename);
 if(datafile.is_open() == false)
 {
 cout << "Error: Failed to open Absorption Coefficient file: " << Filename << endl;
 exit(1);
 }
 string dataline;
 
 std::vector<AbsorbCoeffDataPoint> AbsorbCoeffFileData;
 
 AbsorbCoeffDataPoint DataPoint;
 
 while(getline(datafile, dataline, '\r')) //excel copy/paste gives \r not \n
 {
 stringstream linestream(dataline);
 linestream >> DataPoint.Wavelength >> DataPoint.Mu;
 
 AbsorbCoeffFileData.push_back( DataPoint );
 }
 
 AbsorbCoeffDataPoints.clear();
 
 Delta = (Max-Min)/float(NumDataPoints);
 
 for(float i = 0; i <= NumDataPoints; i++)
 {
 float Wavelength = Min + i*Delta;
 
 AbsorbCoeffDataPoint C = FindAbsorbCoeffDataPointFromFile( Wavelength, &AbsorbCoeffFileData);
 
 
 AbsorbCoeffDataPoints.push_back( C.Mu );
 }
 
 datafile.close();
 }
 
 float GetAbsorbCoeffDataPoint( float Wavelength )
 {
 float fIndex = (Wavelength-Min)/(Delta);
 //return AbsorbCoeffDataPoints[int(fIndex)]; //With a large number of points this is close enough
 
 //Code below linearly interpolates to find the "best" value.
 //Will crash if x = MaxX. Fix either by checking for this or bumping an extra value onto the end
 //of the FormFactorDataPoints array
 
 //bounds checking in case we don't have data for particular energies
 if(int(fIndex) >= (int(AbsorbCoeffDataPoints.size()) - 1) )
 {
 return AbsorbCoeffDataPoints[AbsorbCoeffDataPoints.size() - 1]; //just pull out the last one
 }
 
 if(int(fIndex) < 0 )
 {
 return AbsorbCoeffDataPoints[0];
 }
 
 float weight = fIndex - int(fIndex);
 
 
 
 return lerp(AbsorbCoeffDataPoints[int(fIndex)], AbsorbCoeffDataPoints[int(fIndex) + 1], weight);
 
 
 }
 
 
 float lerp( float start, float end, float weight)
 {
 return start + ((end-start)*weight);
 }
 
 };
 */