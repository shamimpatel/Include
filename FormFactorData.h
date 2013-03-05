#ifndef _FORMFACTORDATA_H_
#define _FORMFACTORDATA_H_


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
using namespace std;

struct FormFactorDataPoint
{
  float x; //sin(theta/2)/lambda
  float FormFactor; //f(x)
};

bool FormFacLssThnComp( FormFactorDataPoint A, FormFactorDataPoint B)
{
  if(A.x < B.x)
    {
      return true;
    }
  else
    {
      return false;
    }
}

class FormFactorData
{
 private:
  float MinX, MaxX;
  int NumDataPoints;
  float XDelta;

  std::vector<float> FormFactorDataPoints;

 private:

  FormFactorDataPoint FindFormFactorDataPointFromFile( float x, std::vector<FormFactorDataPoint> *FormFacData)
  {
    FormFactorDataPoint DummySearch;
    DummySearch.x = x;

    vector<FormFactorDataPoint>::iterator C = lower_bound( FormFacData->begin(),
							   FormFacData->end(),
							   DummySearch,
							   FormFacLssThnComp );
    if( C == FormFacData->begin() )
      {
	//then we're done
	DummySearch = *C;
	//cout << DummySearch.Energy << endl; //can only get first element if it's a match
      }
    else
      {
	float weight = (DummySearch.x - (*(C-1)).x)/((*C).x - (*(C-1)).x); //requested energy is always >= C-1.energy

	//interpolate between values
	DummySearch.FormFactor = lerp( (*(C-1)).FormFactor, (*C).FormFactor, weight);
	//could also interpolate the last one too
      }

    return DummySearch;
  }

 public:
  FormFactorData(float Min, float Max, int NumDataPoints)
  {
    this->MinX = Min;
    this->MaxX = Max;
    this->NumDataPoints = NumDataPoints;
  }


 public:
  void LoadData( const char* Filename )
  {
    ifstream datafile(Filename);
    if(datafile.is_open() == false)
    {
        cout << "Error: Failed to open FormFactor file: " << Filename << endl;
        exit(1);
    }
    string dataline;

    std::vector<FormFactorDataPoint> FormFactorFileData;

    FormFactorDataPoint DataPoint;

    while(getline(datafile, dataline, '\r')) //excel copy/paste gives \r not \n
    {
      stringstream linestream(dataline);
      linestream >> DataPoint.x >> DataPoint.FormFactor;

      FormFactorFileData.push_back( DataPoint );
    }

    FormFactorDataPoints.clear();

    XDelta = (MaxX-MinX)/float(NumDataPoints);

    for(float i = 0; i <= NumDataPoints; i++)
    {
      float x = MinX + i*XDelta;

      FormFactorDataPoint C = FindFormFactorDataPointFromFile( x, &FormFactorFileData);


      FormFactorDataPoints.push_back( C.FormFactor );
    }



    /*int N = 5000;

    for(int i = 0; i <= N; i++)
    {
      float x = lerp(MinX,MaxX, float(i)/float(N));

      float NearF = GetFormFactorDataPoint(x); //what we get from just finding an f(x) near to x

      //float InterpF = FindFormFactorDataPointFromFile( x, &FormFactorFileData).FormFactor;

      //cout << x << '\t' << NearF << "\t" << InterpF << "\t" << fabs(InterpF - NearF)/InterpF  << endl;
      }*/

  }

  float GetFormFactorDataPoint( float x)
  {
    float fIndex = (x-MinX)/(XDelta);
    //return FormFactorDataPoints[int(fIndex)]; //With a large number of points this is close enough

    //Code below linearly interpolates to find the "best" value.
    //Will crash if x = MaxX. Fix either by checking for this or bumping an extra value onto the end
    //of the FormFactorDataPoints array
    //float weight = fIndex - int(fIndex);
    //return lerp(FormFactorDataPoints[int(fIndex)], FormFactorDataPoints[int(fIndex) + 1], weight);

    if(int(fIndex) >= (int(FormFactorDataPoints.size()) - 1) )
    {
      return FormFactorDataPoints[FormFactorDataPoints.size() - 1]; //just pull out the last one
    }

    if(int(fIndex) < 0 )
    {
      return FormFactorDataPoints[0];
    }

    float weight = fIndex - int(fIndex);



    return lerp(FormFactorDataPoints[int(fIndex)], FormFactorDataPoints[int(fIndex) + 1], weight);


  }


  float lerp( float start, float end, float weight)
  {
    return start + ((end-start)*weight);
  }

};


#endif
