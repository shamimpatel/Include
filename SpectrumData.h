//
//  SpectrumData.h
//  TestCode
//
//  Created by Shamim Patel on 06/03/2013.
//  Copyright (c) 2013 Shamim Patel. All rights reserved.
//

#ifndef _SPECTRUMDATA_H_
#define _SPECTRUMDATA_H_

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>

using namespace std;

struct SpectrumDataPoint
{
    float Energy;
    float RelativeIntensity;
};

bool SpectrumDataLssThnComp( SpectrumDataPoint A, SpectrumDataPoint B)
{
    if(A.Energy < B.Energy)
    {
        return true;
    }
    else
    {
        return false;
    }
}

class SpectrumData
{
private:
    float Min, Max;
    int NumDataPoints;
    float Delta;
    
    std::vector<float> SpectrumDataPoints;
    
private:
    
    SpectrumDataPoint FindSpectrumDataPointFromFile( float Energy, std::vector<SpectrumDataPoint> *SpecData)
    {
        SpectrumDataPoint DummySearch;
        DummySearch.Energy = Energy;
        
        vector<SpectrumDataPoint>::iterator C = lower_bound(   SpecData->begin(),
                                                               SpecData->end(),
                                                               DummySearch,
                                                               SpectrumDataLssThnComp );
        if( C == SpecData->begin() )
        {
            //then we're done
            DummySearch = *C;
            //cout << DummySearch.Energy << endl; //can only get first element if it's a match
        }
        else
        {
            float weight = (DummySearch.Energy - (*(C-1)).Energy)/((*C).Energy - (*(C-1)).Energy); //requested energy is always >= C-1.energy
            
            //interpolate between values
            DummySearch.RelativeIntensity = lerp( (*(C-1)).RelativeIntensity, (*C).RelativeIntensity, weight);
            //could also interpolate the last one too
        }
        
        return DummySearch;
    }
    
public:
    
    
    SpectrumData(float Min, float Max, int NumDataPoints)
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
            cout << "Error: Failed to open Spectrum Data file: " << Filename << endl;
            exit(1);
        }
        string dataline;
        
        std::vector<SpectrumDataPoint> SpectrumFileData;
        
        SpectrumDataPoint DataPoint;
        
        while(getline(datafile, dataline, '\r')) //excel copy/paste gives \r not \n
        {
            stringstream linestream(dataline);
            linestream >> DataPoint.Energy >> DataPoint.RelativeIntensity;
           
            if(DataPoint.RelativeIntensity < 0.0)
            {
                DataPoint.RelativeIntensity = 0.0;
            }
           
            SpectrumFileData.push_back( DataPoint );
        }
        
        SpectrumDataPoints.clear();
        
        Delta = (Max-Min)/float(NumDataPoints);
        
        for(float i = 0; i <= NumDataPoints; i++)
        {
            float Energy = Min + i*Delta;
            
            SpectrumDataPoint C = FindSpectrumDataPointFromFile( Energy, &SpectrumFileData);
            
            
            SpectrumDataPoints.push_back( C.RelativeIntensity );
        }
        
        datafile.close();
    }
    
    float GetSpectrumDataPoint( float Energy )
    {
        float fIndex = (Energy-Min)/(Delta);
        //return AbsorbCoeffDataPoints[int(fIndex)]; //With a large number of points this is close enough
        
        //Code below linearly interpolates to find the "best" value.
        //Will crash if x = MaxX. Fix either by checking for this or bumping an extra value onto the end
        //of the FormFactorDataPoints array
        
        //bounds checking in case we don't have data for particular energies
        if(int(fIndex) >= (int(SpectrumDataPoints.size()) - 1) )
        {
            return SpectrumDataPoints[SpectrumDataPoints.size() - 1]; //just pull out the last one
        }
        
        if(int(fIndex) < 0 )
        {
            return SpectrumDataPoints[0];
        }
        
        float weight = fIndex - int(fIndex);
        
        
        
        return lerp(SpectrumDataPoints[int(fIndex)], SpectrumDataPoints[int(fIndex) + 1], weight);
        
        
    }
    
    
    float lerp( float start, float end, float weight)
    {
        return start + ((end-start)*weight);
    }
    
};

#endif
