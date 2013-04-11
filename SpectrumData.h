//
//  SpectrumData.h
//  TestCode
//
//  Created by Shamim Patel on 06/03/2013.
//  Copyright (c) 2013 Shamim Patel. All rights reserved.
//

#ifndef _SPECTRUMDATA_H_
#define _SPECTRUMDATA_H_
/*
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>*/

#include "Integration.h"
#include "DataLoader.h"


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

class SpectrumData : protected DataLoader
{

    
public:
    
    
    SpectrumData(float MinE, float MaxE, int NumDataColumns, int NumDataPoints, const char* Filename) : DataLoader( MinE, MaxE, NumDataColumns, NumDataPoints, Filename)
    {
    }
    
    
public:
      
    float GetSpectrumDataPoint( float Energy )
    {
        return this->GetDataPoint(Energy).ColumnValues[0];
    }
    
    double IntegrateSpectrum( float MinE, float MaxE, float dE)
    {
        double Sum = 0.0;
        int nSteps = (MaxE - MinE)/dE;
        
        double x0 = GetSpectrumDataPoint(MinE);
        double xN = GetSpectrumDataPoint(MaxE);
        
        double ExtraTerms = (x0 + xN)*dE*0.5;
        
        for( int Step = 1; Step < nSteps; Step++)
        {
            double x = MinE + dE*Step;
            Sum += GetSpectrumDataPoint(x)*dE;
        }
        return Sum + ExtraTerms;        
    }
    
    
    void GetEnergyIntervals( int numIntervals, std::vector<std::pair<double, double> > *Intervals)
    {
        Intervals->clear();
        double dE = 0.0001;
        double TotalIntegratedSpectrum = IntegrateSpectrum(Min,Max,0.0001);
        
        cout << "Total Integrated Spectrum:\t" << TotalIntegratedSpectrum << endl;
        
        double IntervalMinimum = Min;
        double E = Min;
        double Sum = 0.0;
        
        while( int(Intervals->size()) < numIntervals )
        {
            double RequiredIntegral = TotalIntegratedSpectrum*double(Intervals->size() + 1)/double(numIntervals);
            
            std::pair<double,double> Interval(IntervalMinimum ,0);
            
            while( (Sum) < RequiredIntegral)
            {
                Sum += GetSpectrumDataPoint(E) * dE;
                E += dE;
                
                if(E > Max)
                {
                    break;
                }                
            }            
            Interval.second = E;
            Intervals->push_back(Interval);
            IntervalMinimum = E;
			cout << "Created " << Intervals->size() << " intervals" << endl;
        }
        
        if( (*Intervals)[Intervals->size() - 1].second > Max)
        {
            (*Intervals)[Intervals->size() - 1].second = Max;
        }
        
    }
    

    
};

#endif
