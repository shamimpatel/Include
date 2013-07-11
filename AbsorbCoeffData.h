#ifndef _ABSORBCOEFFDATA_H_
#define _ABSORBCOEFFDATA_H_


#include "DataLoader.h"
/*
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
*/
class AbsorbCoeffDataEnergy : public DataLoader
{
public:
    
    AbsorbCoeffDataEnergy( float MinEnergy, float MaxEnergy, int NumDataPoints, const char* Filename):
    DataLoader( MinEnergy, MaxEnergy, 1, NumDataPoints, Filename){}
    
    double GetAbsorbCoeffDataPointEnergy( double Energy )
    {
        return GetDataPoint(Energy).ColumnValues[0];
    }
    
};

#endif

