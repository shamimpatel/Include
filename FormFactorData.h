#ifndef _FORMFACTORDATA_H_
#define _FORMFACTORDATA_H_




#include "DataLoader.h"

class FormFactorData : public DataLoader
{
public:
    
    FormFactorData( float Min, float Max, int NumDataColumns, int NumDataPoints, const char* Filename): DataLoader( Min, Max, NumDataColumns, NumDataPoints, Filename)
	{
		
	}
    ~FormFactorData()
	{
		
	}
	
    double GetFormFactorDataPoint( double x)
    {
        return GetDataPoint(x).ColumnValues[0];
    }
    
};


#endif


