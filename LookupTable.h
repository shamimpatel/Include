//
//  LookupTable.h
//  TestCode
//
//  Created by Shamim Patel on 19/06/2013.
//  Copyright (c) 2013 Shamim Patel. All rights reserved.
//

#ifndef _LOOKUPTABLE_H_
#define _LOOKUPTABLE_H_

#include <iostream>
#include <cmath>
#include <vector>

#ifdef USING_MPI
	#include "MPIUtils.h"
#else
	#define MPI_cout cout
#endif


using namespace std;

#define UNSAFE_LOOKUPS
//#define NO_INTERP
//#define USE_ARRAY //In a release mode compile there doesn't seem to be any benefit to using an array instead of a std::vector (they probably get optimised into almost the same code)


#ifndef PI
#define PI (3.14159265)
#endif


#ifdef LOOKUPTABLE_USE_SINGLE_PRECISION
typedef float _TablePrecision;
#else
typedef double _TablePrecision;
#endif

class LookupTable
{
	_TablePrecision Min, Max;
	int NumDataPoints;
	_TablePrecision Delta;
	
#ifndef USE_ARRAY
	std::vector<_TablePrecision> Data;
#else
	_TablePrecision* Data;
#endif
		
	_TablePrecision lerp( _TablePrecision start, _TablePrecision end, _TablePrecision weight)
    {
        return start + ((end-start)*weight);
    }
	
public:
	LookupTable( _TablePrecision Min, _TablePrecision Max, int NumDataPoints, _TablePrecision (*f)(_TablePrecision))
	{
		this->Min = Min;
		this->Max = Max;
		this->NumDataPoints = NumDataPoints;
		Delta = (Max-Min)/double(NumDataPoints-1);
		
#ifdef USE_ARRAY
		Data = (_TablePrecision*)malloc(sizeof(_TablePrecision)*(NumDataPoints+1));
#endif
		
		for( int i = 0; i < NumDataPoints; i++)
		{
			_TablePrecision x = Min + Delta*i;
			
#ifndef USE_ARRAY
			//Sometimes the compiler (GCC appears to do it more often) seems to optimize Data.push_back(x) to just taking a reference to x so all data members are the same.
			//I can't reliably reproduce this but doing the following does seem to clear it up.
			Data.push_back(0.0);
			Data.back() = f(x);
#else
			Data[i] = f(x);
#endif
		}
		
#ifdef UNSAFE_LOOKUPS
		MPI_cout << "Warning: Lookup tables will not check for accesses out of initialization bounds." << endl;
#endif
#ifdef NO_INTERP
		MPI_cout << "Warning: Lookup results are not interpolated!" << endl;
#endif
		
		//Add an extra element to the array so that we avoid any chance of looking up a value outside of the array.
		//The bounds checking should avoid this but I'm not certain it will work everywhere due to precision issues.
		//Without bounds checking is is guaranteed to happen for looking up the maximum value.
		
#ifndef USE_ARRAY
		//Bounds checking will handle this however without them looking up the Maximum value it will attempt to lerp between the last element and whatever is next in memory so adding the same thing again to the array will prevent this.
		_TablePrecision FinalVal = Data.back();
		Data.push_back(0.0);
		Data.back() = FinalVal;
		Data.shrink_to_fit();
#else
		_TablePrecision FinalVal = Data[NumDataPoints-1];
		Data[NumDataPoints] = FinalVal;
#endif
		
	}
	
	~LookupTable()
	{
#ifdef USE_ARRAY
		free(Data);
#endif
	}
	
	_TablePrecision Lookup( _TablePrecision x )
	{
		_TablePrecision fIndex = (x-Min)/Delta;
#ifndef UNSAFE_LOOKUPS
		if( fIndex < 0.0 || fIndex > _TablePrecision(NumDataPoints-1))
		{
			MPI_cout << "Error: Requested point (" << x << ") is out of lookup table range ("<< Min << " --> " << Max << ")" << endl;
			exit(1);
		}
		else if( fIndex == (NumDataPoints-1) ) //not sure about double comparison to int. Does work for small numbers but probably won't do well for bigger ones.
		{
			return Data[NumDataPoints-1];
		}
#endif
		
#ifndef NO_INTERP
		_TablePrecision weight = fIndex - int(fIndex);
		return lerp( Data[int(fIndex)], Data[int(fIndex) +1], weight );
#else
		return Data[ int(fIndex) ];
#endif
	}
};

class CosLookupTable: public LookupTable
{
public:
	CosLookupTable( int NumDataPoints ):LookupTable(0.0, 2.0*PI, NumDataPoints, cos){}
	
	_TablePrecision Lookup( _TablePrecision Theta )
	{
		Theta = fabs(Theta);
		int NumPeriods = int( Theta/(2.0*PI) );
		Theta -= (2.0*PI)*NumPeriods;

		return LookupTable::Lookup(Theta);
	}
};

class SinLookupTable: public LookupTable
{
public:
	SinLookupTable( int NumDataPoints ):LookupTable(0.0, 2.0*PI, NumDataPoints, sin){}
	
	_TablePrecision Lookup( _TablePrecision Theta )
	{
		int NumPeriods = abs(int(Theta/(2.0*PI)));
		if( Theta < 0)
		{
			Theta += (2.0*PI)*(NumPeriods+1);			
		}
		else
		{
			Theta -= (2.0*PI)*NumPeriods;
		}
		return LookupTable::Lookup(Theta);
	}
};


#endif
