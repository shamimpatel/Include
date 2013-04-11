//
//  DataLoader.h
//  TestCode
//
//  Created by Shamim Patel on 25/03/2013.
//  Copyright (c) 2013 Shamim Patel. All rights reserved.
//

#ifndef _DATALOADER_H_
#define _DATALOADER_H_


#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <algorithm>


using namespace std;

struct DataPoint
{
    double x;
    std::vector<double> ColumnValues;
};

bool DataPointLssThnComp( DataPoint A, DataPoint B)
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

class DataLoader
{
protected:
    double Min, Max, Delta;
    int NumDataColumns, NumDataPoints;
    std::string FileName;
    
    std::vector< DataPoint > DataPoints;
    
    void LoadData( const char* Filename, bool bSkipFirstLine)
    {
		this->FileName = Filename;
        ifstream datafile;
        datafile.open(Filename, ios::in);
        
        if(datafile.is_open() == false)
        {
            cout << "Error: Failed to open data file: " << Filename << endl;
            exit(1);
        }
        
        
         //excel copy/paste gives \r not \n so need to check for this
        char NewlineChar = '\n', ReturnChar = '\r';
        char EndLineCharacter = '\0';
		do
		{
			datafile.get(EndLineCharacter);
			if(EndLineCharacter == '\0')
			{
				cout << "Error: Could not find end of line character in datafile: " << Filename << endl;
				exit(1);
			}
			/*if( EndLineCharacter == NewlineChar || EndLineCharacter == ReturnChar )
			{
				break;
			}*/
		} while ( EndLineCharacter != NewlineChar && EndLineCharacter != ReturnChar );
		
		/*
		if( EndLineCharacter != NewlineChar && EndLineCharacter != ReturnChar ) //can hit this if file is not properly formatted
		{
			cout << "Error: Could not find end of line character in datafile: " << Filename << endl;
			exit(1);
		}*/
		
        datafile.seekg(0); //seek back to beginning as we have just read a line of the file.

        
        std::string dataline;
        
        std::vector<DataPoint> DataFile;
        
        DataPoint Data;
        
        double FileMin, FileMax;
        bool bFirstIteration = true;
        
		if( bSkipFirstLine )
		{
			getline(datafile, dataline, EndLineCharacter); //skip ahead one line
		}
		
        while(getline(datafile, dataline, EndLineCharacter))
        {
            stringstream linestream(dataline);
            
            linestream >> Data.x;
            
            if(bFirstIteration)
            {
                FileMin = FileMax = Data.x;
                bFirstIteration = false;
            }
			else
			{
				if( Data.x < DataFile.back().x) //file should be in ascending order
				{
					cout << "Error: Datafile not sorted properly (" << Filename << ")" << endl;
					exit(1);
				}
			}
            
            
            if( Data.x < FileMin)
            {
                FileMin = Data.x;
            }
            
            if( Data.x > FileMax)
            {
                FileMax = Data.x;
            }
            
            
            for(int i = 0; i < NumDataColumns; i++)
            {
                double ColumnValue;
                linestream >> ColumnValue;
				
				//Data.ColumnValues.push_back( ColumnValue );
				
				//doing it this way is important. Pushing the double by itself (as above) does not work! (not sure why, something about copy constructors probably)
                Data.ColumnValues.push_back( 0.0 );
                Data.ColumnValues[i] = ColumnValue;
            }
            
            DataFile.push_back( Data );
        }
        
        
        if (FileMin > Min || FileMax < Max )
        {
            cout << "Error: Datafile (" << Filename << ") does not cover requested bounds: " << Min << " --> " << Max << endl;
			cout << "Detected file bounds: " << FileMin << " --> " << FileMax << endl;
            exit(1);
        }
        
        DataPoints.clear();
        
        Delta = (Max-Min)/double(NumDataPoints);
        
        for(float i = 0; i <= NumDataPoints; i++)
        {
            float x = Min + i*Delta;
            
            DataPoint C = FindDataPointFromFile( x, &DataFile);
            
            DataPoints.push_back( C );
        }
        
        datafile.close();
        /*
        for(int i = 0; i < int(DataFile.size()); i++)
        {
            cout << DataFile[i].x << "\t" << DataFile[i].ColumnValues[0] << endl;
        }
        */
    }
    
    DataPoint FindDataPointFromFile( float x, std::vector<DataPoint> *DataFile)
    {
        DataPoint SearchPoint;
        SearchPoint.x = x;
        
        vector<DataPoint>::iterator C = std::lower_bound( DataFile->begin(),
                                                    DataFile->end(),
                                                    SearchPoint,
                                                    DataPointLssThnComp );
        if( C == DataFile->begin() )
        {
            //then we're done
            SearchPoint = *C;
        }
        else
        {
            float weight = (SearchPoint.x - (*(C-1)).x)/((*C).x - (*(C-1)).x); //requested energy is always >= C-1.x
            
            //interpolate between values
            for(int i = 0; i < NumDataColumns; i++)
            {
                double ColumnValue = lerp( (*(C-1)).ColumnValues[i], (*C).ColumnValues[i], weight);
                SearchPoint.ColumnValues.push_back(ColumnValue);
            }
            
            
            //could also interpolate the last one too
        }
        
        return SearchPoint;
    }
    
    float lerp( float start, float end, float weight)
    {
        return start + ((end-start)*weight);
    }
    
public:
    //Min/Max for checking data
    DataLoader( float Min, float Max, int NumDataColumns, int NumDataPoints, const char* Filename, bool bSkipFirstLine = false)
    {
        if( Min > Max)
        {
            cout << "Error: Minimum > Maximum requested when loading " << Filename << endl;
            exit(1);
        }
        
        this->Min = Min;
        this->Max = Max;
        this->NumDataColumns = NumDataColumns;
        this->NumDataPoints  = NumDataPoints;
        this->FileName = FileName;
        
        LoadData(Filename, bSkipFirstLine);
        
    }
    
    DataPoint GetDataPoint( double x )
    {
        double fIndex = (x-Min)/(Delta);
        if(int(fIndex) >= (int(DataPoints.size()) - 1) )
        {
            if( fIndex > (int(DataPoints.size()) - 1) )
            {
                cout << "Error: Requested datapoint (x = " << x << ") out of range from " << FileName << endl;
                exit(1);
            }           
            
            return DataPoints[DataPoints.size() - 1]; //just pull out the last one
        }
        
        if(int(fIndex) < 0 )
        {
            cout << "Error: Requested datapoint (x = " << x << ") out of range from " << FileName << endl;
            exit(1);
            
            //return DataPoints[0];
        }
        
        float weight = fIndex - int(fIndex);
        
        DataPoint Data;
        Data.x = x;
        
        
        for(int i = 0; i < NumDataColumns; i++)
        {
            double ColumnValue = lerp(DataPoints[int(fIndex)].ColumnValues[i], DataPoints[int(fIndex) + 1].ColumnValues[i], weight);
            Data.ColumnValues.push_back(ColumnValue);
        }
        
        return Data;
    }
    
};

#endif
