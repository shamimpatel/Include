//
//  FluorescenceData.h
//  TestCode
//
//  Created by Shamim Patel on 05/04/2013.
//  Copyright (c) 2013 Shamim Patel. All rights reserved.
//

#ifndef _FLUORESCENCEDATA_H
#define _FLUORESCENCEDATA_H

#include "FileReading.h"
#include "DataLoader.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>


#ifndef __generator_type
	typedef boost::mt19937 base_generator_type;
	#define __generator_type
#endif


class FluorescenceData
{
private:
	struct ShellEmissionData
	{
		int ShellId;
		double Yield;
		int NumLines;
		struct EmissionLine
		{
			double Energy;
			double RelativeProbability;
		};
		std::vector< EmissionLine > EmissionLines;
	};
	
	
	DataLoader* ShellProbabilities;
	boost::variate_generator<base_generator_type&, boost::uniform_real<> > *uni;
	std::vector< ShellEmissionData > ShellDataArray;
	
	
	DataPoint PrePickedPoint;
	
public:
	FluorescenceData( double MinE, double MaxE, int NumDataPoints, const char* LineDataFileName, const char* ShellProbDataFileName,
						boost::variate_generator<base_generator_type&, boost::uniform_real<> > *Generator )
	{
		this->uni = Generator;
		
		ifstream EmissionLineData( LineDataFileName );
		
		
		if(EmissionLineData.is_open() == false)
        {
            cout << "Error: Failed to open data file: " << LineDataFileName << endl;
            exit(1);
        }	
		
		
		char NewlineChar = '\n', ReturnChar = '\r';
		char EndLineCharacter = '\0';
		//while( EmissionLineData >> EndLineCharacter)
		do
		{
			EmissionLineData.get(EndLineCharacter);
			if(EndLineCharacter == '\0')
			{
				cout << "Error: Could not find end of line character in datafile: " << LineDataFileName << endl;
				exit(1);
			}
			/*
			if( EndLineCharacter == NewlineChar || EndLineCharacter == ReturnChar )
			{
				break;
			}
			 */
		} while ( EndLineCharacter != NewlineChar && EndLineCharacter != ReturnChar );
		
	
		EmissionLineData.seekg(0); //seek back to beginning as we have just read a line of the file.
				
		
		
		std::string dataline;
		while(getline(EmissionLineData, dataline, EndLineCharacter))
		{
			stringstream linestream(dataline);
			std::string Val;
			
			ShellEmissionData ShellData;
			
			linestream >> Val;
			IntFromString(Val, ShellData.ShellId);
			
			linestream >> Val;
			DoubleFromString(Val, ShellData.Yield);
			
			linestream >> Val;
			IntFromString(Val, ShellData.NumLines);
			
			if (ShellData.NumLines <= 0)
			{
				cout << "Error: Number of emission lines <= 0" << endl;
				exit(1);
			}
			
			//cout << ShellData.NumLines << endl;
			
			for( int i = 0; i < ShellData.NumLines; i++)
			{
				ShellEmissionData::EmissionLine Line;
				
				linestream >> Val;
				DoubleFromString(Val, Line.Energy);
				
				linestream >> Val;
				DoubleFromString(Val, Line.RelativeProbability);
				
				
				ShellData.EmissionLines.push_back(Line);
			}
			
			ShellDataArray.push_back(ShellData);
		}
		
		EmissionLineData.close();
		
		ShellProbabilities = new DataLoader( MinE, MaxE, int(ShellDataArray.size()), NumDataPoints, ShellProbDataFileName, true );
		
		/*for (auto it = ShellDataArray.begin(); it != ShellDataArray.end(); ++it)
		 {
		 cout << it->ShellId << "\t" << it->Yield << "\t" << it->NumLines << "\t";
		 
		 for( auto it2 = it->EmissionLines.begin(); it2 != it->EmissionLines.end(); ++it2)
		 {
		 cout << it2->Energy << "\t" << it2->RelativeProbability << "\t";
		 }
		 
		 cout << endl;		
		 }*/
	}
	
	void PreSelectEnergy( double Energy )
	{
		PrePickedPoint = ShellProbabilities->GetDataPoint(Energy);
	}
	
	double PickFluorescenceEnergy( double Energy )
	{
		//DataPoint Point = ShellProbabilities->GetDataPoint(Energy);
				
		PrePickedPoint = ShellProbabilities->GetDataPoint(Energy);
		
		return PickFluorescenceEnergy();
	}
	
	double PickFluorescenceEnergy()
	{
		double R = (*uni)();
		double ProbSum = 0.0;
		for( int ShellNum=0; ShellNum < int(ShellDataArray.size()); ShellNum++)
		{
			ProbSum += PrePickedPoint.ColumnValues[ShellNum];
			if(R < ProbSum)
			{
				if((*uni)() < ShellDataArray[ShellNum].Yield)
				{
					
					double R2 = (*uni)();
					double ProbSum2 = 0.0;
					
					for (int EmissionLineNum = 0; EmissionLineNum < ShellDataArray[ShellNum].NumLines; EmissionLineNum++)
					{
						ProbSum2 += ShellDataArray[ShellNum].EmissionLines[EmissionLineNum].RelativeProbability;
						if( R2 <= ProbSum2 )
						{
							return ShellDataArray[ShellNum].EmissionLines[EmissionLineNum].Energy;
						}
					}
				}
				else
				{
					return -1.0;
				}
			}
		}
		return -1.0;
	}
	
	
	
};


#endif
