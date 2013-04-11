#ifndef _POWDERDIFFRACTION_H_
#define _POWDERDIFFRACTION_H_

#include <iostream>
#include "Vector.h"
#include "math.h"
#include "FormFactorData.h"
//#include "LatticePlane.h"
#include "AbsorbCoeffData.h"
//#include <float.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

using namespace std;

#ifndef base_generator_type
  typedef boost::mt19937 base_generator_type;
#endif


int CharToInt( char c )
{
    c = c - '0'; //this works as on an ascii table the integers are one after the other. Subtracting '0' gives the right result. Should be platform independant.
    int i = c;
    return i;
}

struct DiffractionDataLine
{
    float Energy;
    std::vector< float > BraggAngles; //radians
    std::vector< float > ScatterProbs;
    std::vector< float > RockingCurves;
    float ModifiedScatterProb;
};


class PowderDiffraction
{

    boost::variate_generator<base_generator_type&, boost::uniform_real<> > *uni;

 public:
    PowderDiffraction( float LatticeConst, boost::variate_generator<base_generator_type&, boost::uniform_real<> > *Generator )
    {
        this->a0 = LatticeConst;
        this->uni = Generator;
    }


    float a0;
    float MinE;
    float MaxE;
    float DeltaE;
    int   NumLatticePlanes;

    std::vector< DiffractionDataLine > DiffDataLines;

    void LoadData( const char* Filename )
    {
        ifstream Data( Filename );

        std::string S;
        getline( Data, S, '\n');

        stringstream stream(S);


        std::vector<std::string> PlaneNames;

        string s;
        stream >> s; //ditch the first one (just reads "Energy")

        while( stream >> s )
        {
            PlaneNames.push_back( s );
        }

        PlaneNames.erase( --PlaneNames.end() ); //last column is "Sum"

        std::vector< float > LatticeSpacings;

        MPI_cout << "Found the following Lattice Planes:\t";

        for( unsigned int i = 0; i < PlaneNames.size(); i++) //loop through the planes and pull out miller indices
        {
            int h = CharToInt(PlaneNames[i][0]);
            int k = CharToInt(PlaneNames[i][1]);
            int l = CharToInt(PlaneNames[i][2]);

            MPI_cout << h << k << l << "\t";

            float LatticeSpacing = a0 / sqrt( float(h)*float(h) + float(l)*float(l) + float(k)*float(k) );
            LatticeSpacings.push_back( LatticeSpacing );
        }

        MPI_cout << endl;

        while(getline(Data, s, '\n'))
        {
            stringstream linestream(s);

            DiffractionDataLine L;

            linestream >> L.Energy;

            for( int i = 0; i < int(LatticeSpacings.size()); i++)
            {
                float BraggAngle, ScatterProb, RockingCurve;
                linestream >> BraggAngle >> ScatterProb >> RockingCurve;
                L.BraggAngles.push_back( BraggAngle );
                L.ScatterProbs.push_back( ScatterProb );
                L.RockingCurves.push_back( RockingCurve );
            }

            linestream >> L.ModifiedScatterProb; //total scattering probability for a lattice plane (technicaly set of lattice planes since multiplicity is included)

            DiffDataLines.push_back( L );
        }

        Data.close();

        MinE = MaxE = DiffDataLines[0].Energy;

        for(int i = 0; i < int(DiffDataLines.size()); i++ )
        {
            if(DiffDataLines[i].Energy < MinE)
            {
                MinE = DiffDataLines[i].Energy;
            }

            if(DiffDataLines[i].Energy > MaxE)
            {
                MaxE = DiffDataLines[i].Energy;
            }
        }

        DeltaE = (MaxE-MinE)/float(DiffDataLines.size()-1);

        NumLatticePlanes = int(DiffDataLines[0].BraggAngles.size());

        MPI_cout << "DiffractionData: MinE: " << MinE << "\tMaxE: " << MaxE << "\tDeltaE: " << DeltaE << endl;
    }

    float PickBraggScatteringAngle( float Energy, float &RockingCurve )
    {
        float fIndex = (Energy-MinE)/DeltaE;
        int Index = int(fIndex); //Energy Index

        //float R = UniformRand();
        float R = (*uni)();
        float ScatterProbSum = 0.0;
        for( int i = 0; i < NumLatticePlanes; i++) //loop through each lattice plane
        {
            float ScatterProb;
            
            //Get the scattering probability of this Lattice plane
            float weight = fIndex - float(Index);
            if(Index >= (int(DiffDataLines.size()) - 1) )
            {
                ScatterProb = DiffDataLines[DiffDataLines.size() - 1].ScatterProbs[i]; //just pull out the last one
            }
            else if(Index < 0 )
            {
                ScatterProb = DiffDataLines[0].ScatterProbs[i];
            }
            else
            {                
                ScatterProb = lerp(DiffDataLines[Index].ScatterProbs[i], DiffDataLines[Index + 1].ScatterProbs[i], weight);
            }

            ScatterProbSum += ScatterProb; //add it to the total
            
            if( R < ScatterProbSum)
            {
                if(Index >= (int(DiffDataLines.size()) - 1) )
                {
                    RockingCurve = DiffDataLines[DiffDataLines.size() - 1].RockingCurves[i];
                    return DiffDataLines[DiffDataLines.size() - 1].BraggAngles[i]; //just pull out the last one
                }

                if(Index < 0 )
                {
                    RockingCurve = DiffDataLines[0].RockingCurves[i];
                    return DiffDataLines[0].BraggAngles[i];
                }

                RockingCurve = lerp(DiffDataLines[Index].RockingCurves[i], DiffDataLines[Index + 1].RockingCurves[i], weight);
                return lerp(DiffDataLines[Index].BraggAngles[i], DiffDataLines[Index + 1].BraggAngles[i], weight);
            }
        }
        
        //we might get here as after processing the probability of scattering into a lattice plane that we care about (for this particular energy) may not add up to 1;        
        return -1; //Scatter into angle that we don't care about
    }

    float GetModifiedScatterProb( float Energy )
    {
        float fIndex = (Energy-MinE)/DeltaE;

        if(int(fIndex) >= (int(DiffDataLines.size()) - 1) )
        {
            if( fIndex > (int(DiffDataLines.size()) - 1) )
            {
                cout << "Error: Requested datapoint out of range for diffraction probability data." << endl;
                exit(1);
            }            
            
            return DiffDataLines[DiffDataLines.size() - 1].ModifiedScatterProb; //just pull out the last one
        }

        if(int(fIndex) < 0 )
        {
            cout << "Error: Requested datapoint out of range for diffraction probability data." << endl;
            exit(1);
            
            //return DiffDataLines[0].ModifiedScatterProb;
        }

        float weight = fIndex - int(fIndex);

        return lerp(DiffDataLines[int(fIndex)].ModifiedScatterProb, DiffDataLines[int(fIndex) + 1].ModifiedScatterProb, weight);
    }

    float lerp( float start, float end, float weight)
    {
        return start + ((end-start)*weight);
    }

};

#endif
