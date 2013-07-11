#ifndef _LATTICEPLANE_H_
#define _LATTICEPLANE_H_

#include "Vector.h"
#include "XRay.h"
#include <cmath> //included by Vector.h anyway?
#include "FormFactorData.h"
#include "AbsorbCoeffData.h"
#include "Integration.h"
#include <map>
#include "FileReading.h"

//double DebyeFunc( double z );

class LatticePlane
{
private:
	struct _Atom
	{
		int Type;
		Vector Position;
	};
	
	
	struct _AtomData
	{
		//int AtomType;
		//double RealPart;
		//double ImaginaryPart;
		std::string FormFacDataFilename;
		FormFactorData* FormFacData;
		std::string AbsorbCoeffDataFilename;
		AbsorbCoeffDataEnergy* AbsorbCoData;
		double mass_amu;
		double density;
		
		
		_AtomData()
		{
			//AtomType = -1;
			//RealPart = 0.0;
			//ImaginaryPart = 0.0;
		}
	};
	
	static std::vector< _Atom > UnitCell;
	static std::map< int, _AtomData > AtomData;
	
	
public:
	static Vector a1,a2,a3; //Real Lattice Vectors
    static Vector b1,b2,b3; //reciprocal lattice vectors.
    Vector H;
    Vector UnitH; //unit vector normal to lattice plane.
	
    int h,k,l;
    
    int Multiplicity;
    
    //FormFactorData* FormFacData;
    static AbsorbCoeffDataEnergy* AbsorbCoData;
	
    static double UnitCellVol;
    
    double Temperature, DebyeTemperature;
    double mass_amu;
    
    static double DebyeWallerPreFactor;
	
	
	static double DebyeFunc( double z )
	{
		if(z == 0.0)
		{
			return 0.0;
		}
		
		return z/(exp(z)-1);
	}
    
	static double MinE, MaxE;
	
public:
    
    float HalfMagH; //magnitude of H/2
	
    
	static void SetDebyeWallerPreFactor( double PreFactor )
	{
		DebyeWallerPreFactor = PreFactor;
	}
	
	static void SetupUnitCell( const char* Filename, Vector a1, Vector a2, Vector a3, double LatticeConstant,
							  double MinE, double MaxE, AbsorbCoeffDataEnergy* AbsorbCoData )
	{
		LatticePlane::a1 = a1;
		LatticePlane::a2 = a2;
		LatticePlane::a3 = a3;
		
		LatticePlane::UnitCellVol = a2.Cross(a3).Dot(a1);
		
		b1 = (a2.Cross(a3))/UnitCellVol;
		b2 = (a3.Cross(a1))/UnitCellVol;
		b3 = (a1.Cross(a2))/UnitCellVol;
		
		LatticePlane::MinE = MinE;
		LatticePlane::MaxE = MaxE;
		
		
		ifstream datafile;
		
		datafile.open(Filename, ios::in);
		
		if(datafile.is_open() == false)
		{
			cout << "Error: Failed to open data file: " << Filename << endl;
			exit(1);
		}
		
		char NewlineChar = '\n', ReturnChar = '\r';
		char Character = '\0';
		do
		{
			datafile.get(Character);
			if(Character == '\0')
			{
				cout << "Error: Could not find end of line character in datafile: " << Filename << endl;
				exit(1);
			}
		}
		while ( Character != NewlineChar && Character != ReturnChar );
		
		datafile.seekg(0); //seek back to beginning as we have just read a line of the file.
		
		
		std::string dataline;
		//Pull data from file
		while(getline(datafile, dataline, Character))
		{
			stringstream linestream(dataline);
			std::string CommandWord;
			
			linestream >> CommandWord;
			
			if( CommandWord == std::string("Atom") )
			{
				std::string Val;
				linestream >> Val;
				
				int AtomType;
				IntFromString(Val, AtomType);
				
				linestream >> Val;
				Vector Position;
				if(!VectorFromString(Val, Position))
				{
					cout << "Error: Atom string (\"" << Val << "\")not properly formatted" << endl;
					exit(1);
				}
				_Atom Atom;
				Atom.Type = AtomType;
				Atom.Position = Position;
				UnitCell.push_back( Atom );
			}
			
			
			if( CommandWord == std::string("Type") )
			{
				std::string Val;
				linestream >> Val;
				
				int AtomType;
				IntFromString(Val, AtomType);
				
				if (AtomData.count(AtomType) != 0)
				{
					cout << "Error: Duplicate atom type id found (" << AtomType << ")!" << endl;
					exit(1);
				}
				
				linestream >> Val;
				_AtomData AtomicData;
				AtomicData.FormFacDataFilename = Val;
				
				linestream >> Val;
				AtomicData.AbsorbCoeffDataFilename = Val;
				
				linestream >> Val;
				double mass;
				DoubleFromString(Val, mass);
				AtomicData.mass_amu = mass;
				
				linestream >> Val;
				double density;
				DoubleFromString(Val, density);
				AtomicData.density = density;
				
				double MaxFormFactor = 1.0/EnergyToWavelength(MaxE);
				MaxFormFactor += MaxFormFactor*0.02; // Add a 2% tolerance just to cover numerical errors.
				
				FormFactorData* FormFactor = new FormFactorData( 0.0, MaxFormFactor, 1, 10000, AtomicData.FormFacDataFilename.c_str());
				AtomicData.FormFacData = FormFactor;
				
				//double MinWavelength = EnergyToWavelength(MaxE); MinWavelength -= MinWavelength*0.02;
				//double MaxWavelength = EnergyToWavelength(MinE); MaxWavelength += MaxWavelength*0.02;
				
				AbsorbCoeffDataEnergy* AbsorbData = new AbsorbCoeffDataEnergy( MinE, MaxE, 5000, AtomicData.AbsorbCoeffDataFilename.c_str());
				AtomicData.AbsorbCoData = AbsorbData;
				
				
				AtomData.insert(std::pair<int, _AtomData>(AtomType,AtomicData));
			}
		}
		
		
		//simple check to ensure that atom types are valid
		for( std::vector< _Atom >::iterator it = UnitCell.begin(); it != UnitCell.end(); ++it)
		{
			int AtomType = (*it).Type;
			
			if (AtomData.count(AtomType) != 1)
			{
				cout << "Error: Atom type " << (*it).Type << " not defined" << endl;
				exit(1);
			}
		}
		
		LatticePlane::AbsorbCoData = AbsorbCoData;
		
	}
	
	
    LatticePlane( int h, int k, int l,
                 int Multiplicity)
    {
        H = h*b1 + k*b2 + l*b3;
        UnitH = H.Normalized();
        HalfMagH = 0.5f * H.Magnitude();
        this->h = h;
        this->k = k;
        this->l = l;
        this->Multiplicity = Multiplicity;
    }
    

    
    //if the direction (and therefore angle) don't match up with the plane and wavelength then this is meaningless
    Vector BraggReflectXRay( Vector Direction, float Wavelength)
    {
        return Wavelength*H + Direction;
    }
        
    float FindBraggReflectionAngle( float Wavelength )
    {
        float sinTheta = Wavelength * HalfMagH;
        
        if( sinTheta > 1.0f )
        {
            return -1.0f;
            //return negative for no valid reflection
        }
        else
        {
            return asin( Wavelength * HalfMagH );
        }
        //don't need to conditionally take modulus as both wavelength and |H|/2 are +ve
    }
    
    
    float FindGeometricalReflectionAngle( Vector Direction )
    {
        //use sin as we want angle between ray and plane, not between ray and normal.
        float angle = asin( Direction.Dot(UnitH ) );
        
        //fabs to push values <0 to >0
        return fabs(angle);
    }
    
    float CalculatePowderScatter( float Wavelength )
    {
        float sinTheta = Wavelength * HalfMagH;
        
        if( sinTheta > 1.0f || sinTheta <= 0.0f ) //always positive so no need to check for < -1.0f
        {
            //if the Bragg reflection isn't allowed just return 0;
            return 0.0f;
        }
        
        float BraggTheta = asin( sinTheta );
		
		double TotalRealPart = 0.0;
		double TotalComplexPart = 0.0;
		double Energy = WavelengthToEnergy(Wavelength);
		for( std::vector< _Atom >::iterator it = UnitCell.begin(); it != UnitCell.end(); ++it)
		{
			_AtomData AtomicData = AtomData.find( it->Type )->second;
			double f_real = AtomicData.FormFacData->GetFormFactorDataPoint(HalfMagH);			
			double AbsorbCoeff = AtomicData.AbsorbCoData->GetAbsorbCoeffDataPointEnergy( Energy );
			//4pi*Avagadro number*electronradius*hbar*c / 1keV = 42080.32 cm^2
			double f2ConversionFactor = 42080.32/AtomicData.mass_amu;
			double f_im = ( double(1E8) * (AbsorbCoeff/AtomicData.density) * (12.39842/Wavelength) )/(f2ConversionFactor);
			
			
			double argument = 2.0*PI*(double(h) * it->Position.x  + double(k) * it->Position.y + double(l) * it->Position.z);
			double RealCoeff = cos(argument);
			double ComplexCoeff = sin(argument);
			
			//may be worth re-adding this if f_real or f_im are particularly large? Usually after squaring the structure factor the precision errors that this is meant to correct disappear anyway
			/*if( RealCoeff < 10E-7) //precision errors can produce very small numbers instead of zero
			 {
			 RealCoeff = 0.0;
			 }
			 if( ComplexCoeff < 10E-7)
			 {
			 ComplexCoeff = 0.0;
			 }*/
			
			TotalRealPart += RealCoeff*f_real - ComplexCoeff*f_im;
			TotalComplexPart += ComplexCoeff*f_real + RealCoeff*f_im;
		}
		
        //mu is in A^-1
        double mu = AbsorbCoData->GetAbsorbCoeffDataPointEnergy( Energy );
		
        /*
		 Old way of doing it:
		 float f1 = FormFacData->GetFormFactorDataPoint( HalfMagH );
		 
		 
		 double density = 16.69f; // in g/cm^3
		 //1E8 as mu is in A^-1. Want cm^-1 which then gets converted to cm^2/g as it is divided by the density
		 //the 232.54 comes from 4pi*Avagadro number*electronradius*hbar*c / 1keV * Atomic mass (in amu)  (in square cm [not m])
		 double f2ConversionFactor = 42080.32/180.948;
		 
		 double f2 = ( double(1E8) * (mu/density) * (12.39842/Wavelength) )/(f2ConversionFactor);
		 
		 
		 //double ModFSquare = 4.0*(f1*f1 + f2*f2); // always +ve so no need for fabs
		 */
		
		double ModFSquare = TotalRealPart*TotalRealPart + TotalComplexPart*TotalComplexPart;
		
		//cout << ModFSquare - ModFSquare2 << endl;
        
        float Cos2Theta = cos(2*BraggTheta);
		
        float r_e = 2.818E-5; //classical electron radius in A
		
		double DebyeWallerFactor = CalculateDebyeWallerFactor(Wavelength);
		
        return (DebyeWallerFactor*Multiplicity*r_e*r_e*Wavelength*Wavelength*Wavelength*ModFSquare*(1+(Cos2Theta*Cos2Theta)))/(16*UnitCellVol*UnitCellVol*mu*sinTheta);
    }
    
    
    //TODO: Handle what happens at Theta ~= pi/2 (is that even physical?)
    float GaussianScherrerWidth( float Wavelength, float CrystalLength ) //This gives a FWHM
    {
        float BraggTheta = FindBraggReflectionAngle( Wavelength );
        
        if( BraggTheta < 0.0f )
        {
            return 0.0f; //return a width of zero if no bragg reflection available.
        }
        else
        {
            //NB: braggangle not 2*braggangle AKA scattering angle
            //Taken from warren Page 253.
            //Factor of 2.35 comes from FWHM->gaussian sigma.
            return (0.94f*Wavelength)/(2.35482*CrystalLength*cos(BraggTheta));
        }
    }
    
    float CalculateDebyeWallerFactor( float Wavelength )
    {
        double BraggTheta = FindBraggReflectionAngle( Wavelength );
        
        if( BraggTheta < 0.0f )
        {
            return 0.0f; //return an intensity of zero.
        }
        
        double SinTheta_Lambda = sin(BraggTheta)/Wavelength;
        
		double TwoM = DebyeWallerPreFactor*SinTheta_Lambda*SinTheta_Lambda;
        
        return exp( -1.0*TwoM );
    }
    
    
    Vector GeometricalReflectXRay( Vector Direction )
    {
        return -2.0*( Direction.Dot( UnitH ) * UnitH ) + Direction;
    }
	
	//See Eq 11.77, page 190 of Warren for where this comes from.
	static double CalculateDebyeWallerPreFactor(double Temp, double DebyeTemp, double mass_amu)
	{
		double x = DebyeTemp/Temp;
		double IntegralPart = Compute1DIntegral(DebyeFunc, 0, x, 0.00001)*(1/x) + 0.25*x;
		double PreFactor = (IntegralPart*22998.4*Temp)/(mass_amu*DebyeTemp*DebyeTemp);
		return PreFactor;
	}	
};

double LatticePlane::DebyeWallerPreFactor = 0.0;
double LatticePlane::UnitCellVol = 0.0;
Vector LatticePlane::a1;
Vector LatticePlane::a2;
Vector LatticePlane::a3;
Vector LatticePlane::b1;
Vector LatticePlane::b2;
Vector LatticePlane::b3;
double LatticePlane::MinE = 0.0;
double LatticePlane::MaxE = 0.0;
std::vector< LatticePlane::_Atom > LatticePlane::UnitCell;
std::map< int, LatticePlane::_AtomData > LatticePlane::AtomData;
AbsorbCoeffDataEnergy* LatticePlane::AbsorbCoData = NULL;






#endif
