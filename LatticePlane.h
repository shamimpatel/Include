#ifndef _LATTICEPLANE_H_
#define _LATTICEPLANE_H_

#include "Vector.h"
#include "XRay.h"
#include <cmath> //included by vector.h anyway?
#include "FormFactorData.h"
#include "AbsorbCoeffData.h"
#include "Integration.h"

const float ElectronRadius_PI = 8.9697826E-6; //in Angstroms


class LatticePlane
{
public:
    //reciprocal lattice vectors.
    Vector b1,b2,b3;
    Vector H;
    Vector UnitH; //unit vector normal to lattice plane.
    
        
    int h,k,l;
    
    int Multiplicity;
    
    FormFactorData* FormFacData;
    AbsorbCoeffData* AbsorbCoData;
            
    float UnitCellVol;
    
    double Temperature, DebyeTemperature;
    double mass_amu;
    
    double DebyeWallerPreFactor;
    
public:
    
    float HalfMagH; //magnitude of H/2
        
    
    LatticePlane( Vector b1, Vector b2, Vector b3, int h, int k, int l,
                 FormFactorData* FormFacData, float UnitCellVol, AbsorbCoeffData* AbsorbCoData,
                 int Multiplicity, double DebyeWallerPreFactor)
    {
        this->b1 = b1;
        this->b2 = b2;
        this->b3 = b3;
        H = h*b1 + k*b2 + l*b3;
        UnitH = H.Normalized();
        HalfMagH = 0.5f * H.Magnitude();
        this->FormFacData = FormFacData;
        this->h = h;
        this->k = k;
        this->l = l;  
        this->UnitCellVol = UnitCellVol;
        this->AbsorbCoData = AbsorbCoData;
        this->Multiplicity = Multiplicity;
        
        this->DebyeWallerPreFactor = DebyeWallerPreFactor;
    }
    
    Vector BraggReflectXRay( XRay Ray )
    {
        return BraggReflectXRay( Ray.Direction, Ray.Wavelength);
    }
    
    //if the direction (and therefore angle) don't match up with the plane and wavelength then this is meaningless
    Vector BraggReflectXRay( Vector Direction, float Wavelength)
    {
        return Wavelength*H + Direction;
    }
    
    float FindBraggReflectionAngle( XRay Ray )
    {
        return FindBraggReflectionAngle( Ray.Wavelength );
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
    
    float FindGeometricalReflectionAngle( XRay Ray )
    {
        return FindGeometricalReflectionAngle( Ray.Direction );
    }
    
    float FindGeometricalReflectionAngle( Vector Direction )
    {
        
        //use sin as we want angle between ray and plane, not between ray and normal.
        float angle = asin( Direction.Dot(UnitH ) );
        
        //fabs to push values <0 to >0
        return fabs(angle);
    }
    
    float CalculatePowderScatter( float Wavelength) //mu in A^-1
    {
        float sinTheta = Wavelength * HalfMagH;
        
        if( sinTheta > 1.0f || sinTheta <= 0.0f ) //always positive so no need to check for < -1.0f
        {
            //if the Bragg reflection isn't allowed just return 0;
            return 0.0f;
        }
        
        float BraggTheta = asin( sinTheta );
        
        float f1 = FormFacData->GetFormFactorDataPoint( HalfMagH );
        
        float r_e = 2.818E-5; //classical electron radius in A
                
        //mu is in A^-1
        float mu = AbsorbCoData->GetAbsorbCoeffDataPoint( Wavelength );
        
        float density = 16.69f; // in g/cm^3
        
                     
        //1E8 as mu is in A^-1. Want cm^-1 which then gets converted to cm^2/g as it is divided by the density
        //the 232.54 comes from 4pi*Avagadro*electronradius*hbar*c / 1keV * Atomic mass (in amu)
        float f2 = ( float(1E8) * (mu/density) * (12.39842f/Wavelength) )/(232.54f);
        
        
        float ModFSquare = 4*(f1*f1 + f2*f2); // always +ve so no need for fabs
        
        float Cos2Theta = cos(2*BraggTheta);
        
        return (Multiplicity*r_e*r_e*Wavelength*Wavelength*Wavelength*ModFSquare*(1+(Cos2Theta*Cos2Theta)))/(16*UnitCellVol*UnitCellVol*mu*sinTheta);
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

};

double DebyeFunc( double z )
{
    if(z == 0.0)
    {
        return 0.0;
    }
    
    return z/(exp(z)-1);
}

//See Eq 11.77, page 190 of Warren for where this comes from.
double CalculateDebyeWallerPreFactor(double Temp, double DebyeTemp, double mass_amu)
{
    double x = DebyeTemp/Temp;
    double IntegralPart = Compute1DIntegral(DebyeFunc, 0, x, 0.00001)*(1/x) + 0.25*x;
    double PreFactor = (IntegralPart*22998.4*Temp)/(mass_amu*DebyeTemp*DebyeTemp);
    return PreFactor;
}



#endif
