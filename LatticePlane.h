#ifndef _LATTICEPLANE_H_
#define _LATTICEPLANE_H_

#include "Vector.h"
#include "XRay.h"
#include <cmath> //included by vector.h anyway?
#include "FormFactorData.h"
#include "AbsorbCoeffData.h"

const float ElectronRadius_PI = 8.9697826E-6; //in Angstroms


class LatticePlane
{
 public:
  //reciprocal lattice vectors.
  Vector b1,b2,b3;
  Vector H;
  Vector UnitH; //unit vector normal to lattice plane.

  float NegElectronRadius_PI_V; // -1* electron radius/(pi * unit cell vol)

  int h,k,l;

  int Multiplicity;

  FormFactorData* FormFacData;
  AbsorbCoeffData* AbsorbCoData;

  float F0Real;
  //float F0Im;

  float UnitCellVol;



 public:

  float HalfMagH; //magnitude of H/2


  LatticePlane( Vector b1, Vector b2, Vector b3, int h, int k, int l,
		FormFactorData* FormFacData, float UnitCellVol, AbsorbCoeffData* AbsorbCoData,
		int Multiplicity)
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
    this->NegElectronRadius_PI_V = (-1.0f*ElectronRadius_PI)/UnitCellVol;
    this->F0Real = FormFacData->GetFormFactorDataPoint( 0.0f );
    this->UnitCellVol = UnitCellVol;
    this->AbsorbCoData = AbsorbCoData;
    this->Multiplicity = Multiplicity;
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

    float F1 = FormFacData->GetFormFactorDataPoint( HalfMagH );

    float r_e = 2.818E-5; //classical electron radius in A

    //float M = 6; //multiplicity. TODO: actually calculate this

    //mu is in A^-1
    float mu = AbsorbCoData->GetAbsorbCoeffDataPoint( Wavelength );

    float density = 16.69f; // in g/cm^3

    //1E8 as mu is in A^-1. Want cm^-1 which then gets converted to cm^2/g
    float F2 = ( 1E8 * (mu/density) * (12.39842f/Wavelength) )/(232.54f);

    float ModFSquare = F1*F1 + F2*F2; // always +ve so no need for fabs

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


  Vector GeometricalReflectXRay( Vector Direction )
  {
    return -2*( Direction.Dot( UnitH ) * UnitH ) + Direction;
  }

};


#endif
