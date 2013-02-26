#ifndef _XRAY_H_
#define _XRAY_H_


#include "Vector.h"


float EnergyToWavelength( float Energy );
float WavelengthToEnergy( float Wavelength );



class XRay
{
 public:
  Vector Source, Direction;
  float Wavelength;
  XRay( Vector Source, Vector Direction, float Wavelength)
  {
    this->Source = Source;
    this->Direction = Direction;
    this->Wavelength = Wavelength;
  }
};



float EnergyToWavelength( float Energy )
{
  //Energy in keV. Wavelength in Angstroms
  return 12.39842/Energy; 
}


float WavelengthToEnergy( float Wavelength )
{
  return 12.39842/Wavelength;  
}



#endif
