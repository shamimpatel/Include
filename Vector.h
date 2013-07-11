#ifndef _VECTOR_H_
#define _VECTOR_H_

#include <cmath>
#include <iostream>
#include <stdio.h>
using namespace std;


#ifndef PI
#define PI (3.14159265)
#endif



#ifdef VECTOR_USE_SINGLE_PRECISION
	typedef float _VectorPrecision;
#else
	typedef double _VectorPrecision;
#endif

_VectorPrecision Deg2Rad( _VectorPrecision Degrees )
{
    return Degrees * (PI/180.0);
}

_VectorPrecision Rad2Deg( _VectorPrecision Radians )
{
    return Radians * (180.0/PI);
}


_VectorPrecision Rad2Deg2( _VectorPrecision Radians )
{
    if( Radians <= 0.0f)
    {
        return Radians;
    }
    else
    {
        return Radians * (180.0f/PI);
    }
}

_VectorPrecision Deg2Rad2( _VectorPrecision Degrees )
{
    if( Degrees <= 0.0f)
    {
        return Degrees;
    }
    else
    {
        return Degrees * (PI/180.0f);
    }
}



class Vector
{
public:
    _VectorPrecision x,y,z;

    Vector ();
    Vector ( _VectorPrecision X, _VectorPrecision Y, _VectorPrecision Z);
    Vector ( _VectorPrecision r, _VectorPrecision theta, _VectorPrecision phi, bool bSpherical)
    {
        x = r*sin(theta)*cos(phi);
        y = r*sin(theta)*sin(phi);
        z = r*cos(theta);
    }

    // slightly more efficient way to create a spherical unit vector.
    //Need the extra input so that there's no compiler error.
    //(float,float,bool) is apparently equivalent to (float,float,float)
    //might be a GCC thing?
    Vector ( _VectorPrecision theta, _VectorPrecision phi, bool bSpherical, bool bUnit)
    {
        x = sin(theta)*cos(phi);
        y = sin(theta)*sin(phi);
        z = cos(theta);
    }

    _VectorPrecision GetTheta()
    {
        return acos(z/Magnitude());
    }

    void SetTheta( _VectorPrecision theta )
    {
        _VectorPrecision r = Magnitude();
        _VectorPrecision phi = GetPhi();

        x = r*sin(theta)*cos(phi);
        y = r*sin(theta)*sin(phi);
        z = r*cos(theta);
    }

    _VectorPrecision GetPhi()
    {
        return atan2(y,x);
    }

    void SetPhi( _VectorPrecision phi)
    {
        _VectorPrecision r = Magnitude();
        _VectorPrecision theta = GetTheta();

        x = r*sin(theta)*cos(phi);
        y = r*sin(theta)*sin(phi);
        z = r*cos(theta);
    }

    _VectorPrecision Dot( const Vector V )
    {
        return x*V.x + y*V.y + z*V.z;
    }

    Vector Cross ( const Vector V)
    {
        return Vector( y*V.z - z*V.y , z*V.x - x*V.z , x*V.y - y*V.x);
    }

    _VectorPrecision Magnitude()
    {
        return sqrt(x*x + y*y + z*z);
    }

    Vector Normalized()
    {
        _VectorPrecision Mag = this->Magnitude();
        if(Mag == 0.0)
        {
            return Vector(0,0,0); //avoid the divide by 0 error
        }

        return this->Divide(Mag);
    }
	
	void SelfNormalize()
	{
		_VectorPrecision Mag = this->Magnitude();
        if(Mag == 0.0)
        {
            this->x = 0.0;
			this->y = 0.0;
			this->z = 0.0;
        }
		else
		{
			this->x = this->x / Mag;
			this->y = this->y / Mag;
			this->z = this->z / Mag;
		}
	}

    Vector Scale( const _VectorPrecision Scale )
    {
        Vector Out;
        Out.x = x*Scale;
        Out.y = y*Scale;
        Out.z = z*Scale;
        return Out;
    }

    Vector operator+ (Vector V)
    {
        return this->Add(V);
    }

    Vector operator- (Vector V)
    {
        return this->Subtract(V);
    }

    Vector operator* (_VectorPrecision Scale)
    {
        return this->Scale( Scale );
    }


    Vector Divide( _VectorPrecision Denominator )
    {
        Vector Out;
        Out.x = x/Denominator;
        Out.y = y/Denominator;
        Out.z = z/Denominator;
        return Out;
    }

    Vector operator/ (_VectorPrecision Denominator)
    {
        return this->Divide( Denominator );
    }


    Vector Add( Vector V)
    {
        Vector Out;
        Out.x = x + V.x;
        Out.y = y + V.y;
        Out.z = z + V.z;
        return Out;
    }

    Vector Subtract( Vector V)
    {
        Vector Out;
        Out.x = x - V.x;
        Out.y = y - V.y;
        Out.z = z - V.z;
        return Out;
    }

    void Print(bool newline = true)
    {
        if(newline)
        {
            printf("%f\t%f\t%f\n", x,y,z);
        }
        else
        {
            printf("%f\t%f\t%f", x,y,z);
        }
    }

};

Vector operator*(_VectorPrecision Scale, Vector V)
{
    return V*Scale;
}

Vector::Vector(_VectorPrecision x, _VectorPrecision y, _VectorPrecision z)
{
    this->x = x;
    this->y = y;
    this->z = z;
}

Vector::Vector(void)
{
    x = 0;
    y = 0;
    z = 0;
}

//MUST be a unit vector.
Vector ScatterUnitVector( Vector In, _VectorPrecision Theta, _VectorPrecision Phi)
{
    Vector Out;

    _VectorPrecision CT = cos(Theta);
    _VectorPrecision ST = sin(Theta);
    _VectorPrecision CP = cos(Phi);
    _VectorPrecision SP = sin(Phi);

    //http://en.wikipedia.org/wiki/Monte_Carlo_method_for_photon_transport#Step_3:_Absorption_and_scattering
    //And: NUCLEAR SCIENCE AND ENGINEERING: 131, 132–136 ~1999 Direction Cosines and Polarization Vectors
    //for Monte Carlo Photon Scattering

    _VectorPrecision zz = sqrt(1.0-In.z*In.z);

    if( fabs(zz) <= 0.001 ) //if |z|~1
    {
        Out.x = ST*CP;
        Out.y = ST*SP;
        Out.z = CT*In.z/fabs(In.z);
        return Out;
    }

    Out.x = In.x*CT + (1.0/zz)*ST*(In.x*In.z*CP-In.y*SP);
    Out.y = In.y*CT + (1.0/zz)*ST*(In.y*In.z*CP+In.x*SP);
    Out.z = In.z*CT - CP*ST*zz;

    return Out;
}

Vector RotateX( Vector V, _VectorPrecision Angle)
{
    Vector VOut;

    _VectorPrecision CosTheta = cos(Angle);
    _VectorPrecision SinTheta = sin(Angle);

    VOut.x = V.x;
    VOut.y = CosTheta*V.y - SinTheta*V.z;
    VOut.z = SinTheta*V.y + CosTheta*V.z;

    return VOut;
}

Vector RotateY( Vector V, _VectorPrecision Angle)
{
    Vector VOut;

    _VectorPrecision CosTheta = cos(Angle);
    _VectorPrecision SinTheta = sin(Angle);

    VOut.x = CosTheta*V.x + SinTheta*V.z;
    VOut.y = V.y;
    VOut.z = CosTheta*V.z - SinTheta*V.x;

    return VOut;
}

Vector RotateZ( Vector V, _VectorPrecision Angle)
{
    Vector VOut;

    _VectorPrecision CosTheta = cos(Angle);
    _VectorPrecision SinTheta = sin(Angle);

    VOut.x = CosTheta*V.x - SinTheta*V.y;
    VOut.y = SinTheta*V.x + CosTheta*V.y;
    VOut.z = V.z;

    return VOut;
}


Vector RotationAxis( Vector V, Vector Axis, _VectorPrecision Angle)
{
    Vector VOut;

    _VectorPrecision CosTheta = cos(Angle);
    _VectorPrecision SinTheta = sin(Angle);

    VOut.x = (CosTheta + Axis.x*Axis.x*(1.0-CosTheta))*V.x          + (Axis.x*Axis.y*(1.0-CosTheta) - Axis.z*SinTheta)*V.y +  (Axis.x*Axis.z*(1.0-CosTheta) + Axis.y*SinTheta)*V.z;
    VOut.y = (Axis.y*Axis.x*(1.0 - CosTheta) + Axis.z*SinTheta)*V.x + (CosTheta + Axis.y*Axis.y*(1.0 - CosTheta))*V.y      +  (V.y*V.z*(1.0-CosTheta) - Axis.x*SinTheta)*V.z;
    VOut.z = (Axis.z*Axis.x*(1.0-CosTheta)-Axis.y*SinTheta)*V.x      + (Axis.z*Axis.y*(1.0-CosTheta) + Axis.x*SinTheta)*V.y +  (CosTheta + Axis.z*Axis.z*(1.0-CosTheta))*V.z;

    return VOut;
}

bool RayPlaneIntersect(Vector PlaneUnitNormal, _VectorPrecision DistToPlane, //these define the plane
                       Vector RaySource, Vector RayDirection, // these define the ray
                       _VectorPrecision &DistanceFromSource, Vector &Intersection) //outputs
{
    _VectorPrecision LineDotNorm = RayDirection.Dot(PlaneUnitNormal);
    if(LineDotNorm == 0.0)
    {
        return false;
    }

    DistanceFromSource = (DistToPlane - RaySource.Dot(PlaneUnitNormal))/LineDotNorm;
    Intersection = RaySource + DistanceFromSource*RayDirection;

    return true;
}

//maps theta between -pi-->pi and phi between -pi/2 and +pi/2
void GetAlternateThetaPhi(Vector V, _VectorPrecision &Theta, _VectorPrecision &Phi)
{
    Phi = V.GetPhi();
    Theta = V.GetTheta();

    if(V.x < 0.0)
    {
        Theta = -1.0*Theta;
        Phi = Phi - PI;
    }

    return;
}

void CreatePlane( Vector PlaneOrigin, Vector PlaneNormal, Vector &PlaneNormalOut, _VectorPrecision &DistToPlane)
{
    PlaneNormalOut = PlaneNormal;
    
    if(PlaneNormalOut.Magnitude() != 1.0)
    {
        PlaneNormalOut = PlaneNormalOut.Normalized();
    }
    
    if(PlaneOrigin.Dot(PlaneNormalOut) < 0.0)
    {
        //cout << "Input Normal points towards Origin" << endl;
        //need to check like this or 0.0 --> -0.0
        //-0.0 with the trig functions goes to the other end of the output range (0-->2pi) which can cause problems
        if(PlaneNormalOut.x != 0.0)
        {
            PlaneNormalOut.x *= -1.0;
        }
        if(PlaneNormalOut.y != 0.0)
        {
            PlaneNormalOut.y *= -1.0;
        }
        if(PlaneNormalOut.z != 0.0)
        {
            PlaneNormalOut.z *= -1.0;
        }
    
    }
    
    DistToPlane = PlaneNormalOut.Dot(PlaneOrigin);
}

//TODO: Come up with a new name for this
//Function to take a vector from coordinate frame A to B where in A the vector 0,0,1 gets rotated to (theta,phi,1 [spherical coordinates])
//theta and phi describe the new frame, Vin is the input.
Vector TransformToNewFrame( const Vector Vin, const _VectorPrecision theta, const _VectorPrecision phi)
{
    _VectorPrecision CT = cos(theta);
    _VectorPrecision ST = sin(theta);
    _VectorPrecision CP = cos(phi);
    _VectorPrecision SP = sin(phi);
    Vector VOut( CP*CT*Vin.x - SP*Vin.y + CP*ST*Vin.z,
                SP*CT*Vin.x + CP*Vin.y + SP*ST*Vin.z,
                CT*Vin.z - ST*Vin.x);
    
    return VOut;
}

#endif
