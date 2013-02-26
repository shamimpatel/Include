#ifndef _VECTOR_H_
#define _VECTOR_H_

#include <cmath>
#include <iostream>
#include <stdio.h>
using namespace std;


#define PI (3.14159265)

float Deg2Rad( float Degrees )
{
    return Degrees * (PI/180.0f);
}

float Rad2Deg( float Radians )
{
    return Radians * (180.0f/PI);
}


float Rad2Deg2( float Radians )
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

float Deg2Rad2( float Degrees )
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
    double x,y,z;

    Vector ();
    Vector ( double X, double Y, double Z);
    Vector ( double r, double theta, double phi, bool bSpherical)
    {
        x = r*sin(theta)*cos(phi);
        y = r*sin(theta)*sin(phi);
        z = r*cos(theta);
    }

    // slightly more efficient way to create a spherical unit vector.
    //Need the extra input so that there's no compiler error.
    //(float,float,bool) is apparently equivalent to (float,float,float)
    //might be a GNU thing?
    Vector ( double theta, double phi, bool bSpherical, bool bUnit)
    {
        x = sin(theta)*cos(phi);
        y = sin(theta)*sin(phi);
        z = cos(theta);
    }

    float GetTheta()
    {
        return acos(z/Magnitude());
    }

    void SetTheta( float theta )
    {
        float r = Magnitude();
        float phi = GetPhi();

        x = r*sin(theta)*cos(phi);
        y = r*sin(theta)*sin(phi);
        z = r*cos(theta);
    }

    float GetPhi()
    {
        return atan2(y,x);
    }

    void SetPhi( float phi)
    {
        float r = Magnitude();
        float theta = GetTheta();

        x = r*sin(theta)*cos(phi);
        y = r*sin(theta)*sin(phi);
        z = r*cos(theta);
    }

    double Dot( const Vector V )
    {
        return x*V.x + y*V.y + z*V.z;
    }

    Vector Cross ( const Vector V)
    {
        return Vector( y*V.z - z*V.y , z*V.x - x*V.z , x*V.y - y*V.x);
    }

    double Magnitude()
    {
        return sqrt(x*x + y*y + z*z);
    }

    Vector Normalized()
    {
        double Mag = this->Magnitude();
        if(Mag == 0.0)
        {
            return Vector(0,0,0); //avoid the divide by 0 error
        }

        return this->Scale(double(1.0)/Mag);
    }

    Vector Scale( const double Scale )
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

    Vector operator* (double Scale)
    {
        return this->Scale( Scale );
    }


    Vector Divide( double Denominator )
    {
        Vector Out;
        Out.x = x/Denominator;
        Out.y = y/Denominator;
        Out.z = z/Denominator;
        return Out;
    }

    Vector operator/ (double Denominator)
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

Vector operator*(double Scale, Vector V)
{
    return V*Scale;
}

Vector::Vector(double x, double y, double z)
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
Vector ScatterUnitVector( Vector In, double Theta, double Phi)
{
    Vector Out;

    double CT = cos(Theta);
    double ST = sin(Theta);
    double CP = cos(Phi);
    double SP = sin(Phi);

    //http://en.wikipedia.org/wiki/Monte_Carlo_method_for_photon_transport#Step_3:_Absorption_and_scattering
    //And: NUCLEAR SCIENCE AND ENGINEERING: 131, 132–136 ~1999 Direction Cosines and Polarization Vectors
    //for Monte Carlo Photon Scattering

    double zz = sqrt(1.0-In.z*In.z);

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

Vector RotateX( Vector V, double Angle)
{
    Vector VOut;

    float CosTheta = cos(Angle);
    float SinTheta = sin(Angle);

    VOut.x = V.x;
    VOut.y = CosTheta*V.y - SinTheta*V.z;
    VOut.z = SinTheta*V.y + CosTheta*V.z;

    return VOut;
}

Vector RotateY( Vector V, double Angle)
{
    Vector VOut;

    float CosTheta = cos(Angle);
    float SinTheta = sin(Angle);

    VOut.x = CosTheta*V.x + SinTheta*V.z;
    VOut.y = V.y;
    VOut.z = CosTheta*V.z - SinTheta*V.x;

    return VOut;
}

Vector RotateZ( Vector V, double Angle)
{
    Vector VOut;

    double CosTheta = cos(Angle);
    double SinTheta = sin(Angle);

    VOut.x = CosTheta*V.x - SinTheta*V.y;
    VOut.y = SinTheta*V.x + CosTheta*V.y;
    VOut.z = V.z;

    return VOut;
}


Vector RotationAxis( Vector V, Vector Axis, double Angle)
{
    Vector VOut;

    double CosTheta = cos(Angle);
    double SinTheta = sin(Angle);

    VOut.x = (CosTheta + Axis.x*Axis.x*(1.0-CosTheta))*V.x          + (Axis.x*Axis.y*(1.0-CosTheta) - Axis.z*SinTheta)*V.y +  (Axis.x*Axis.z*(1.0-CosTheta) + Axis.y*SinTheta)*V.z;
    VOut.y = (Axis.y*Axis.x*(1.0 - CosTheta) + Axis.z*SinTheta)*V.x + (CosTheta + Axis.y*Axis.y*(1.0 - CosTheta))*V.y      +  (V.y*V.z*(1.0-CosTheta) - Axis.x*SinTheta)*V.z;
    VOut.z = (Axis.z*Axis.x*(1.0-CosTheta)-Axis.y*SinTheta)*V.x      + (Axis.z*Axis.y*(1.0-CosTheta) + Axis.x*SinTheta)*V.y +  (CosTheta + Axis.z*Axis.z*(1.0-CosTheta))*V.z;

    return VOut;
}

bool RayPlaneIntersect(Vector PlaneUnitNormal, double DistToPlane,
                       Vector RaySource, Vector RayDirection,
                       double &DistanceFromSource, Vector &Intersection)
{
    float LineDotNorm = RayDirection.Dot(PlaneUnitNormal);
    if(LineDotNorm == 0.0)
    {
        return false;
    }

    DistanceFromSource = (DistToPlane - RaySource.Dot(PlaneUnitNormal))/LineDotNorm;
    Intersection = RaySource + DistanceFromSource*RayDirection;

    return true;
}

//maps theta between -pi-->pi and phi between -pi/2 and +pi/2
void GetAlternateThetaPhi(Vector V, double &Theta, double &Phi)
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


#endif
