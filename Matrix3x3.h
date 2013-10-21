//
//  Matrix3x3.h
//  TestCode
//
//  Created by Shamim Patel on 17/10/2013.
//  Copyright (c) 2013 Shamim Patel. All rights reserved.
//

#ifndef _MATRIX3X3_
#define _MATRIX3X3_


class Matrix3x3
{
public:
	Matrix3x3()
	{
		Data[0] = 0.0;
		Data[1] = 0.0;
		Data[2] = 0.0;
		Data[3] = 0.0;
		Data[4] = 0.0;
		Data[5] = 0.0;
		Data[6] = 0.0;
		Data[7] = 0.0;
		Data[8] = 0.0;
	}
	
	
	Matrix3x3( double m0, double m1, double m2,
			   double m3, double m4, double m5,
			   double m6, double m7, double m8)
	{
		Data[0] = m0; Data[1] = m1; Data[2] = m2;
		Data[3] = m3; Data[4] = m4; Data[5] = m5;
		Data[6] = m6; Data[7] = m7; Data[8] = m8;
	}
	
	
	double Data[9];
	
	
	Matrix3x3 operator*(const Matrix3x3 Other) const
	{
		Matrix3x3 Out;
		Out.Data[0] = Data[0]*Other.Data[0] + Data[1]*Other.Data[3] + Data[2]*Other.Data[6];
		Out.Data[1] = Data[0]*Other.Data[1] + Data[1]*Other.Data[4] + Data[2]*Other.Data[7];
		Out.Data[2] = Data[0]*Other.Data[2] + Data[1]*Other.Data[5] + Data[2]*Other.Data[8];
		
		Out.Data[3] = Data[3]*Other.Data[0] + Data[4]*Other.Data[3] + Data[5]*Other.Data[6];
		Out.Data[4] = Data[3]*Other.Data[1] + Data[4]*Other.Data[4] + Data[5]*Other.Data[7];
		Out.Data[5] = Data[3]*Other.Data[2] + Data[4]*Other.Data[5] + Data[5]*Other.Data[8];
		
		Out.Data[6] = Data[6]*Other.Data[0] + Data[7]*Other.Data[3] + Data[8]*Other.Data[6];
		Out.Data[7] = Data[6]*Other.Data[1] + Data[7]*Other.Data[4] + Data[8]*Other.Data[7];
		Out.Data[8] = Data[6]*Other.Data[2] + Data[7]*Other.Data[5] + Data[8]*Other.Data[8];
		return Out;
	}
	
	Matrix3x3 GetTransposed() const
	{
		Matrix3x3 Out;
		
		Out.Data[0] = Data[0]; Out.Data[1] = Data[3]; Out.Data[2] = Data[6];
		Out.Data[3] = Data[1]; Out.Data[4] = Data[4]; Out.Data[5] = Data[7];
		Out.Data[6] = Data[2]; Out.Data[7] = Data[5]; Out.Data[8] = Data[8];
		
		return Out;
	}
	
	void Print()
	{
		cout << Data[0] << '\t' << Data[1] << '\t' << Data[2] << '\t' << endl;
		cout << Data[3] << '\t' << Data[4] << '\t' << Data[5] << '\t' << endl;
		cout << Data[6] << '\t' << Data[7] << '\t' << Data[8] << '\t' << endl;
	}
	
};


#endif //_MATRIX3X3_
