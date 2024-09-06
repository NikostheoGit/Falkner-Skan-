#include<iostream>
#include<cmath>
#include<ostream>
#include<fstream>
#include<string>
#ifndef MVECTOR_H // the 'include guard'
#define MVECTOR_H // see C++ Primer Sec. 2.9.2

#include <vector>
using namespace std;



// Class that represents a mathematical vector
class MVector
{
public:
	// constructors
	MVector() {}
	explicit MVector(int n) : v(n) {}
	MVector(int n, double x) : v(n, x) {}

	// access element (lvalue)
	double &operator[](int index) { return v[index]; }

	// access element (rvalue)
	double operator[](int index) const { return v[index]; }

	int size() const { return v.size(); } // number of elements
	void push_back(double x) { v.push_back(x); }
private:
	std::vector<double> v;
};

#endif

// Operator overload for "scalar * vector"
inline MVector operator*(const double& lhs, const MVector& rhs)
{
	MVector temp(rhs);
	for (int i = 0; i<temp.size(); i++) temp[i] *= lhs;
	return temp;
}
// Operator overload for vector / scalar
inline MVector operator/(const MVector& rhs, const double& lhs)
{
	if (lhs == 0) cout << "The denominator is 0" << std::endl;

	else {
		MVector temp(rhs);
		for (int i = 0; i < temp.size(); i++) temp[i] /= lhs;
		return temp;
	}
}
// Operator overload for "vector+vector"
inline MVector operator+(const MVector& rhs, const MVector& lhs)
{
	if (rhs.size() != lhs.size()) cout << "These vectors have different size" << std::endl;

	else {
		MVector temp(rhs), temp1(lhs);
		for (int i = 0; i < temp.size(); i++)
			temp[i] += temp1[i];

		return temp;
	}
}

// Operator overload for "vector-vector"
inline MVector operator-(const MVector& rhs, const MVector& lhs)
{
	if (rhs.size() != lhs.size()) cout << "These vectors have different size" << std::endl;

	else {
		MVector temp(rhs), temp1(lhs);
		for (int i = 0; i < temp.size(); i++)
			temp[i] -= temp1[i];
		return temp;
	}
}

// Operator overload for "vector*vector"
double operator*(const MVector& rhs, const MVector& lhs)

{
	if (rhs.size() != lhs.size()) cout << "These vectors have different size" << std::endl;


	else
	{
		MVector temp(rhs), temp2(lhs);
		double sum = 0.0;
		for (int i = 0; i < temp.size(); i++) {
			sum += temp[i] * temp2[i];

		}
		return sum;
	}
}

ostream& operator<<(ostream& os, const MVector& v)
{
	int n = v.size();
	//os << "(";
	for (int i = 0; i < n; i++) {
		// write the code to output the vector components here.
		os << v[i] << " ";
	}

	//os << ")";
	return os;
}
struct MFunction
{
	virtual MVector operator()(const double& x,
		const MVector& y) = 0;
};
class FunctionF1 : public MFunction
{
public:
	MVector operator()(const double& x, const MVector& y)
	{
		MVector temp(2);
		temp[0] = y[0] + x*y[1];
		temp[1] = x*y[0] - y[1];
		return temp;
	}
};

class FunctionF2 : public MFunction
{
public:
	MVector operator()(const double& x, const MVector& y)
	{
		MVector temp(2);
		temp[0] = x;
		temp[1] = y[1];
		return temp;
	}
};

class FunctionF3 : public MFunction
{
public:
	MVector operator()(const double& x, const MVector& y)
	{
		MVector temp(2);
		temp[0] = y[1];
		temp[1] = (32 + 2 * pow(x, 3) - y[0] * y[1]) / 8;
		return temp;
	}
};



// Declaration for an Euler scheme ODE solver

int EulerSolve(int steps, double a, double b, MVector &y, MFunction &f)
{
	ofstream outFile("Euler88.txt");
	double x, h = (b - a) / steps;
	MVector yy;


	for (int i = 0; i <= steps; i++) {

		x = a + i*h;
		yy = y;
		y = yy + h*f(x, y);
		outFile.precision(16);
		outFile.width(10);  outFile << yy << "\t";
		outFile.width(10);  outFile << x << endl;



	}

	cout << yy << endl;
	return 0;
}

// Declaration for Midpoint scheme ODE solver
int Midpoint(int steps, double a, double b, MVector &y, MFunction &f)
{
	ofstream outFile("Midpoint88.txt");
	double x, h = (b - a) / steps;
	MVector yy;

	for (int i = 0; i <= steps; i++) {

		x = a + i*h;
		yy = y;
		y = yy + h*f(x + h / 2, yy + h*f(x, yy) / 2);
		outFile.precision(16);
		outFile.width(10);  outFile << yy << "\t";
		outFile.width(10);  outFile << x << endl;
	}
	cout << yy << endl;
	return 0;


}

//Declaration for Runge Kutta scheme ODE solver

int RungeKutta(int steps, double a, double b, MVector &y, MFunction &f, string filename = "")
{
	{
		// At the beginning of the function...
		bool writeToFile = filename.size()>0;
		std::ofstream myFile;
		if (writeToFile)
		{
			myFile.open(filename.c_str());
			if (!myFile)
			{
				std::cout << "Could not open " << filename << std::endl;
				writeToFile = false;
			}
		}


		double x, h = (b - a) / steps;
		MVector yy, k1, k2, k3, k4;
		if (writeToFile)
		{
			myFile.precision(8);
			myFile.width(10);	myFile << y << "\t";
			myFile.width(10);	myFile << a << endl;
		}


		for (int i = 0; i < steps; i++) {

			x = a + i*h;
			yy = y;
			k1 = f(x, y);
			k2 = f(x + h / 2, y + (h / 2)*k1);
			k3 = f(x + h / 2, y + (h / 2)*k2);
			k4 = f(x + h, y + h*k3);
			y = yy + (h / 6)*(k1 + 2 * k2 + 2 * k3 + k4);

			if (writeToFile) //wherever we want to write something
			{

				myFile.width(10);			myFile << y << "\t";
				myFile.width(10);	myFile << x + h << endl;
			}


		}

		cout.precision(5);

		// At the end of the function...
		if (writeToFile)
		{
			myFile.close();
		}
	}
	return 0;
}



class Eqn1p5Derivs : public MFunction //Class Representing the funcition derivatives (1.2.4)
{
public:
	// constructor to initialise kappa
	Eqn1p5Derivs() { kappa = 1.0; }
	MVector operator()(const double& x, const MVector& y)
	{
		MVector temp(4);
		temp[0] = y[1];
		temp[1] = -kappa*y[1] - x*y[0];
		temp[2] = y[3];
		temp[3] = -kappa*y[3] - x*y[2];
		return temp;
	}
	void SetKappa(double k) { kappa = k; } // change kappa
private:
	double kappa; // class member variable, accessible within
				  // all Eqn1p5Derivs member functions
};

class FunctionF4 : public MFunction
{
public:
	MVector operator()(const double& x, const MVector& y)
	{

		MVector temp(4);
		temp[0] = y[1];
		temp[1] = (32 + 2 * pow(x, 3) - y[0] * y[1]) / 8;
		temp[2] = y[3];
		temp[3] = -(y[0] * y[3]) / 8 - (y[1] * y[2]) / 8;
		return temp;

	}
};

class FunctionF5 : public MFunction // Last page Function
{
public:
	MVector operator()(const double& x, const MVector& y)
	{
		double b = 1.0 / 2.0;
		MVector temp(3);
		temp[0] = y[1];//=y'
		temp[1] = y[2];//=y"
		temp[2] = -y[0] * y[2] - b*(1 - pow(y[1], 2));
		return temp;
	}
};


class FunctionF6 : public MFunction // Last page Function
{
public:
	FunctionF6() { beta = 0.5; }
	MVector operator()(const double& x, const MVector& y)
	{
		MVector temp(6);
		temp[0] = y[1];//=y'
		temp[1] = y[2];//=y"
		temp[2] = -y[0] * y[2] - beta*(1 - pow(y[1], 2));
		temp[3] = y[4]; //z2
		temp[4] = y[5];  //z3
		temp[5] = -y[2] * y[3] + 2 * beta*y[1] * y[4] - y[0] * y[5];

		return temp;
	}
	void SetBeta(double b) { beta = b; } // change kappa
private:
	double beta;
};

class FunctionF7 : public MFunction // 
{
public:

	MVector operator()(const double& x, const MVector& y)
	{
		double b; //try for b=0,1/2,1
		MVector temp(6);
		temp[0] = y[1];//=y'
		temp[1] = y[2];//=y"
		temp[2] = -y[0] * y[2] - b*(1 - pow(y[1], 2));
		temp[3] = y[4]; //z2
		temp[4] = y[5];  //z3
		temp[5] = -y[2] * y[3] + 2 * b*y[1] * y[4] - y[0] * y[5];

		return temp;
	}
};

double FalknerSkan(double b, double k)
{


	FunctionF6 f; //last Shit
	f.SetBeta(b);
	int maxNewtonSteps = 100;
	double tol = pow(10.0, -8);
	for (int i = 0; i < maxNewtonSteps; i++)
	{

		MVector y(6);

		y[0] = 0; y[1] = 0; y[2] = k; y[3] = 0; y[4] = 0; y[5] = 1.0;
		RungeKutta(1000, 0, 50, y, f, "poutsa.txt"); // solve IVP
		double phi = y[1] - 1.0; // calculate residual
		double phidash = y[4]; // 'Jacobian' phidash = z_1(x=1)
		cout << i << endl;
		if (std::abs(phi) <= tol) break; // exit if converged/		
		k -= phi / phidash;
		//if (std::abs(phi) > tol){
		//	maxNewtonSteps = maxNewtonSteps + 1;
	}
	//cout << k << "   " << b << endl;
	return k;
}



int main()
{



	//double h = 0.1, x = 0.5;

	//	MVector u;
	//	u.push_back(1);
	//	u.push_back(2);




	//	cout << "u=" << u << endl;


	// MVector v, y(2); // initialise y with 2 elements
	//		FunctionF1 f; // f has order 2 by definition
	//		y[0] = 1.4; y[1] = -5.7; // assign element values in y
	//		v = f(2., y); // evaluate function f as required

	//		std::cout << "v=" << v << "; y=" << y << std::endl;
	//		v = u + f(2.0, y);
	//		std::cout << "v=" << v << "; y=" << y << std::endl;

	//		v = u + h*f(x,u+h*y);

	//		std::cout << "v=" << v << "; y=" << y << std::endl;


//		MVector v, y(2); // initialise y with 2 elements
//		FunctionF2 f;

//			y[0] = 0; y[1] = 1;


//	EulerSolve(1000, 0, 1, y, f); // If f and y where diferent numbers then we would have different size ventors to multiplay and so we would have an error
//	y[0] = 0; y[1] = 1;
//			Midpoint(1000, 0, 1, y, f);
//			y[0] = 0; y[1] = 1;
//		RungeKutta(1000, 0, 1, y, f,"rk88.txt");



//	MVector v, y(2); // initialise y with 2 elements
//	FunctionF3 f;

//	y[0] = 17; y[1] = 1;


//	RungeKutta(1000, 1.0, 3.0, y, f, "Rk88.txt");
//	y[0] = 17; y[1] = 1;
//	EulerSolve(1000, 1.0, 3.0, y, f); // If f and y where diferent numbers then we would have different size ventors to multiplay and so we would have an error
//	y[0] = 17; y[1] = 1;
//	Midpoint(1000, 1.0, 3.0, y, f);




	//	Eqn1p5Derivs f;
	//	int maxNewtonSteps = 100;
	//	double guess = 0;
	//	double tol = 1e-8;
	//	for (int i = 0; i < maxNewtonSteps; i++)
	//	{
	//		MVector y(4);
	//		// y[0] = y, y[1] = dy/dx, y[2] = Z_1, y[3] = Z_2
	//		y[0] = 0; y[1] = guess; y[2] = 0.0; y[3] = 1.0;
	//		RungeKutta(100, 0.0, 1.0, y, f); // solve IVP
	//		double phi = y[0] - 1; // calculate residual
	//		double phidash = y[2]; // 'Jacobian' phidash = z_1(x=1)
	//		if (std::abs(phi) < tol) break; // exit if converged
	//		guess -= phi / phidash; // apply newton step
	//	}


	//	return 0;


	//	FunctionF4 f;
	//	int maxNewtonSteps = 1;
	//		double guess = 0;
	//		double tol = pow(10.0,-8);
	//		for (int i = 0; i < maxNewtonSteps; i++)
	//		{

	//			MVector y(4);
	//		//y[0] = y, y[1] = dy / dx, y[2] = Z_1, y[3] = Z_2
	//			y[0] = 17.0; y[1] = guess; y[2] = 0.0; y[3] = 1.0;

	//            RungeKutta(100, 1.0, 3.0, y, f); // solve IVP
	//			double phi = y[0] - 43.0/3.0; // calculate residual
	//			double phidash = y[2]; // 'Jacobian' phidash = z_1(x=1)
	//			if (std::abs(phi) <= tol) break; // exit if converged
	//            guess -= phi / phidash;
	//			if (std::abs(phi) > tol){
	//				maxNewtonSteps = maxNewtonSteps + 1;
	//}





	//		}


	//	MVector u,w; //VECTOR*VECTOR TEST
	//	u.push_back(2);
	//	u.push_back(2);
	//	u.push_back(2);
	//	w.push_back(3);
	//   w.push_back(4);
	//
	//	cout << u << w << endl;
	//	cout << u+w << endl; 

	//   MVector v, y(3); // initialise y with 2 elements
	//	FunctionF5 f;
	//	double guess = 0.92;
	//	y[0] = 0; y[1] = 0; y[2] = guess;
	//   RungeKutta(8000, 0, 100, y, f, "poutsa2.txt");

	//	MVector v, y(6);
	//	FunctionF6 f;
	//   double guess = 0.92;

	//	y[0] = 0; y[1] = 0; y[2] = guess; y[3] = 0; y[4] = 0; y[5] = 1;
	//	RungeKutta(8000, 0, 100, y, f, "poutsa1.txt");

	//	FunctionF6 f; //last Shit
	//	int maxNewtonSteps = 100;
	//	double guess = 0.92;
	//	double tol = pow(10.0, -8);
	//	for (int i = 0; i < maxNewtonSteps; i++)
	//	{

	//	MVector y(6);

	//		y[0] = 0; y[1] = 0; y[2] = guess; y[3] = 0; y[4] = 0; y[5] = 1.0;

	//		RungeKutta(8000, 0, 50, y, f, "poutsa.txt"); // solve IVP
	//		double phi = y[1] - 1.0; // calculate residual
	//		double phidash = y[4]; // 'Jacobian' phidash = z_1(x=1)
	//		cout << i << endl;
	//		if (std::abs(phi) <= tol) break; // exit if converged/		
	//	guess -= phi / phidash;
	//		//if (std::abs(phi) > tol){
	//		//	maxNewtonSteps = maxNewtonSteps + 1;
	//		}
	FalknerSkan(0.5, 0.92);
	ofstream outFile("Falkner88.txt");

	double step = 0.0001;
	double out;
	double Equation = FalknerSkan(0.5, 0.92);
	for (double b = 0.000; b <= 1.0; b = b + step)
	{
		out = FalknerSkan(b, Equation);
		Equation = out;
		cout << out << "   " << b << endl;
		outFile.precision(10);
		outFile.width(10);  outFile << out << "\t";
		outFile.width(10);  outFile << b << endl;
	}

	return 0;
}