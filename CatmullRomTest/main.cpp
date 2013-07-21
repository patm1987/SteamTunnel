#include <iostream>

#include "Curves/CatmullRom.h"

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
	cout << "Catmull Rom Test" << endl;

	CatmullRom curve;
	curve.AddControlPoint(Vector3f(0.f, 0.f, 0.f));
	curve.AddControlPoint(Vector3f(0.f, 1.f, 0.f));

	cout << "Some Samples" << endl;
	cout << curve.Evaluate(0.f) << endl << endl;
	cout << curve.Evaluate(.01f) << endl << endl;
	cout << curve.Evaluate(.25f) << endl << endl;
	cout << curve.Evaluate(.5f) << endl << endl;
	cout << curve.Evaluate(.75f) << endl << endl;
	cout << curve.Evaluate(.99f) << endl << endl;
	cout << curve.Evaluate(1.f) << endl << endl;

	return 0;
}
