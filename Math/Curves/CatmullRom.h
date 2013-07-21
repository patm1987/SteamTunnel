#ifndef _CATMULL_ROM_H
#define _CATMULL_ROM_H

#include <vector>

#include "Eigen/Dense"

class CatmullRom
{
private:
	static const float kDefaultTightness;

public:
	CatmullRom(void);
	~CatmullRom(void);

	void AddControlPoint(const Eigen::Vector3f &controlPoint);
	void SetTightness(float tightness);

	Eigen::Vector3f Evaluate(float t);

private:
	Eigen::Vector3f EvaluateSubset(int p0, int p1, float amount);
	Eigen::Vector3f ComputeTangent(int index);

	std::vector<Eigen::Vector3f> m_controlPoints;
	float m_tightness;
};

#endif//_CATMULL_ROM_H
