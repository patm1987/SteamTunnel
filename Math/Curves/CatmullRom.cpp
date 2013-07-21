#include "CatmullRom.h"

#include <cassert>

const float CatmullRom::kDefaultTightness = .5f;

CatmullRom::CatmullRom(void) :
	m_controlPoints(),
	m_tightness(kDefaultTightness)
{
}


CatmullRom::~CatmullRom(void)
{
}

void CatmullRom::AddControlPoint(const Eigen::Vector3f &controlPoint)
{
	m_controlPoints.push_back(controlPoint);
}

void CatmullRom::SetTightness(float tightness)
{
	assert(tightness >= 0.f && tightness <= 1.f);
	m_tightness = tightness;
}

Eigen::Vector3f CatmullRom::Evaluate(float t)
{
	if (m_controlPoints.size() <= 0)
	{
		return Eigen::Vector3f(0.f, 0.f, 0.f);
	}
	else if (m_controlPoints.size() == 1)
	{
		return m_controlPoints[0];
	}
	else
	{
		// figure out which control points we're between
		float amount = t / float(m_controlPoints.size() - 1);
		int startIndex = int(amount);
		int endIndex = startIndex + 1;
		amount = amount - float(startIndex);

		// TODO: use "near" rather than equals to (floating point error);
		if(amount == 0.f)
		{
			return m_controlPoints[startIndex];
		}
		else if (amount == 1.f)
		{
			return m_controlPoints[endIndex];
		}
		else
		{
			return EvaluateSubset(startIndex, endIndex, amount);
		}
	}
}

Eigen::Vector3f CatmullRom::EvaluateSubset(int p0, int p1, float amount)
{
	assert(p0 >= 0 && p0 < int(m_controlPoints.size()));
	assert(p1 >= 0 && p1 < int(m_controlPoints.size()));
	assert(amount >= 0.f && amount <= 1.f);

	const Eigen::Vector3f &point0 = m_controlPoints[p0];
	const Eigen::Vector3f &point1 = m_controlPoints[p1];
	Eigen::Vector3f tangent0 = ComputeTangent(p0);
	Eigen::Vector3f tangent1 = ComputeTangent(p1);

	float t = amount;
	float t2 = t * t;
	float t3 = t2 * t;
	Eigen::Vector4f hermiteVector(
		2.f * t3 - 3.f * t2 + 1.f,
		-2.f * t3 + 3.f * t2,
		t3 - 2.f * t2 + t,
		t3 - t2);

	return hermiteVector[0] * point0
		+ hermiteVector[1] * point1
		+ hermiteVector[2] * tangent0
		+ hermiteVector[3] * tangent1;
}

Eigen::Vector3f CatmullRom::ComputeTangent(int index)
{
	// todo: leaving
	// have to replace tangent computation with 2* (pi+1 pi) and 2* (pi-1 pi) for edge vertices (not just same vertex twice!)
	// see notes
	if (index <= 0)
	{
		return m_controlPoints[0];
	}
	else if (index >= int(m_controlPoints.size()) - 1)
	{
		return m_controlPoints[m_controlPoints.size() - 1];
	}
	else
	{
		return m_tightness * (m_controlPoints[index + 1] - m_controlPoints[index - 1]);
	}
}
