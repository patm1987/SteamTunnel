#include "CatmullRom.h"

#include <cassert>

//! @brief	how tight is the curve by default
const float CatmullRom::kDefaultTightness = .5f;

CatmullRom::CatmullRom(void) :
	m_controlPoints(),
	m_tightness(kDefaultTightness)
{
}


CatmullRom::~CatmullRom(void)
{
}

/*!
 * @brief	adds the given control point to the spline
 * @param	controlPoint the control point to add
 */
void CatmullRom::AddControlPoint(const Eigen::Vector3f &controlPoint)
{
	m_controlPoints.push_back(controlPoint);
}

/*!
 * @brief	configures the @a tightness of the catmull rom hermite spline
 */
void CatmullRom::SetTightness(float tightness)
{
	assert(tightness >= 0.f && tightness <= 1.f);
	m_tightness = tightness;
}

/*!
 * @brief	evaluates this curve at @a t
 * @param	t the amoun along the curve to evaluate. t e [0, count - 1]
 */
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

/*!
 * @brief	evaluates the curve between two control points
 * @param	p0 the start of the curve subset
 * @param	p1 the end of the curve subset
 * @param	amount a value between 0 and 1 describing how far along p0 and p1
 *			we need to evaluate
 */
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

/*!
 * @breif	compute the tangent for the given @a index
 * @param	index the index for which we'll compute the tangent
 * @note	to attempt to maintain a smooth curve, endpoints are currently
 *			computed as p_0 = 2* (p_1 - p_0) and p_i = 2*(p_i - p_i-1)
 */
Eigen::Vector3f CatmullRom::ComputeTangent(int index)
{
	if (index <= 0)
	{
		return 2.f * (m_controlPoints[1] - m_controlPoints[0]);
	}
	else if (index >= int(m_controlPoints.size()) - 1)
	{
		int lastIndex = m_controlPoints.size() - 1;
		return 2.f * (m_controlPoints[lastIndex] - m_controlPoints[lastIndex - 1]);
	}
	else
	{
		return m_tightness * (m_controlPoints[index + 1] - m_controlPoints[index - 1]);
	}
}
