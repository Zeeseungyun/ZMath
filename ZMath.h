#pragma once
#include "Common/ZOperatorOverloading.h"

namespace Zee {
namespace Math {
namespace Dbl {
	constexpr double Epsilon = 2.2204460492503131e-015;//*DBL_EPSILON is 2.2204460492503131e-016
	constexpr double Dbl_Max = 1.7976931348623158e+308; //DBL_MAX;
	constexpr double Dbl_Min = 2.2250738585072014e-308; //DBL_MIN;

	constexpr double PI = 3.141592653589793238;
	constexpr double PI_Half = PI * 0.5;
	constexpr double PI_Quad = PI * 0.25;
	constexpr double PI_Invert = 1.0f / PI;
	constexpr double DegreeToRadian = PI / 180.0;
	constexpr double RadianToDegree = 180.0 * PI_Invert;
}//namespace Dbl

namespace Flt {
	constexpr float Epsilon = 1.192092896e-06F;//*FLT_EPSILON is 1.192092896e-07F.
	constexpr float Flt_Max = 3.402823466e+38F; //FLT_MAX
	constexpr float Flt_Min = 1.175494351e-38F; //FLT_MIN

	constexpr float PI = 3.14159265f;
	constexpr float PI_Half = PI * 0.5f;
	constexpr float PI_Quad = PI * 0.25f;
	constexpr float PI_Invert = 1.0f / PI;
	constexpr float DegreeToRadian = PI / 180.0f;
	constexpr float RadianToDegree = 180.0f * PI_Invert;
}//namespace Flt 
	using namespace Flt;

	template<typename T>
	T degreeToRadian(const T& v) {
		return v * DegreeToRadian;
	}

	template<typename T>
	T radianToDegree(const T& v) {
		return v * RadianToDegree;
	}

}//namespace Math 
	using vec4 = Math::vec4_t<float>;
	using vec3 = Math::vec3_t<float>;
	using vec2 = Math::vec2_t<float>;
	
	using mat4 = Math::mat4x4_t<float>;
	using mat3 = Math::mat3x3_t<float>;
	using mat2 = Math::mat3x3_t<float>;
	
	using quat = Math::quat_base<float>;
}//namespace Zee 