#pragma once
#include "Common/ZOperatorOverloading.h"
#include <numeric>

namespace Zee {
namespace Math {
namespace Dbl {
	constexpr double Epsilon = std::numeric_limits<double>::epsilon();
	constexpr double Dbl_Max = std::numeric_limits<double>::max(); 
	constexpr double Dbl_Min = std::numeric_limits<double>::min(); 
	constexpr double PI = 3.141592653589793238;
	constexpr double PI_Half = PI * 0.5;
	constexpr double PI_Quad = PI * 0.25;
	constexpr double PI_Invert = 1.0 / PI;
	constexpr double DegreeToRadian = PI / 180.0;
	constexpr double RadianToDegree = 180.0 * PI_Invert;
}//namespace Dbl

namespace Flt {
	
	constexpr float Epsilon = std::numeric_limits<float>::epsilon();
	constexpr float Flt_Max = std::numeric_limits<float>::max(); 
	constexpr float Flt_Min = std::numeric_limits<float>::min(); 
	
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
	using mat2 = Math::mat2x2_t<float>;
	
	using quat = Math::quat_base<float>;
}//namespace Zee 