#pragma once
namespace Zee {
namespace Math {
	template<int Row, int Col> struct bools;

	template<typename ScalarType, int Dim> struct vec_base;
	template<typename ScalarType, int Row, int Col> struct mat_base;
	template<typename ScalarType> struct quat_base;
	template<typename ScalarType> struct complex_base;

	template<typename ScalarType> using vec2_t = vec_base<ScalarType, 2>;
	template<typename ScalarType> using vec3_t = vec_base<ScalarType, 3>;
	template<typename ScalarType> using vec4_t = vec_base<ScalarType, 4>;

	template<typename ScalarType> using mat2x2_t = mat_base<ScalarType, 2, 2>;
	template<typename ScalarType> using mat3x3_t = mat_base<ScalarType, 3, 3>;
	template<typename ScalarType> using mat4x4_t = mat_base<ScalarType, 4, 4>;

	template<typename ScalarType> using quat_t    = quat_base<ScalarType>;
	template<typename ScalarType> using complex_t = complex_base<ScalarType>;
}//namespace Math
}//namespace Zee