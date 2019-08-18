#pragma once
#include "ZTypeTraits.h"
#include <cmath>
#include <cassert>
#ifdef Z_MATH_STATIC_ASSERT
#define Z_MATH_ASSERT(exp, explain) static_assert((exp),explain)

#else
#define Z_MATH_ASSERT(exp, explain) assert((exp)&&explain)
#endif

#define Z_MATH_ASSERT_EQUIVALENT Z_MATH_ASSERT(0, __FUNCTION__ " : must equivalent type.")
#define Z_MATH_ASSERT_SQUARE Z_MATH_ASSERT(0, __FUNCTION__ " : matrix must square type.")

namespace Zee {
namespace Math {
	template<typename T, typename U>
	T& assign(T& dst, const U& src) {
		static_assert(is_arithmetic<T>::value && is_arithmetic<U>::value, __FUNCTION__ " : must arithmetic type.");
		return dst = static_cast<T>(src);
	}

	template<typename T, typename U>
	bool is_equal(const T& lhs, const U& rhs) {
		static_assert(is_arithmetic<T>::value && is_arithmetic<U>::value, __FUNCTION__ " : must arithmetic type.");
		using PromotionType = promotion_t<T, U>;
		return static_cast<PromotionType>(lhs) == static_cast<PromotionType>(rhs);
	}

	template<typename T, typename U>
	bool is_not_equal(const T& lhs, const U& rhs) {
		static_assert(is_arithmetic<T>::value && is_arithmetic<U>::value, __FUNCTION__ " : must arithmetic type.");
		using PromotionType = promotion_t<T, U>;
		return static_cast<PromotionType>(lhs) != static_cast<PromotionType>(rhs);
	}

	template<typename T, typename U , typename EpsT>
	bool is_near_equal(const T& lhs, const U& rhs , const EpsT& e) {
		static_assert(is_arithmetic<T>::value && is_arithmetic<U>::value, __FUNCTION__ " : must arithmetic type.");
		static_assert(is_floating_point<EpsT>::value, __FUNCTION__ "EpsT must floating point type.");
		return std::abs(static_cast<EpsT>(lhs) - static_cast<EpsT>(rhs)) < e;
	}

	template<typename T, typename U, typename EpsT>
	bool is_not_near_equal(const T& lhs, const U& rhs, const EpsT& e) {
		return !is_near_equal(lhs,rhs, e);
	}

	template<typename T>
	bool is_zero(const T& v) {
		static_assert(is_arithmetic<T>::value, __FUNCTION__ " : must arithmetic type.");
		return v == T(0);
	}

	template<typename T>
	bool is_not_zero(const T& v) {
		static_assert(is_arithmetic<T>::value, __FUNCTION__ " : must arithmetic type.");
		return v != T(0);
	}

	template<typename T, typename EpsT>
	bool is_near_zero(const T& v , const EpsT& e) {
		static_assert(is_arithmetic<T>::value, __FUNCTION__ " : must arithmetic type.");
		static_assert(is_floating_point<EpsT>::value, __FUNCTION__ "EpsT must floating point type.");
		using PromotionType = promotion_t<T, EpsT>;
		return static_cast<PromotionType>(std::abs(v)) < e;
	}

	template<typename T, typename EpsT>
	bool is_not_near_zero(const T& v, const EpsT& e) {
		return !is_near_zero(v, e);
	}

	template<typename T, typename U>
	bool is_less(const T& lhs, const U& rhs) {
		static_assert(is_arithmetic<T>::value && is_arithmetic<U>::value, __FUNCTION__ " : must arithmetic type.");
		using PromotionType = promotion_t<T, U>;
		return static_cast<PromotionType>(lhs) < static_cast<PromotionType>(rhs);
	}

	template<typename T, typename U>
	bool is_less_equal(const T& lhs, const U& rhs) {
		static_assert(is_arithmetic<T>::value && is_arithmetic<U>::value, __FUNCTION__ " : must arithmetic type.");
		using PromotionType = promotion_t<T, U>;
		return static_cast<PromotionType>(lhs) <= static_cast<PromotionType>(rhs);
	}

	template<typename T, typename U>
	bool is_greater(const T& lhs, const U& rhs) {
		static_assert(is_arithmetic<T>::value && is_arithmetic<U>::value, __FUNCTION__ " : must arithmetic type.");
		using PromotionType = promotion_t<T, U>;
		return static_cast<PromotionType>(lhs) > static_cast<PromotionType>(rhs);
	}

	template<typename T, typename U>
	bool is_greater_equal(const T& lhs, const U& rhs) {
		static_assert(is_arithmetic<T>::value && is_arithmetic<U>::value, __FUNCTION__ " : must arithmetic type.");
		using PromotionType = promotion_t<T, U>;
		return static_cast<PromotionType>(lhs) >= static_cast<PromotionType>(rhs);
	}

	template<typename RetT, typename T>
	RetT abs(const T& v) {
		static_assert(is_arithmetic<T>::value, __FUNCTION__ " : must arithmetic type.");
		static_assert(is_arithmetic<RetT>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return static_cast<RetT>(std::abs(v));
	}

	template<typename RetT, typename T>
	RetT minus(const T& v) {
		static_assert(is_arithmetic<T>::value, __FUNCTION__ " : must arithmetic type.");
		static_assert(is_arithmetic<RetT>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return static_cast<RetT>(-v);
	}

	template<typename RetT, typename T>
	RetT plus(const T& v) {
		static_assert(is_arithmetic<T>::value, __FUNCTION__ " : must arithmetic type.");
		static_assert(is_arithmetic<RetT>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return static_cast<RetT>(v);
	}

	template<typename RetT, typename T, typename U>
	RetT add(const T& lhs, const U& rhs) {
		static_assert(is_arithmetic<T>::value && is_arithmetic<U>::value, __FUNCTION__ " : must arithmetic type.");
		static_assert(is_arithmetic<RetT>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return static_cast<RetT>(static_cast<promotion_t<RetT, T, U>>(lhs) + static_cast<promotion_t<RetT, T, U>>(rhs));
	}

	template<typename RetT, typename T, typename U>
	RetT sub(const T& lhs, const U& rhs) {
		static_assert(is_arithmetic<T>::value && is_arithmetic<U>::value, __FUNCTION__ " : must arithmetic type.");
		static_assert(is_arithmetic<RetT>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return static_cast<RetT>(static_cast<promotion_t<RetT, T, U>>(lhs) - static_cast<promotion_t<RetT, T, U>>(rhs));
	}

	template<typename RetT, typename T, typename U>
	RetT mul(const T& lhs, const U& rhs) {
		static_assert(is_arithmetic<T>::value && is_arithmetic<U>::value, __FUNCTION__ " : must arithmetic type.");
		static_assert(is_arithmetic<RetT>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return static_cast<RetT>(static_cast<promotion_t<RetT, T, U>>(lhs) * static_cast<promotion_t<RetT, T, U>>(rhs));
	}

	template<typename RetT, typename T, typename U>
	RetT div(const T& lhs, const U& rhs) {
		static_assert(is_arithmetic<T>::value && is_arithmetic<U>::value, __FUNCTION__ " : must arithmetic type.");
		static_assert(is_arithmetic<RetT>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return static_cast<RetT>(static_cast<promotion_t<RetT, T, U>>(lhs) / static_cast<promotion_t<RetT, T, U>>(rhs));
	}

	template<typename RetT, typename T, typename U>
	RetT min(const T& v, const U& range) {
		static_assert(is_arithmetic<T>::value && is_arithmetic<U>::value, __FUNCTION__ " : must arithmetic type.");
		static_assert(is_arithmetic<RetT>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return static_cast<RetT>(v) < static_cast<RetT>(range) ? static_cast<RetT>(v) : static_cast<RetT>(range);
	}

	template<typename RetT, typename T, typename U>
	RetT max(const T& v, const U& range) {
		static_assert(is_arithmetic<T>::value && is_arithmetic<U>::value, __FUNCTION__ " : must arithmetic type.");
		static_assert(is_arithmetic<RetT>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return static_cast<RetT>(v) > static_cast<RetT>(range) ? static_cast<RetT>(v) : static_cast<RetT>(range);
	}

	template<typename RetT, typename T, typename U>
	RetT clamp(const T& v, const U& begin, const U& end) {
		static_assert(is_arithmetic<T>::value && is_arithmetic<U>::value, __FUNCTION__ " : must arithmetic type.");
		static_assert(is_arithmetic<RetT>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return
			static_cast<RetT>(v) < static_cast<RetT>(begin) ?
			static_cast<RetT>(begin) :
			static_cast<RetT>(v) > static_cast<RetT>(end) ?
			static_cast<RetT>(end) :
			static_cast<RetT>(v);
	}

	template<typename RetT, typename T>
	RetT saturate(const T& v) {
		static_assert(is_arithmetic<T>::value, __FUNCTION__ " : must arithmetic type.");
		static_assert(is_arithmetic<RetT>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return
			static_cast<RetT>(v) < static_cast<RetT>(0) ?
			static_cast<RetT>(0) :
			static_cast<RetT>(v) > static_cast<RetT>(1) ?
			static_cast<RetT>(1) :
			static_cast<RetT>(v);
	}

	template<typename RetT, typename T, typename U, typename FloatT>
	RetT lerp(const T& v0, const U& v1, const FloatT& t) {
		static_assert(is_arithmetic<T>::value && is_arithmetic<U>::value, __FUNCTION__ " : must arithmetic type.");
		static_assert(is_arithmetic<RetT>::value, __FUNCTION__ " : Return type must arithmetic type.");
		static_assert(is_floating_point<FloatT>::value, __FUNCTION__ "FloatT must floating point type.");
		using PromotionType = promotion_t<RetT, FloatT, T, U>;

		//return v0 * (1 - t) + v1 * t;
		return static_cast<RetT>(
			static_cast<PromotionType>(v0) * (static_cast<PromotionType>(1) - static_cast<PromotionType>(t)) + 
			static_cast<PromotionType>(v1) * static_cast<PromotionType>(t)
		);
	}

namespace Details {
	template<typename T>
	void hermite_support(const T& t, T(&exp)[4]) {
		const auto t2 = t * t;
		const auto t3 = t2 * t;
		const auto a0 = 2 * t3;
		const auto a1 = 3 * t2;

		exp[0] = a0 - a1 + 1;
		exp[1] = t3 - 2 * t2 + t;
		exp[2] = a1 - a0;
		exp[3] = t3 - t2;
	}

	template<typename RetT, typename T , typename U>
	RetT hermite_calculate(const T(&exp)[4], const U& v0, const U& tan0, const U& v1, const U& tan1) {
		return static_cast<RetT>(
			static_cast<T>(v0) * exp[0] +
			static_cast<T>(tan0) * exp[1] +
			static_cast<T>(v1) * exp[2] +
			static_cast<T>(tan1) * exp[3]
		);
	}

	template<typename T>
	void catmullRom_support(const T& t, T(&exp)[4]) {
		const auto t2 = t * t;
		const auto t3 = t2 * t;
		const auto a0 = t3 * 3;

		exp[0] = -t3 + 2 * t2 - t;
		exp[1] = a0 - 5 * t2 + 2;
		exp[2] = -a0 + 4 * t2 + t;
		exp[3] = t3 - t2;
	}

	template<typename RetT, typename T, typename U>
	RetT catmullRom_calculate(const T(&exp)[4], const U& v0, const U& v1, const U& v2, const U& v3) {
		return static_cast<T>((
			static_cast<T>(v0) * exp[0] +
			static_cast<T>(v1) * exp[1] +
			static_cast<T>(v2) * exp[2] +
			static_cast<T>(v3) * exp[3]
		) * static_cast<T>(0.5));
	}
}//namespace Details

	template<typename RetT, typename T, typename FloatT>
	RetT hermite(const T& v0, const T& tan0, const T& v1, const T& tan1, const FloatT& t) {
		static_assert(is_arithmetic<T>::value, __FUNCTION__ " : must arithmetic type.");
		static_assert(is_arithmetic<RetT>::value, __FUNCTION__ " : Return type must arithmetic type.");
		static_assert(is_floating_point<FloatT>::value, __FUNCTION__ "FloatT must floating point type.");

		using PromotionType = promotion_t<RetT, FloatT, T>;

		PromotionType exp[4];
		Details::hermite_support(static_cast<PromotionType>(t), exp);
		return Details::hermite_calculate<RetT, PromotionType, T>(exp, v0, tan0, v1, tan1);
	}

	//ret : v0 ~ v1 (t)
	template<typename RetT, typename T, typename FloatT>
	RetT hermite(const T& v0, const T& v1, const T& v2, const FloatT& t) {
		static_assert(is_arithmetic<T>::value, __FUNCTION__ " : must arithmetic type.");
		static_assert(is_arithmetic<RetT>::value, __FUNCTION__ " : Return type must arithmetic type.");
		static_assert(is_floating_point<FloatT>::value, __FUNCTION__ "FloatT must floating point type.");
		return hermite<RetT, T, FloatT>(v0, v1 - v0, v1, v2 - v1, t);
	}

	//ret : v1 ~ v2 (t)
	template<typename RetT, typename T, typename FloatT>
	RetT catmullRom(const T& v0, const T& v1, const T& v2, const T& v3 ,const FloatT& t) {
		static_assert(is_arithmetic<T>::value, __FUNCTION__ " : must arithmetic type.");
		static_assert(is_arithmetic<RetT>::value, __FUNCTION__ " : Return type must arithmetic type.");
		static_assert(is_floating_point<FloatT>::value, __FUNCTION__ "FloatT must floating point type.");

		using PromotionType = promotion_t<RetT, FloatT, T>;

		PromotionType exp[4];
		Details::catmullRom_support(static_cast<PromotionType>(t), exp);
		return Details::catmullRom_calculate<RetT, PromotionType, T>(exp, v0, v1, v2, v3);
	}

	//input range v0, v1, t : 0 ~ 1.
	//ret : 0 ~ 1
	template<typename RetT, typename T, typename U, typename FloatT>
	RetT smoothStep(const T& v0, const T& v1, const FloatT& t) {
		static_assert(is_arithmetic<T>::value, __FUNCTION__ " : must arithmetic type.");
		static_assert(is_arithmetic<RetT>::value, __FUNCTION__ " : Return type must arithmetic type.");
		static_assert(is_floating_point<FloatT>::value, __FUNCTION__ "FloatT must floating point type.");

		using PromotionType = promotion_t<RetT, FloatT, T>;

		const PromotionType t2 = saturate<PromotionType, PromotionType>(div<PromotionType, FloatT, FloatT>(t - static_cast<FloatT>(v0), static_cast<FloatT>(v1) - static_cast<FloatT>(v0)));
		return static_cast<RetT>(
			t2 * t2 *(3 - (t2 * 2))
		);
	}
#pragma region(assign)
	template<typename RetType, typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol>
	auto assign(mat_base<LeftType, LeftRow, LeftCol>& lhs, const mat_base<RightType, RightRow, RightCol>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return lhs;
	}

	template<typename RetType, typename LeftType, typename RightType, int LeftDim, int RightDim>
	auto assign(vec_base<LeftType, LeftDim>& lhs, const vec_base<RightType, RightDim>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return lhs;
	}

	//n x m matrix
	template<typename T, typename U, int Row, int Col>
	mat_base<T, Row, Col>& assign(mat_base<T, Row, Col>& dst, const mat_base<U, Row, Col>& src) {
		for (int r = 0; r != Row; ++r) {
			assign<T, U, Col>(dst[r], src[r]);
		}
		return dst;
	}

	//n dim vector
	template<typename T, typename U, int Dim>
	vec_base<T, Dim>& assign(vec_base<T, Dim>& dst, const vec_base<U, Dim>& src) {
		for (int d = 0; d != Dim; ++d) {
			assign(dst[d], src[d]);
		}
		return dst;
	}

	//4 x 4 matrix
	template<typename T, typename U>
	mat4x4_t<T>& assign(mat4x4_t<T>& dst, const mat4x4_t<U>& src) {
		assign(dst[0], src[0]);
		assign(dst[1], src[1]);
		assign(dst[2], src[2]);
		assign(dst[3], src[3]);
		return dst;
	}

	//3 x 3 matrix
	template<typename T, typename U>
	mat3x3_t<T>& assign(mat3x3_t<T>& dst, const mat3x3_t<U>& src) {
		assign(dst[0], src[0]);
		assign(dst[1], src[1]);
		assign(dst[2], src[2]);
		return dst;
	}

	//2 x 2 matrix
	template<typename T, typename U>
	mat2x2_t<T>& assign(mat2x2_t<T>& dst, const mat2x2_t<U>& src) {
		assign(dst[0], src[0]);
		assign(dst[1], src[1]);
		return dst;
	}

	//4 dim vector
	template<typename T, typename U>
	vec4_t<T>& assign(vec4_t<T>& dst, const vec4_t<U>& src) {
		assign(dst[0], src[0]);
		assign(dst[1], src[1]);
		assign(dst[2], src[2]);
		assign(dst[3], src[3]);
		return dst;
	}

	//3 dim vector
	template<typename T, typename U>
	vec3_t<T>& assign(vec3_t<T>& dst, const vec3_t<U>& src) {
		assign(dst[0], src[0]);
		assign(dst[1], src[1]);
		assign(dst[2], src[2]);
		return dst;
	}

	//2 dim vector
	template<typename T, typename U>
	vec2_t<T>& assign(vec2_t<T>& dst, const vec2_t<U>& src) {
		assign(dst[0], src[0]);
		assign(dst[1], src[1]);
		return dst;
	}

	//quat
	template<typename T, typename U>
	quat_base<T>& assign(quat_base<T>& dst, const quat_base<U>& src) {
		assign(dst[0], src[0]);
		assign(dst[1], src[1]);
		assign(dst[2], src[2]);
		assign(dst[3], src[3]);
		return dst;
	}

	//complex
	template<typename T, typename U>
	complex_base<T>& assign(complex_base<T>& dst, const complex_base<U>& src) {
		assign(dst[0], src[0]);
		assign(dst[1], src[1]);
		return dst;
	}
#pragma endregion

//binary
#pragma region(is_equal)

	template<typename LeftType, typename RightType, int LeftRow, int LeftCol , int RightRow, int RightCol>
	bool is_equal(const mat_base<LeftType, LeftRow, LeftCol>& lhs, const mat_base<RightType, RightRow, RightCol>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int LeftDim, int RightDim>
	bool is_equal(const vec_base<LeftType, LeftDim>& lhs, const vec_base<RightType, RightDim>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_equal(const complex_base<LeftType>& c, const vec_base<RightType, Dim>& v) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_equal(const vec_base<LeftType, Dim>& v, const complex_base<RightType>& c) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_equal(const quat_base<LeftType>& q, const vec_base<RightType, Dim>& v) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_equal(const vec_base<LeftType, Dim>& v, const quat_base<RightType>& q) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	//n x m matrix
	template<typename LeftType, typename RightType, int Row, int Col>
	bool is_equal(const mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		for (int i = 0; i != Row; ++i) {
			if (is_not_equal(lhs[i], rhs[i])) 
				return false;
		}
		return true;
	}

	//n dim vector 
	template<typename LeftType, typename RightType, int Dim>
	bool is_equal(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		for (int i = 0; i != Dim; ++i) {
			if (is_not_equal(lhs[i], rhs[i]))
				return false;
		}
		return true;
	}

	//4 x 4 matrix
	template<typename LeftType, typename RightType>
	bool is_equal(const mat4x4_t<LeftType>& lhs, const mat4x4_t<RightType>& rhs) {
		return is_equal(lhs[0], rhs[0]) && is_equal(lhs[1], rhs[1]) && is_equal(lhs[2], rhs[2]) && is_equal(lhs[3], rhs[3]);
	}

	//3 x 3 matrix
	template<typename LeftType, typename RightType>
	bool is_equal(const mat3x3_t<LeftType>& lhs, const mat3x3_t<RightType>& rhs) {
		return is_equal(lhs[0], rhs[0]) && is_equal(lhs[1], rhs[1]) && is_equal(lhs[2], rhs[2]);
	}

	//2 x 2 matrix
	template<typename LeftType, typename RightType>
	bool is_equal(const mat2x2_t<LeftType>& lhs, const mat2x2_t<RightType>& rhs) {
		return is_equal(lhs[0], rhs[0]) && is_equal(lhs[1], rhs[1]);
	}

	//4 dim vector 
	template<typename LeftType, typename RightType>
	bool is_equal(const vec4_t<LeftType>& lhs, const vec4_t<RightType>& rhs) {
		return is_equal(lhs[0], rhs[0]) && is_equal(lhs[1], rhs[1]) && is_equal(lhs[2], rhs[2]) && is_equal(lhs[3], rhs[3]);
	}

	//3 dim vector 
	template<typename LeftType, typename RightType>
	bool is_equal(const vec3_t<LeftType>& lhs, const vec3_t<RightType>& rhs) {
		return is_equal(lhs[0], rhs[0]) && is_equal(lhs[1], rhs[1]) && is_equal(lhs[2], rhs[2]);
	}

	//2 dim vector 
	template<typename LeftType, typename RightType>
	bool is_equal(const vec2_t<LeftType>& lhs, const vec2_t<RightType>& rhs) {
		return is_equal(lhs[0], rhs[0]) && is_equal(lhs[1], rhs[1]);
	}

	//quat
	template<typename LeftType, typename RightType>
	bool is_equal(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		return is_equal(lhs[0], rhs[0]) && is_equal(lhs[1], rhs[1]) && is_equal(lhs[2], rhs[2]) && is_equal(lhs[3], rhs[3]);
	}

	//complex
	template<typename LeftType, typename RightType>
	bool is_equal(const complex_base<LeftType>& lhs, const complex_base<RightType>& rhs) {
		return is_equal(lhs[0], rhs[0]) && is_equal(lhs[1], rhs[1]);
	}

	//vec4t quat 
	template<typename LeftType, typename RightType>
	bool is_equal(const vec4_t<LeftType>& v, const quat_base<RightType>& q) {
		return is_equal(v[0], q[0]) && is_equal(v[1], q[1]) && is_equal(v[2], q[2]) && is_equal(v[3], q[3]);
	}

	//quat vec4t
	template<typename LeftType, typename RightType>
	bool is_equal(const quat_base<LeftType>& q, const vec4_t<RightType>& v) {
		return is_equal(v, q);
	}

	//vec2t complex
	template<typename LeftType, typename RightType>
	bool is_equal(const vec2_t<LeftType>& v, const complex_base<RightType>& c) {
		return is_equal(v[0], c[0]) && is_equal(v[1], c[1]);
	}

	//complex vec2t
	template<typename LeftType, typename RightType>
	bool is_equal(const complex_base<LeftType>& c , const vec2_t<RightType> & v) {
		return is_equal(v, c);
	}

#pragma endregion

//binary
#pragma region(is_not_equal)

	template<typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol>
	bool is_not_equal(const mat_base<LeftType, LeftRow, LeftCol>& lhs, const mat_base<RightType, RightRow, RightCol>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int LeftDim, int RightDim>
	bool is_not_equal(const vec_base<LeftType, LeftDim>& lhs, const vec_base<RightType, RightDim>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_not_equal(const complex_base<LeftType>& c, const vec_base<RightType, Dim>& v) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_not_equal(const vec_base<LeftType, Dim>& v, const complex_base<RightType>& c) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_not_equal(const quat_base<LeftType>& q, const vec_base<RightType, Dim>& v) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_not_equal(const vec_base<LeftType, Dim>& v, const quat_base<RightType>& q) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	//n x m matrix
	template<typename LeftType, typename RightType, int Row, int Col>
	bool is_not_equal(const mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		for (int i = 0; i != Row; ++i) {
			if (is_equal(lhs[i], rhs[i]))
				return false;
		}
		return true;
	}

	//n dim vector 
	template<typename LeftType, typename RightType, int Dim>
	bool is_not_equal(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		for (int i = 0; i != Dim; ++i) {
			if (is_equal(lhs[i], rhs[i]))
				return false;
		}
		return true;
	}

	//4 x 4 matrix
	template<typename LeftType, typename RightType>
	bool is_not_equal(const mat4x4_t<LeftType>& lhs, const mat4x4_t<RightType>& rhs) {
		return is_not_equal(lhs[0], rhs[0]) && is_not_equal(lhs[1], rhs[1]) && is_not_equal(lhs[2], rhs[2]) && is_not_equal(lhs[3], rhs[3]);
	}

	//3 x 3 matrix
	template<typename LeftType, typename RightType>
	bool is_not_equal(const mat3x3_t<LeftType>& lhs, const mat3x3_t<RightType>& rhs) {
		return is_not_equal(lhs[0], rhs[0]) && is_not_equal(lhs[1], rhs[1]) && is_not_equal(lhs[2], rhs[2]);
	}

	//2 x 2 matrix
	template<typename LeftType, typename RightType>
	bool is_not_equal(const mat2x2_t<LeftType>& lhs, const mat2x2_t<RightType>& rhs) {
		return is_not_equal(lhs[0], rhs[0]) && is_not_equal(lhs[1], rhs[1]);
	}

	//4 dim vector 
	template<typename LeftType, typename RightType>
	bool is_not_equal(const vec4_t<LeftType>& lhs, const vec4_t<RightType>& rhs) {
		return is_not_equal(lhs[0], rhs[0]) && is_not_equal(lhs[1], rhs[1]) && is_not_equal(lhs[2], rhs[2]) && is_not_equal(lhs[3], rhs[3]);
	}

	//3 dim vector 
	template<typename LeftType, typename RightType>
	bool is_not_equal(const vec3_t<LeftType>& lhs, const vec3_t<RightType>& rhs) {
		return is_not_equal(lhs[0], rhs[0]) && is_not_equal(lhs[1], rhs[1]) && is_not_equal(lhs[2], rhs[2]);
	}

	//2 dim vector 
	template<typename LeftType, typename RightType>
	bool is_not_equal(const vec2_t<LeftType>& lhs, const vec2_t<RightType>& rhs) {
		return is_not_equal(lhs[0], rhs[0]) && is_not_equal(lhs[1], rhs[1]);
	}

	//quat
	template<typename LeftType, typename RightType>
	bool is_not_equal(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		return is_not_equal(lhs[0], rhs[0]) && is_not_equal(lhs[1], rhs[1]) && is_not_equal(lhs[2], rhs[2]) && is_not_equal(lhs[3], rhs[3]);
	}

	//complex
	template<typename LeftType, typename RightType>
	bool is_not_equal(const complex_base<LeftType>& lhs, const complex_base<RightType>& rhs) {
		return is_not_equal(lhs[0], rhs[0]) && is_not_equal(lhs[1], rhs[1]);
	}

	//vec4t quat 
	template<typename LeftType, typename RightType>
	bool is_not_equal(const vec4_t<LeftType>& v, const quat_base< RightType>& q) {
		return is_not_equal(v[0], q[0]) && is_not_equal(v[1], q[1]) && is_not_equal(v[2], q[2]) && is_not_equal(v[3], q[3]);
	}

	//quat vec4t
	template<typename LeftType, typename RightType>
	bool is_not_equal(const quat_base<LeftType>& q, const vec4_t<RightType>& v) {
		return is_not_equal(v, q);
	}

	//vec2t complex
	template<typename LeftType, typename RightType>
	bool is_not_equal(const vec2_t<LeftType>& v, const complex_base<RightType>& c) {
		return is_not_equal(v[0], c[0]) && is_not_equal(v[1], c[1]);
	}

	//complex vec2t
	template<typename LeftType, typename RightType>
	bool is_not_equal(const complex_base<RightType>& c, const vec2_t<RightType>& v) {
		return is_not_equal(v, c);
	}

#pragma endregion

//binary
#pragma region(is_near_equal)
	template<typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol, typename EpsT>
	bool is_near_equal(const mat_base<LeftType, LeftRow, LeftCol>& lhs, const mat_base<RightType, RightRow, RightCol>& rhs, const EpsT& e) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int LeftDim, int RightDim, typename EpsT>
	bool is_near_equal(const vec_base<LeftType, LeftDim>& lhs, const vec_base<RightType, RightDim>& rhs, const EpsT& e) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim, typename EpsT>
	bool is_near_equal(const vec_base<LeftType, Dim>& v, const complex_base<RightType>& c, const EpsT& e) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim, typename EpsT>
	bool is_near_equal(const complex_base<LeftType>& c, const vec_base<RightType, Dim>& v, const EpsT& e) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim, typename EpsT>
	bool is_near_equal(const vec_base<LeftType, Dim>& v, const quat_base<RightType>& q, const EpsT& e) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim, typename EpsT>
	bool is_near_equal(const quat_base<LeftType>& q, const vec_base<RightType, Dim>& v, const EpsT& e) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	//n x m matrix
	template<typename LeftType, typename RightType, int Row, int Col, typename EpsT>
	bool is_near_equal(const mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs, const EpsT& e) {
		for (int i = 0; i != Row; ++i) {
			if (is_near_not_equal(lhs[i], rhs[i], e)) return false;
		}
		return true;
	}

	//n dim vector 
	template<typename LeftType, typename RightType, int Dim, typename EpsT>
	bool is_near_equal(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs, const EpsT& e) {
		for (int i = 0; i != Dim; ++i) {
			if (is_near_not_equal(lhs[i], rhs[i]))
				return false;
		}
		return true;
	}

	//4 x 4 matrix
	template<typename LeftType, typename RightType, typename EpsT>
	bool is_near_equal(const mat4x4_t<LeftType>& lhs, const mat4x4_t<RightType>& rhs, const EpsT& e) {
		return is_near_equal(lhs[0], rhs[0], e) && is_near_equal(lhs[1], rhs[1], e) && is_near_equal(lhs[2], rhs[2], e) && is_near_equal(lhs[3], rhs[3], e);
	}

	//3 x 3 matrix
	template<typename LeftType, typename RightType, typename EpsT>
	bool is_near_equal(const mat3x3_t<LeftType>& lhs, const mat3x3_t<RightType>& rhs, const EpsT& e) {
		return is_near_equal(lhs[0], rhs[0], e) && is_near_equal(lhs[1], rhs[1], e) && is_near_equal(lhs[2], rhs[2], e);
	}

	//2 x 2 matrix
	template<typename LeftType, typename RightType, typename EpsT>
	bool is_near_equal(const mat2x2_t<LeftType>& lhs, const mat2x2_t<RightType>& rhs, const EpsT& e) {
		return is_near_equal(lhs[0], rhs[0], e) && is_near_equal(lhs[1], rhs[1], e);
	}

	//4 dim vector 
	template<typename LeftType, typename RightType, typename EpsT>
	bool is_near_equal(const vec4_t<LeftType>& lhs, const vec4_t<RightType>& rhs, const EpsT& e) {
		return is_near_equal(lhs[0], rhs[0], e) && is_near_equal(lhs[1], rhs[1], e) && is_near_equal(lhs[2], rhs[2], e) && is_near_equal(lhs[3], rhs[3], e);
	}

	//3 dim vector 
	template<typename LeftType, typename RightType, typename EpsT>
	bool is_near_equal(const vec3_t<LeftType>& lhs, const vec3_t<RightType>& rhs, const EpsT& e) {
		return is_near_equal(lhs[0], rhs[0], e) && is_near_equal(lhs[1], rhs[1], e) && is_near_equal(lhs[2], rhs[2], e);
	}

	//2 dim vector 
	template<typename LeftType, typename RightType, typename EpsT>
	bool is_near_equal(const vec2_t<LeftType>& lhs, const vec2_t<RightType>& rhs, const EpsT& e) {
		return is_near_equal(lhs[0], rhs[0], e) && is_near_equal(lhs[1], rhs[1], e);
	}

	//quat
	template<typename LeftType, typename RightType, typename EpsT>
	bool is_near_equal(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs, const EpsT& e) {
		return is_near_equal(lhs[0], rhs[0], e) && is_near_equal(lhs[1], rhs[1], e) && is_near_equal(lhs[2], rhs[2], e) && is_near_equal(lhs[3], rhs[3], e);
	}

	//complex
	template<typename LeftType, typename RightType, typename EpsT>
	bool is_near_equal(const complex_base<LeftType>& lhs, const complex_base<RightType>& rhs, const EpsT& e) {
		return is_near_equal(lhs[0], rhs[0], e) && is_near_equal(lhs[1], rhs[1], e);
	}

	//vec4t quat 
	template<typename LeftType, typename RightType, typename EpsT>
	bool is_near_equal(const vec4_t<LeftType> & v, const quat_base< RightType>& q, const EpsT& e) {
		return is_near_equal(v[0], q[0], e) && is_near_equal(v[1], q[1], e) && is_near_equal(v[2], q[2], e) && is_near_equal(v[3], q[3], e);
	}

	//quat vec4t
	template<typename LeftType, typename RightType, typename EpsT>
	bool is_near_equal(const quat_base<LeftType>& q, const vec4_t<RightType>& v, const EpsT& e) {
		return is_near_equal(v, q, e);
	}

	//vec2t complex
	template<typename LeftType, typename RightType, typename EpsT>
	bool is_near_equal(const vec2_t<LeftType>& v, const complex_base<RightType>& c, const EpsT& e) {
		return is_near_equal(v[0], c[0], e) && is_near_equal(v[1], c[1], e);
	}

	//complex vec2t
	template<typename LeftType, typename RightType, typename EpsT>
	bool is_near_equal(const complex_base<LeftType>& c, const vec2_t<RightType>& v, const EpsT& e) {
		return is_near_equal(v, c, e);
	}

#pragma endregion

//binary
#pragma region(is_not_near_equal)
	template<typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol, typename EpsT>
	bool is_not_near_equal(const mat_base<LeftType, LeftRow, LeftCol>& lhs, const mat_base<RightType, RightRow, RightCol>& rhs, const EpsT& e) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int LeftDim, int RightDim, typename EpsT>
	bool is_not_near_equal(const vec_base<LeftType, LeftDim>& lhs, const vec_base<RightType, RightDim>& rhs, const EpsT& e) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim, typename EpsT>
	bool is_not_near_equal(const vec_base<LeftType, Dim>& v, const complex_base<RightType>& c, const EpsT& e) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim, typename EpsT>
	bool is_not_near_equal(const complex_base<LeftType>& c, const vec_base<RightType, Dim>& v, const EpsT& e) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim, typename EpsT>
	bool is_not_near_equal(const vec_base<LeftType, Dim>& v, const quat_base<RightType>& q, const EpsT& e) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim, typename EpsT>
	bool is_not_near_equal(const quat_base<LeftType>& q, const vec_base<RightType, Dim>& v, const EpsT& e) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	//n x m matrix
	template<typename LeftType, typename RightType, int Row, int Col, typename EpsT>
	bool is_not_near_equal(const mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs, const EpsT& e) {
		for (int i = 0; i != Row; ++i) {
			if (is_near_equal(lhs[i], rhs[i], e)) return false;
		}
		return true;
	}

	//n dim vector 
	template<typename LeftType, typename RightType, int Dim, typename EpsT>
	bool is_not_near_equal(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs, const EpsT& e) {
		for (int i = 0; i != Dim; ++i) {
			if (is_near_equal(lhs[i], rhs[i], e))
				return false;
		}
		return true;
	}

	//4 x 4 matrix
	template<typename LeftType, typename RightType, typename EpsT>
	bool is_not_near_equal(const mat4x4_t<LeftType>& lhs, const mat4x4_t<RightType>& rhs, const EpsT& e) {
		return is_not_near_equal(lhs[0], rhs[0], e) && is_not_near_equal(lhs[1], rhs[1], e) && is_not_near_equal(lhs[2], rhs[2], e) && is_not_near_equal(lhs[3], rhs[3], e);
	}

	//3 x 3 matrix
	template<typename LeftType, typename RightType, typename EpsT>
	bool is_not_near_equal(const mat3x3_t<LeftType>& lhs, const mat3x3_t<RightType>& rhs, const EpsT& e) {
		return is_not_near_equal(lhs[0], rhs[0], e) && is_not_near_equal(lhs[1], rhs[1], e) && is_not_near_equal(lhs[2], rhs[2], e);
	}

	//2 x 2 matrix
	template<typename LeftType, typename RightType, typename EpsT>
	bool is_not_near_equal(const mat2x2_t<LeftType>& lhs, const mat2x2_t<RightType>& rhs, const EpsT& e) {
		return is_not_near_equal(lhs[0], rhs[0], e) && is_not_near_equal(lhs[1], rhs[1], e);
	}

	//4 dim vector 
	template<typename LeftType, typename RightType, typename EpsT>
	bool is_not_near_equal(const vec4_t<LeftType>& lhs, const vec4_t<RightType>& rhs, const EpsT& e) {
		return is_not_near_equal(lhs[0], rhs[0], e) && is_not_near_equal(lhs[1], rhs[1], e) && is_not_near_equal(lhs[2], rhs[2], e) && is_not_near_equal(lhs[3], rhs[3], e);
	}

	//3 dim vector 
	template<typename LeftType, typename RightType, typename EpsT>
	bool is_not_near_equal(const vec3_t<LeftType>& lhs, const vec3_t<RightType>& rhs, const EpsT& e) {
		return is_not_near_equal(lhs[0], rhs[0], e) && is_not_near_equal(lhs[1], rhs[1], e) && is_not_near_equal(lhs[2], rhs[2], e);
	}

	//2 dim vector 
	template<typename LeftType, typename RightType, typename EpsT>
	bool is_not_near_equal(const vec2_t<LeftType>& lhs, const vec2_t<RightType>& rhs, const EpsT& e) {
		return is_not_near_equal(lhs[0], rhs[0], e) && is_not_near_equal(lhs[1], rhs[1], e);
	}

	//quat
	template<typename LeftType, typename RightType, typename EpsT>
	bool is_not_near_equal(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs, const EpsT& e) {
		return is_not_near_equal(lhs[0], rhs[0], e) && is_not_near_equal(lhs[1], rhs[1], e) && is_not_near_equal(lhs[2], rhs[2], e) && is_not_near_equal(lhs[3], rhs[3], e);
	}

	//complex
	template<typename LeftType, typename RightType, typename EpsT>
	bool is_not_near_equal(const complex_base<LeftType>& lhs, const complex_base<RightType>& rhs, const EpsT& e) {
		return is_not_near_equal(lhs[0], rhs[0], e) && is_not_near_equal(lhs[1], rhs[1], e);
	}

	//vec4t quat 
	template<typename LeftType, typename RightType, typename EpsT>
	bool is_not_near_equal(const vec4_t<LeftType>& v, const quat_base< RightType>& q, const EpsT& e) {
		return is_not_near_equal(v[0], q[0], e) && is_not_near_equal(v[1], q[1], e) && is_not_near_equal(v[2], q[2], e) && is_not_near_equal(v[3], q[3], e);
	}

	//quat vec4t
	template<typename LeftType, typename RightType, typename EpsT>
	bool is_not_near_equal(const quat_base<LeftType>& q, const vec4_t<RightType>& v, const EpsT& e) {
		return is_not_near_equal(v, q, e);
	}

	//vec2t complex
	template<typename LeftType, typename RightType, typename EpsT>
	bool is_not_near_equal(const vec2_t<LeftType>& v, const complex_base<RightType>& c, const EpsT& e) {
		return is_not_near_equal(v[0], c[0], e) && is_not_near_equal(v[1], c[1], e);
	}

	//complex vec2t
	template<typename LeftType, typename RightType, typename EpsT>
	bool is_not_near_equal(const complex_base<LeftType>& c, const vec2_t<RightType>& v, const EpsT& e) {
		return is_not_near_equal(v, c, e);
	}

#pragma endregion

//unary
#pragma region(is_zero)
	//n x m matrix
	template<typename Type, int Row, int Col>
	bool is_zero(const mat_base<Type, Row, Col>& v) {
		for (int i = 0; i != Row; ++i) {
			if (!is_zero(v[i])) return false;
		}
		return true;
	}

	//n dim vector 
	template<typename Type, int Dim>
	bool is_zero(const vec_base<Type, Dim>& v) {
		for (int i = 0; i != Dim; ++i) {
			if (!is_zero(v[i]))	return false;
		}
		return true;
	}

	//4 x 4 matrix
	template<typename Type>
	bool is_zero(const mat4x4_t<Type>& v) {
		return is_zero(v[0]) && is_zero(v[1]) && is_zero(v[2]) && is_zero(v[3]);
	}

	//3 x 3 matrix
	template<typename Type>
	bool is_zero(const mat3x3_t<Type>& v) {
		return is_zero(v[0]) && is_zero(v[1]) && is_zero(v[2]);
	}

	//2 x 2 matrix
	template<typename Type>
	bool is_zero(const mat2x2_t<Type>& v) {
		return is_zero(v[0]) && is_zero(v[1]);
	}

	//4 dim vector 
	template<typename Type>
	bool is_zero(const vec4_t<Type>& v) {
		return is_zero(v[0]) && is_zero(v[1]) && is_zero(v[2]) && is_zero(v[3]);
	}

	//3 dim vector 
	template<typename Type>
	bool is_zero(const vec3_t<Type>& v) {
		return is_zero(v[0]) && is_zero(v[1]) && is_zero(v[2]);
	}

	//2 dim vector 
	template<typename Type>
	bool is_zero(const vec2_t<Type>& v) {
		return is_zero(v[0]) && is_zero(v[1]);
	}

	//quat
	template<typename Type>
	bool is_zero(const quat_base<Type>& v) {
		return is_zero(v[0]) && is_zero(v[1]) && is_zero(v[2]) && is_zero(v[3]);
	}

	//complex
	template<typename Type>
	bool is_zero(const complex_base<Type>& v) {
		return is_zero(v[0]) && is_zero(v[1]);
	}

#pragma endregion

//unary
#pragma region(is_not_zero)
	//n x m matrix
	template<typename Type, int Row, int Col>
	bool is_not_zero(const mat_base<Type, Row, Col>& v) {
		for (int i = 0; i != Row; ++i) {
			if (is_zero(v[i])) 
				return false;
		}
		return true;
	}

	//n dim vector 
	template<typename Type, int Dim>
	bool is_not_zero(const vec_base<Type, Dim>& v) {
		for (int i = 0; i != Dim; ++i) {
			if (is_zero(v[i]))	
				return false;
		}
		return true;
	}

	//4 x 4 matrix
	template<typename Type>
	bool is_not_zero(const mat4x4_t<Type>& v) {
		return is_not_zero(v[0]) && is_not_zero(v[1]) && is_not_zero(v[2]) && is_not_zero(v[3]);
	}

	//3 x 3 matrix
	template<typename Type>
	bool is_not_zero(const mat3x3_t<Type>& v) {
		return is_not_zero(v[0]) && is_not_zero(v[1]) && is_not_zero(v[2]);
	}

	//2 x 2 matrix
	template<typename Type>
	bool is_not_zero(const mat2x2_t<Type>& v) {
		return is_not_zero(v[0]) && is_not_zero(v[1]);
	}

	//4 dim vector 
	template<typename Type>
	bool is_not_zero(const vec4_t<Type>& v) {
		return is_not_zero(v[0]) && is_not_zero(v[1]) && is_not_zero(v[2]) && is_not_zero(v[3]);
	}

	//3 dim vector 
	template<typename Type>
	bool is_not_zero(const vec3_t<Type>& v) {
		return is_not_zero(v[0]) && is_not_zero(v[1]) && is_not_zero(v[2]);
	}

	//2 dim vector 
	template<typename Type>
	bool is_not_zero(const vec2_t<Type>& v) {
		return is_not_zero(v[0]) && is_not_zero(v[1]);
	}

	//quat
	template<typename Type>
	bool is_not_zero(const quat_base<Type>& v) {
		return is_not_zero(v[0]) && is_not_zero(v[1]) && is_not_zero(v[2]) && is_not_zero(v[3]);
	}

	//complex
	template<typename Type>
	bool is_not_zero(const complex_base<Type>& v) {
		return is_not_zero(v[0]) && is_not_zero(v[1]);
	}

#pragma endregion

//unary
#pragma region(is_near_zero)
	//n x m matrix
	template<typename Type, int Row, int Col, typename EpsT>
	bool is_near_zero(const mat_base<Type, Row, Col>& v, const EpsT& e) {
		for (int i = 0; i != Row; ++i) {
			if (is_not_near_zero(v[i], e)) return false;
		}
		return true;
	}

	//n dim vector 
	template<typename Type, int Dim, typename EpsT>
	bool is_near_zero(const vec_base<Type, Dim>& v, const EpsT& e) {
		for (int i = 0; i != Dim; ++i) {
			if (is_not_near_zero(v[i], e)) return false;
		}
		return true;
	}

	//4 x 4 matrix
	template<typename Type, typename EpsT>
	bool is_near_zero(const mat4x4_t<Type>& v, const EpsT& e) {
		return is_near_zero(v[0], e) && is_near_zero(v[1], e) && is_near_zero(v[2], e) && is_near_zero(v[3], e);
	}

	//3 x 3 matrix
	template<typename Type, typename EpsT>
	bool is_near_zero(const mat3x3_t<Type>& v, const EpsT& e) {
		return is_near_zero(v[0], e) && is_near_zero(v[1], e) && is_near_zero(v[2], e);
	}

	//2 x 2 matrix
	template<typename Type, typename EpsT>
	bool is_near_zero(const mat2x2_t<Type>& v, const EpsT& e) {
		return is_near_zero(v[0], e) && is_near_zero(v[1], e);
	}

	//4 dim vector 
	template<typename Type, typename EpsT>
	bool is_near_zero(const vec4_t<Type>& v, const EpsT& e) {
		return is_near_zero(v[0], e) && is_near_zero(v[1], e) && is_near_zero(v[2], e) && is_near_zero(v[3], e);
	}

	//3 dim vector 
	template<typename Type, typename EpsT>
	bool is_near_zero(const vec3_t<Type>& v, const EpsT& e) {
		return is_near_zero(v[0], e) && is_near_zero(v[1], e) && is_near_zero(v[2], e);
	}

	//2 dim vector 
	template<typename Type, typename EpsT>
	bool is_near_zero(const vec2_t<Type>& v, const EpsT& e) {
		return is_near_zero(v[0], e) && is_near_zero(v[1], e);
	}

	//quat
	template<typename Type, typename EpsT>
	bool is_near_zero(const quat_base<Type>& v, const EpsT& e) {
		return is_near_zero(v[0], e) && is_near_zero(v[1], e) && is_near_zero(v[2], e) && is_near_zero(v[3], e);
	}

	//complex
	template<typename Type, typename EpsT>
	bool is_near_zero(const complex_base<Type>& v, const EpsT& e) {
		return is_near_zero(v[0], e) && is_near_zero(v[1], e);
	}

#pragma endregion

//unary
#pragma region(is_not_near_zero)
	//n x m matrix
	template<typename Type, int Row, int Col, typename EpsT>
	bool is_not_near_zero(const mat_base<Type, Row, Col>& v, const EpsT& e) {
		for (int i = 0; i != Row; ++i) {
			if (is_near_zero(v[i], e)) return false;
		}
		return true;
	}

	//n dim vector 
	template<typename Type, int Dim, typename EpsT>
	bool is_not_near_zero(const vec_base<Type, Dim>& v, const EpsT& e) {
		for (int i = 0; i != Dim; ++i) {
			if (is_near_zero(v[i], e)) return false;
		}
		return true;
	}

	//4 x 4 matrix
	template<typename Type, typename EpsT>
	bool is_not_near_zero(const mat4x4_t<Type>& v, const EpsT& e) {
		return is_not_near_zero(v[0], e) && is_not_near_zero(v[1], e) && is_not_near_zero(v[2], e) && is_not_near_zero(v[3], e);
	}

	//3 x 3 matrix
	template<typename Type, typename EpsT>
	bool is_not_near_zero(const mat3x3_t<Type>& v, const EpsT& e) {
		return is_not_near_zero(v[0], e) && is_not_near_zero(v[1], e) && is_not_near_zero(v[2], e);
	}

	//2 x 2 matrix
	template<typename Type, typename EpsT>
	bool is_not_near_zero(const mat2x2_t<Type>& v, const EpsT& e) {
		return is_not_near_zero(v[0], e) && is_not_near_zero(v[1], e);
	}

	//4 dim vector 
	template<typename Type, typename EpsT>
	bool is_not_near_zero(const vec4_t<Type>& v, const EpsT& e) {
		return is_not_near_zero(v[0], e) && is_not_near_zero(v[1], e) && is_not_near_zero(v[2], e) && is_not_near_zero(v[3], e);
	}

	//3 dim vector 
	template<typename Type, typename EpsT>
	bool is_not_near_zero(const vec3_t<Type>& v, const EpsT& e) {
		return is_not_near_zero(v[0], e) && is_not_near_zero(v[1], e) && is_not_near_zero(v[2], e);
	}

	//2 dim vector 
	template<typename Type, typename EpsT>
	bool is_not_near_zero(const vec2_t<Type>& v, const EpsT& e) {
		return is_not_near_zero(v[0], e) && is_not_near_zero(v[1], e);
	}

	//quat
	template<typename Type, typename EpsT>
	bool is_not_near_zero(const quat_base<Type>& v, const EpsT& e) {
		return is_not_near_zero(v[0], e) && is_not_near_zero(v[1], e) && is_not_near_zero(v[2], e) && is_not_near_zero(v[3], e);
	}

	//complex
	template<typename Type, typename EpsT>
	bool is_not_near_zero(const complex_base<Type>& v, const EpsT& e) {
		return is_not_near_zero(v[0], e) && is_not_near_zero(v[1], e);
	}

#pragma endregion

//binary
#pragma region(is_less)

	template<typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol>
	bool is_less(const mat_base<LeftType, LeftRow, LeftCol>& lhs, const mat_base<RightType, RightRow, RightCol>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int LeftDim, int RightDim>
	bool is_less(const vec_base<LeftType, LeftDim>& lhs, const vec_base<RightType, RightDim>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_less(const complex_base<LeftType>& c, const vec_base<RightType, Dim>& v) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_less(const vec_base<LeftType, Dim>& v, const complex_base<RightType>& c) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_less(const quat_base<LeftType>& q, const vec_base<RightType, Dim>& v) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_less(const vec_base<LeftType, Dim>& v, const quat_base<RightType>& q) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	//n x m matrix
	template<typename LeftType, typename RightType, int Row, int Col>
	bool is_less(const mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		for (int i = 0; i != Row; ++i) {
			if (!is_less(lhs[i], rhs[i]))
				return false;
		}
		return true;
	}

	//n dim vector 
	template<typename LeftType, typename RightType, int Dim>
	bool is_less(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		for (int i = 0; i != Dim; ++i) {
			if (!is_less(lhs[i], rhs[i]))
				return false;
		}
		return true;
	}

	//4 x 4 matrix
	template<typename LeftType, typename RightType>
	bool is_less(const mat4x4_t<LeftType>& lhs, const mat4x4_t<RightType>& rhs) {
		return is_less(lhs[0], rhs[0]) && is_less(lhs[1], rhs[1]) && is_less(lhs[2], rhs[2]) && is_less(lhs[3], rhs[3]);
	}

	//3 x 3 matrix
	template<typename LeftType, typename RightType>
	bool is_less(const mat3x3_t<LeftType>& lhs, const mat3x3_t<RightType>& rhs) {
		return is_less(lhs[0], rhs[0]) && is_less(lhs[1], rhs[1]) && is_less(lhs[2], rhs[2]);
	}

	//2 x 2 matrix
	template<typename LeftType, typename RightType>
	bool is_less(const mat2x2_t<LeftType>& lhs, const mat2x2_t<RightType>& rhs) {
		return is_less(lhs[0], rhs[0]) && is_less(lhs[1], rhs[1]);
	}

	//4 dim vector 
	template<typename LeftType, typename RightType>
	bool is_less(const vec4_t<LeftType>& lhs, const vec4_t<RightType>& rhs) {
		return is_less(lhs[0], rhs[0]) && is_less(lhs[1], rhs[1]) && is_less(lhs[2], rhs[2]) && is_less(lhs[3], rhs[3]);
	}

	//3 dim vector 
	template<typename LeftType, typename RightType>
	bool is_less(const vec3_t<LeftType>& lhs, const vec3_t<RightType>& rhs) {
		return is_less(lhs[0], rhs[0]) && is_less(lhs[1], rhs[1]) && is_less(lhs[2], rhs[2]);
	}

	//2 dim vector 
	template<typename LeftType, typename RightType>
	bool is_less(const vec2_t<LeftType>& lhs, const vec2_t<RightType>& rhs) {
		return is_less(lhs[0], rhs[0]) && is_less(lhs[1], rhs[1]);
	}

	//quat
	template<typename LeftType, typename RightType>
	bool is_less(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		return is_less(lhs[0], rhs[0]) && is_less(lhs[1], rhs[1]) && is_less(lhs[2], rhs[2]) && is_less(lhs[3], rhs[3]);
	}

	//complex
	template<typename LeftType, typename RightType>
	bool is_less(const complex_base<LeftType>& lhs, const complex_base<RightType>& rhs) {
		return is_less(lhs[0], rhs[0]) && is_less(lhs[1], rhs[1]);
	}

	//vec4t quat 
	template<typename LeftType, typename RightType>
	bool is_less(const vec4_t<LeftType> & v, const quat_base< RightType>& q) {
		return is_less(v[0], q[0]) && is_less(v[1], q[1]) && is_less(v[2], q[2]) && is_less(v[3], q[3]);
	}

	//quat vec4t
	template<typename LeftType, typename RightType>
	bool is_less(const quat_base<LeftType>& q, const vec4_t<RightType> & v) {
		return is_less(v, q);
	}

	//vec2t complex
	template<typename LeftType, typename RightType>
	bool is_less(const vec2_t<LeftType>& v, const complex_base<RightType>& c) {
		return is_less(v[0], c[0]) && is_less(v[1], c[1]);
	}

	//complex vec2t
	template<typename LeftType, typename RightType>
	bool is_less(const complex_base<LeftType>& c, const vec2_t<RightType>& v) {
		return is_less(v, c);
	}

#pragma endregion

//binary
#pragma region(is_less_equal)

	template<typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol>
	bool is_less_equal(const mat_base<LeftType, LeftRow, LeftCol>& lhs, const mat_base<RightType, RightRow, RightCol>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int LeftDim, int RightDim>
	bool is_less_equal(const vec_base<LeftType, LeftDim>& lhs, const vec_base<RightType, RightDim>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_less_equal(const complex_base<LeftType>& c, const vec_base<RightType, Dim>& v) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_less_equal(const vec_base<LeftType, Dim>& v, const complex_base<RightType>& c) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_less_equal(const quat_base<LeftType>& q, const vec_base<RightType, Dim>& v) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_less_equal(const vec_base<LeftType, Dim>& v, const quat_base<RightType>& q) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	//n x m matrix
	template<typename LeftType, typename RightType, int Row, int Col>
	bool is_less_equal(const mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		for (int i = 0; i != Row; ++i) {
			if (!is_less_equal(lhs[i], rhs[i]))
				return false;
		}
		return true;
	}

	//n dim vector 
	template<typename LeftType, typename RightType, int Dim>
	bool is_less_equal(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		for (int i = 0; i != Dim; ++i) {
			if (!is_less_equal(lhs[i], rhs[i]))
				return false;
		}
		return true;
	}

	//4 x 4 matrix
	template<typename LeftType, typename RightType>
	bool is_less_equal(const mat4x4_t<LeftType>& lhs, const mat4x4_t<RightType>& rhs) {
		return is_less_equal(lhs[0], rhs[0]) && is_less_equal(lhs[1], rhs[1]) && is_less_equal(lhs[2], rhs[2]) && is_less_equal(lhs[3], rhs[3]);
	}

	//3 x 3 matrix
	template<typename LeftType, typename RightType>
	bool is_less_equal(const mat3x3_t<LeftType>& lhs, const mat3x3_t<RightType>& rhs) {
		return is_less_equal(lhs[0], rhs[0]) && is_less_equal(lhs[1], rhs[1]) && is_less_equal(lhs[2], rhs[2]);
	}

	//2 x 2 matrix
	template<typename LeftType, typename RightType>
	bool is_less_equal(const mat2x2_t<LeftType>& lhs, const mat2x2_t<RightType>& rhs) {
		return is_less_equal(lhs[0], rhs[0]) && is_less_equal(lhs[1], rhs[1]);
	}

	//4 dim vector 
	template<typename LeftType, typename RightType>
	bool is_less_equal(const vec4_t<LeftType>& lhs, const vec4_t<RightType>& rhs) {
		return is_less_equal(lhs[0], rhs[0]) && is_less_equal(lhs[1], rhs[1]) && is_less_equal(lhs[2], rhs[2]) && is_less_equal(lhs[3], rhs[3]);
	}

	//3 dim vector 
	template<typename LeftType, typename RightType>
	bool is_less_equal(const vec3_t<LeftType>& lhs, const vec3_t<RightType>& rhs) {
		return is_less_equal(lhs[0], rhs[0]) && is_less_equal(lhs[1], rhs[1]) && is_less_equal(lhs[2], rhs[2]);
	}

	//2 dim vector 
	template<typename LeftType, typename RightType>
	bool is_less_equal(const vec2_t<LeftType>& lhs, const vec2_t<RightType>& rhs) {
		return is_less_equal(lhs[0], rhs[0]) && is_less_equal(lhs[1], rhs[1]);
	}

	//quat
	template<typename LeftType, typename RightType>
	bool is_less_equal(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		return is_less_equal(lhs[0], rhs[0]) && is_less_equal(lhs[1], rhs[1]) && is_less_equal(lhs[2], rhs[2]) && is_less_equal(lhs[3], rhs[3]);
	}

	//complex
	template<typename LeftType, typename RightType>
	bool is_less_equal(const complex_base<LeftType>& lhs, const complex_base<RightType>& rhs) {
		return is_less_equal(lhs[0], rhs[0]) && is_less_equal(lhs[1], rhs[1]);
	}

	//vec4t quat 
	template<typename LeftType, typename RightType>
	bool is_less_equal(const vec4_t<LeftType>& v, const quat_base< RightType>& q) {
		return is_less_equal(v[0], q[0]) && is_less_equal(v[1], q[1]) && is_less_equal(v[2], q[2]) && is_less_equal(v[3], q[3]);
	}

	//quat vec4t
	template<typename LeftType, typename RightType>
	bool is_less_equal(const quat_base<LeftType>& q, const vec4_t<RightType>& v) {
		return is_less_equal(v, q);
	}

	//vec2t complex
	template<typename LeftType, typename RightType>
	bool is_less_equal(const vec2_t<LeftType>& v, const complex_base<RightType>& c) {
		return is_less_equal(v[0], c[0]) && is_less_equal(v[1], c[1]);
	}

	//complex vec2t
	template<typename LeftType, typename RightType>
	bool is_less_equal(const complex_base<LeftType>& c, const vec2_t<RightType>& v) {
		return is_less_equal(v, c);
	}

#pragma endregion

//binary
#pragma region(is_greater)

	template<typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol>
	bool is_greater(const mat_base<LeftType, LeftRow, LeftCol>& lhs, const mat_base<RightType, RightRow, RightCol>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int LeftDim, int RightDim>
	bool is_greater(const vec_base<LeftType, LeftDim>& lhs, const vec_base<RightType, RightDim>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_greater(const complex_base<LeftType>& c, const vec_base<RightType, Dim>& v) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_greater(const vec_base<LeftType, Dim>& v, const complex_base<RightType>& c) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_greater(const quat_base<LeftType>& q, const vec_base<RightType, Dim>& v) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_greater(const vec_base<LeftType, Dim>& v, const quat_base<RightType>& q) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	//n x m matrix
	template<typename LeftType, typename RightType, int Row, int Col>
	bool is_greater(const mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		for (int i = 0; i != Row; ++i) {
			if (!is_greater(lhs[i], rhs[i]))
				return false;
		}
		return true;
	}

	//n dim vector 
	template<typename LeftType, typename RightType, int Dim>
	bool is_greater(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		for (int i = 0; i != Dim; ++i) {
			if (!is_greater(lhs[i], rhs[i]))
				return false;
		}
		return true;
	}

	//4 x 4 matrix
	template<typename LeftType, typename RightType>
	bool is_greater(const mat4x4_t<LeftType>& lhs, const mat4x4_t<RightType>& rhs) {
		return is_greater(lhs[0], rhs[0]) && is_greater(lhs[1], rhs[1]) && is_greater(lhs[2], rhs[2]) && is_greater(lhs[3], rhs[3]);
	}

	//3 x 3 matrix
	template<typename LeftType, typename RightType>
	bool is_greater(const mat3x3_t<LeftType>& lhs, const mat3x3_t<RightType>& rhs) {
		return is_greater(lhs[0], rhs[0]) && is_greater(lhs[1], rhs[1]) && is_greater(lhs[2], rhs[2]);
	}

	//2 x 2 matrix
	template<typename LeftType, typename RightType>
	bool is_greater(const mat2x2_t<LeftType>& lhs, const mat2x2_t<RightType>& rhs) {
		return is_greater(lhs[0], rhs[0]) && is_greater(lhs[1], rhs[1]);
	}

	//4 dim vector 
	template<typename LeftType, typename RightType>
	bool is_greater(const vec4_t<LeftType>& lhs, const vec4_t<RightType>& rhs) {
		return is_greater(lhs[0], rhs[0]) && is_greater(lhs[1], rhs[1]) && is_greater(lhs[2], rhs[2]) && is_greater(lhs[3], rhs[3]);
	}

	//3 dim vector 
	template<typename LeftType, typename RightType>
	bool is_greater(const vec3_t<LeftType>& lhs, const vec3_t<RightType>& rhs) {
		return is_greater(lhs[0], rhs[0]) && is_greater(lhs[1], rhs[1]) && is_greater(lhs[2], rhs[2]);
	}

	//2 dim vector 
	template<typename LeftType, typename RightType>
	bool is_greater(const vec2_t<LeftType>& lhs, const vec2_t<RightType>& rhs) {
		return is_greater(lhs[0], rhs[0]) && is_greater(lhs[1], rhs[1]);
	}

	//quat
	template<typename LeftType, typename RightType>
	bool is_greater(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		return is_greater(lhs[0], rhs[0]) && is_greater(lhs[1], rhs[1]) && is_greater(lhs[2], rhs[2]) && is_greater(lhs[3], rhs[3]);
	}

	//complex
	template<typename LeftType, typename RightType>
	bool is_greater(const complex_base<LeftType>& lhs, const complex_base<RightType>& rhs) {
		return is_greater(lhs[0], rhs[0]) && is_greater(lhs[1], rhs[1]);
	}

	//vec4t quat 
	template<typename LeftType, typename RightType>
	bool is_greater(const vec4_t<LeftType>& v, const quat_base<RightType>& q) {
		return is_greater(v[0], q[0]) && is_greater(v[1], q[1]) && is_greater(v[2], q[2]) && is_greater(v[3], q[3]);
	}

	//quat vec4t
	template<typename LeftType, typename RightType>
	bool is_greater(const quat_base<LeftType>& q, const vec4_t<RightType>& v) {
		return is_greater(v, q);
	}

	//vec2t complex
	template<typename LeftType, typename RightType>
	bool is_greater(const vec2_t<LeftType> & v, const complex_base<RightType>& c) {
		return is_greater(v[0], c[0]) && is_greater(v[1], c[1]);
	}

	//complex vec2t
	template<typename LeftType, typename RightType>
	bool is_greater(const complex_base<LeftType>& c, const vec2_t<RightType>& v) {
		return is_greater(v, c);
	}

#pragma endregion

//binary
#pragma region(is_greater_equal)

	template<typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol>
	bool is_greater_equal(const mat_base<LeftType, LeftRow, LeftCol>& lhs, const mat_base<RightType, RightRow, RightCol>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int LeftDim, int RightDim>
	bool is_greater_equal(const vec_base<LeftType, LeftDim>& lhs, const vec_base<RightType, RightDim>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_greater_equal(const complex_base<LeftType>& c, const vec_base<RightType, Dim>& v) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_greater_equal(const vec_base<LeftType, Dim>& v, const complex_base<RightType>& c) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_greater_equal(const quat_base<LeftType>& q, const vec_base<RightType, Dim>& v) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_greater_equal(const vec_base<LeftType, Dim>& v, const quat_base<RightType>& q) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	//n x m matrix
	template<typename LeftType, typename RightType, int Row, int Col>
	bool is_greater_equal(const mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		for (int i = 0; i != Row; ++i) {
			if (!is_greater_equal(lhs[i], rhs[i]))
				return false;
		}
		return true;
	}

	//n dim vector 
	template<typename LeftType, typename RightType, int Dim>
	bool is_greater_equal(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		for (int i = 0; i != Dim; ++i) {
			if (!is_greater_equal(lhs[i], rhs[i]))
				return false;
		}
		return true;
	}

	//4 x 4 matrix
	template<typename LeftType, typename RightType>
	bool is_greater_equal(const mat4x4_t<LeftType>& lhs, const mat4x4_t<RightType>& rhs) {
		return is_greater_equal(lhs[0], rhs[0]) && is_greater_equal(lhs[1], rhs[1]) && is_greater_equal(lhs[2], rhs[2]) && is_greater_equal(lhs[3], rhs[3]);
	}

	//3 x 3 matrix
	template<typename LeftType, typename RightType>
	bool is_greater_equal(const mat3x3_t<LeftType>& lhs, const mat3x3_t<RightType>& rhs) {
		return is_greater_equal(lhs[0], rhs[0]) && is_greater_equal(lhs[1], rhs[1]) && is_greater_equal(lhs[2], rhs[2]);
	}

	//2 x 2 matrix
	template<typename LeftType, typename RightType>
	bool is_greater_equal(const mat2x2_t<LeftType>& lhs, const mat2x2_t<RightType>& rhs) {
		return is_greater_equal(lhs[0], rhs[0]) && is_greater_equal(lhs[1], rhs[1]);
	}

	//4 dim vector 
	template<typename LeftType, typename RightType>
	bool is_greater_equal(const vec4_t<LeftType>& lhs, const vec4_t<RightType>& rhs) {
		return is_greater_equal(lhs[0], rhs[0]) && is_greater_equal(lhs[1], rhs[1]) && is_greater_equal(lhs[2], rhs[2]) && is_greater_equal(lhs[3], rhs[3]);
	}

	//3 dim vector 
	template<typename LeftType, typename RightType>
	bool is_greater_equal(const vec3_t<LeftType>& lhs, const vec3_t<RightType>& rhs) {
		return is_greater_equal(lhs[0], rhs[0]) && is_greater_equal(lhs[1], rhs[1]) && is_greater_equal(lhs[2], rhs[2]);
	}

	//2 dim vector 
	template<typename LeftType, typename RightType>
	bool is_greater_equal(const vec2_t<LeftType>& lhs, const vec2_t<RightType>& rhs) {
		return is_greater_equal(lhs[0], rhs[0]) && is_greater_equal(lhs[1], rhs[1]);
	}

	//quat
	template<typename LeftType, typename RightType>
	bool is_greater_equal(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		return is_greater_equal(lhs[0], rhs[0]) && is_greater_equal(lhs[1], rhs[1]) && is_greater_equal(lhs[2], rhs[2]) && is_greater_equal(lhs[3], rhs[3]);
	}

	//complex
	template<typename LeftType, typename RightType>
	bool is_greater_equal(const complex_base<LeftType>& lhs, const complex_base<RightType>& rhs) {
		return is_greater_equal(lhs[0], rhs[0]) && is_greater_equal(lhs[1], rhs[1]);
	}

	//vec4t quat 
	template<typename LeftType, typename RightType>
	bool is_greater_equal(const vec4_t<LeftType>& v, const quat_base<RightType>& q) {
		return is_greater_equal(v[0], q[0]) && is_greater_equal(v[1], q[1]) && is_greater_equal(v[2], q[2]) && is_greater_equal(v[3], q[3]);
	}

	//quat vec4t
	template<typename LeftType, typename RightType>
	bool is_greater_equal(const quat_base<LeftType>& q, const vec4_t<RightType>& v) {
		return is_greater_equal(v, q);
	}

	//vec2t complex
	template<typename LeftType, typename RightType>
	bool is_greater_equal(const vec2_t<LeftType>& v, const complex_base<RightType>& c) {
		return is_greater_equal(v[0], c[0]) && is_greater_equal(v[1], c[1]);
	}

	//complex vec2t
	template<typename LeftType, typename RightType>
	bool is_greater_equal(const complex_base<LeftType>& c, const vec2_t<RightType>& v) {
		return is_greater_equal(v, c);
	}

#pragma endregion

//binary
#pragma region(is_equal_each)
	template<typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol>
	bool is_equal_each(const mat_base<LeftType, LeftRow, LeftCol>& lhs, const mat_base<RightType, RightRow, RightCol>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int LeftDim, int RightDim>
	bool is_equal_each(const vec_base<LeftType, LeftDim>& lhs, const vec_base<RightType, RightDim>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_equal_each(const complex_base<LeftType>& c, const vec_base<RightType, Dim>& v) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_equal_each(const vec_base<LeftType, Dim>& v, const complex_base<RightType>& c) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_equal_each(const quat_base<LeftType>& q, const vec_base<RightType, Dim>& v) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Row, int Col>
	bools<Row, Col> is_equal_each(const mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		bools<Row, Col> ret;
		for (int row = 0; row != Row; ++row) {
			ret[row] = is_equal_each<LeftType, RightType, Col>(lhs[row], rhs[row]);
		}
		return ret;
	}

	template<typename LeftType, typename RightType, int Dim>
	bools<Dim, 0> is_equal_each(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		bools<Dim, 0> ret;
		for (int d = 0; d != Dim; ++d) {
			ret[d] = is_equal<LeftType, RightType>(lhs[d], rhs[d]);
		}
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<4, 0> is_equal_each(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		bools<4, 0> ret;
		ret[0] = is_equal<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_equal<LeftType, RightType>(lhs[1], rhs[1]);
		ret[2] = is_equal<LeftType, RightType>(lhs[2], rhs[2]);
		ret[3] = is_equal<LeftType, RightType>(lhs[3], rhs[3]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<2, 0> is_equal_each(const complex_base<LeftType>& lhs, const complex_base<RightType>& rhs) {
		bools<2, 0> ret;
		ret[0] = is_equal<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_equal<LeftType, RightType>(lhs[1], rhs[1]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<4, 0> is_equal_each(const quat_base<LeftType>& lhs, const vec4_t<RightType>& rhs) {
		bools<4, 0> ret;
		ret[0] = is_equal<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_equal<LeftType, RightType>(lhs[1], rhs[1]);
		ret[2] = is_equal<LeftType, RightType>(lhs[2], rhs[2]);
		ret[3] = is_equal<LeftType, RightType>(lhs[3], rhs[3]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<4, 0> is_equal_each(const vec4_t<LeftType>& lhs, const quat_base<RightType>& rhs) {
		bools<4, 0> ret;
		ret[0] = is_equal<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_equal<LeftType, RightType>(lhs[1], rhs[1]);
		ret[2] = is_equal<LeftType, RightType>(lhs[2], rhs[2]);
		ret[3] = is_equal<LeftType, RightType>(lhs[3], rhs[3]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<2, 0> is_equal_each(const complex_base<LeftType>& lhs, const vec2_t<RightType>& rhs) {
		bools<2, 0> ret;
		ret[0] = is_equal<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_equal<LeftType, RightType>(lhs[1], rhs[1]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<2, 0> is_equal_each(const vec2_t<LeftType>& lhs, const complex_base<RightType>& rhs) {
		bools<2, 0> ret;
		ret[0] = is_equal<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_equal<LeftType, RightType>(lhs[1], rhs[1]);
		return ret;
	}
#pragma endregion

//binary
#pragma region(is_not_equal_each)
	template<typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol>
	bool is_not_equal_each(const mat_base<LeftType, LeftRow, LeftCol>& lhs, const mat_base<RightType, RightRow, RightCol>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int LeftDim, int RightDim>
	bool is_not_equal_each(const vec_base<LeftType, LeftDim>& lhs, const vec_base<RightType, RightDim>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_not_equal_each(const complex_base<LeftType>& c, const vec_base<RightType, Dim>& v) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_not_equal_each(const vec_base<LeftType, Dim>& v, const complex_base<RightType>& c) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_not_equal_each(const quat_base<LeftType>& q, const vec_base<RightType, Dim>& v) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Row, int Col>
	bools<Row, Col> is_not_equal_each(const mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		bools<Row, Col> ret;
		for (int row = 0; row != Row; ++row) {
			ret[row] = is_not_equal_each<LeftType, RightType, Col>(lhs[row], rhs[row]);
		}
		return ret;
	}

	template<typename LeftType, typename RightType, int Dim>
	bools<Dim, 0> is_not_equal_each(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		bools<Dim, 0> ret;
		for (int d = 0; d != Dim; ++d) {
			ret[d] = is_not_equal<LeftType, RightType>(lhs[d], rhs[d]);
		}
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<4, 0> is_not_equal_each(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		bools<4, 0> ret;
		ret[0] = is_not_equal<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_not_equal<LeftType, RightType>(lhs[1], rhs[1]);
		ret[2] = is_not_equal<LeftType, RightType>(lhs[2], rhs[2]);
		ret[3] = is_not_equal<LeftType, RightType>(lhs[3], rhs[3]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<2, 0> is_not_equal_each(const complex_base<LeftType>& lhs, const complex_base<RightType>& rhs) {
		bools<2, 0> ret;
		ret[0] = is_not_equal<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_not_equal<LeftType, RightType>(lhs[1], rhs[1]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<4, 0> is_not_equal_each(const quat_base<LeftType>& lhs, const vec4_t<RightType>& rhs) {
		bools<4, 0> ret;
		ret[0] = is_not_equal<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_not_equal<LeftType, RightType>(lhs[1], rhs[1]);
		ret[2] = is_not_equal<LeftType, RightType>(lhs[2], rhs[2]);
		ret[3] = is_not_equal<LeftType, RightType>(lhs[3], rhs[3]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<4, 0> is_not_equal_each(const vec4_t<LeftType>& lhs, const quat_base<RightType>& rhs) {
		bools<4, 0> ret;
		ret[0] = is_not_equal<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_not_equal<LeftType, RightType>(lhs[1], rhs[1]);
		ret[2] = is_not_equal<LeftType, RightType>(lhs[2], rhs[2]);
		ret[3] = is_not_equal<LeftType, RightType>(lhs[3], rhs[3]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<2, 0> is_not_equal_each(const complex_base<LeftType>& lhs, const vec2_t<RightType>& rhs) {
		bools<2, 0> ret;
		ret[0] = is_not_equal<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_not_equal<LeftType, RightType>(lhs[1], rhs[1]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<2, 0> is_not_equal_each(const vec2_t<LeftType>& lhs, const complex_base<RightType>& rhs) {
		bools<2, 0> ret;
		ret[0] = is_not_equal<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_not_equal<LeftType, RightType>(lhs[1], rhs[1]);
		return ret;
	}

#pragma endregion

//binary
#pragma region(is_near_equal_each)
	template<typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol, typename EpsT>
	bool is_near_equal_each(const mat_base<LeftType, LeftRow, LeftCol>& lhs, const mat_base<RightType, RightRow, RightCol>& rhs, const EpsT& e) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int LeftDim, int RightDim, typename EpsT>
	bool is_near_equal_each(const vec_base<LeftType, LeftDim>& lhs, const vec_base<RightType, RightDim>& rhs, const EpsT& e) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim, typename EpsT>
	bool is_near_equal_each(const complex_base<LeftType>& c, const vec_base<RightType, Dim>& v, const EpsT& e) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim, typename EpsT>
	bool is_near_equal_each(const vec_base<LeftType, Dim>& v, const complex_base<RightType>& c, const EpsT& e) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim, typename EpsT>
	bool is_near_equal_each(const quat_base<LeftType>& q, const vec_base<RightType, Dim>& v, const EpsT& e) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Row, int Col, typename EpsT>
	bools<Row, Col> is_near_equal_each(const mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs, const EpsT& e) {
		bools<Row, Col> ret;
		for (int row = 0; row != Row; ++row) {
			ret[row] = is_near_equal_each<LeftType, RightType, Col>(lhs[row], rhs[row], e);
		}
		return ret;
	}

	template<typename LeftType, typename RightType, int Dim, typename EpsT>
	bools<Dim, 0> is_near_equal_each(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs, const EpsT& e) {
		bools<Dim, 0> ret;
		for (int d = 0; d != Dim; ++d) {
			ret[d] = is_near_equal<LeftType, RightType>(lhs[d], rhs[d], e);
		}
		return ret;
	}

	template<typename LeftType, typename RightType, typename EpsT>
	bools<4, 0> is_near_equal_each(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs, const EpsT& e) {
		bools<4, 0> ret;
		ret[0] = is_near_equal<LeftType, RightType>(lhs[0], rhs[0], e);
		ret[1] = is_near_equal<LeftType, RightType>(lhs[1], rhs[1], e);
		ret[2] = is_near_equal<LeftType, RightType>(lhs[2], rhs[2], e);
		ret[3] = is_near_equal<LeftType, RightType>(lhs[3], rhs[3], e);
		return ret;
	}

	template<typename LeftType, typename RightType, typename EpsT>
	bools<2, 0> is_near_equal_each(const complex_base<LeftType>& lhs, const complex_base<RightType>& rhs, const EpsT& e) {
		bools<2, 0> ret;
		ret[0] = is_near_equal<LeftType, RightType>(lhs[0], rhs[0], e);
		ret[1] = is_near_equal<LeftType, RightType>(lhs[1], rhs[1], e);
		return ret;
	}

	template<typename LeftType, typename RightType, typename EpsT>
	bools<4, 0> is_near_equal_each(const quat_base<LeftType>& lhs, const vec4_t<RightType>& rhs, const EpsT& e) {
		bools<4, 0> ret;
		ret[0] = is_near_equal<LeftType, RightType>(lhs[0], rhs[0], e);
		ret[1] = is_near_equal<LeftType, RightType>(lhs[1], rhs[1], e);
		ret[2] = is_near_equal<LeftType, RightType>(lhs[2], rhs[2], e);
		ret[3] = is_near_equal<LeftType, RightType>(lhs[3], rhs[3], e);
		return ret;
	}

	template<typename LeftType, typename RightType, typename EpsT>
	bools<4, 0> is_near_equal_each(const vec4_t<LeftType>& lhs, const quat_base<RightType>& rhs, const EpsT& e) {
		bools<4, 0> ret;
		ret[0] = is_near_equal<LeftType, RightType>(lhs[0], rhs[0], e);
		ret[1] = is_near_equal<LeftType, RightType>(lhs[1], rhs[1], e);
		ret[2] = is_near_equal<LeftType, RightType>(lhs[2], rhs[2], e);
		ret[3] = is_near_equal<LeftType, RightType>(lhs[3], rhs[3], e);
		return ret;
	}

	template<typename LeftType, typename RightType, typename EpsT>
	bools<2, 0> is_near_equal_each(const complex_base<LeftType>& lhs, const vec2_t<RightType>& rhs, const EpsT& e) {
		bools<2, 0> ret;
		ret[0] = is_near_equal<LeftType, RightType>(lhs[0], rhs[0], e);
		ret[1] = is_near_equal<LeftType, RightType>(lhs[1], rhs[1], e);
		return ret;
	}

	template<typename LeftType, typename RightType, typename EpsT>
	bools<2, 0> is_near_equal_each(const vec2_t<LeftType>& lhs, const complex_base<RightType>& rhs, const EpsT& e) {
		bools<2, 0> ret;
		ret[0] = is_near_equal<LeftType, RightType>(lhs[0], rhs[0], e);
		ret[1] = is_near_equal<LeftType, RightType>(lhs[1], rhs[1], e);
		return ret;
	}
#pragma endregion

//binary
#pragma region(is_not_near_equal_each)
	template<typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol, typename EpsT>
	bool is_not_near_equal_each(const mat_base<LeftType, LeftRow, LeftCol>& lhs, const mat_base<RightType, RightRow, RightCol>& rhs, const EpsT& e) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int LeftDim, int RightDim, typename EpsT>
	bool is_not_near_equal_each(const vec_base<LeftType, LeftDim>& lhs, const vec_base<RightType, RightDim>& rhs, const EpsT& e) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim, typename EpsT>
	bool is_not_near_equal_each(const complex_base<LeftType>& c, const vec_base<RightType, Dim>& v, const EpsT& e) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim, typename EpsT>
	bool is_not_near_equal_each(const vec_base<LeftType, Dim>& v, const complex_base<RightType>& c, const EpsT& e) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim, typename EpsT>
	bool is_not_near_equal_each(const quat_base<LeftType>& q, const vec_base<RightType, Dim>& v, const EpsT& e) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Row, int Col, typename EpsT>
	bools<Row, Col> is_not_near_equal_each(const mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs, const EpsT& e) {
		bools<Row, Col> ret;
		for (int row = 0; row != Row; ++row) {
			ret[row] = is_not_near_equal_each<LeftType, RightType, Col>(lhs[row], rhs[row], e);
		}
		return ret;
	}

	template<typename LeftType, typename RightType, int Dim, typename EpsT>
	bools<Dim, 0> is_not_near_equal_each(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs, const EpsT& e) {
		bools<Dim, 0> ret;
		for (int d = 0; d != Dim; ++d) {
			ret[d] = is_not_near_equal<LeftType, RightType>(lhs[d], rhs[d], e);
		}
		return ret;
	}

	template<typename LeftType, typename RightType, typename EpsT>
	bools<4, 0> is_not_near_equal_each(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs, const EpsT& e) {
		bools<4, 0> ret;
		ret[0] = is_not_near_equal<LeftType, RightType>(lhs[0], rhs[0], e);
		ret[1] = is_not_near_equal<LeftType, RightType>(lhs[1], rhs[1], e);
		ret[2] = is_not_near_equal<LeftType, RightType>(lhs[2], rhs[2], e);
		ret[3] = is_not_near_equal<LeftType, RightType>(lhs[3], rhs[3], e);
		return ret;
	}

	template<typename LeftType, typename RightType, typename EpsT>
	bools<2, 0> is_not_near_equal_each(const complex_base<LeftType>& lhs, const complex_base<RightType>& rhs, const EpsT& e) {
		bools<2, 0> ret;
		ret[0] = is_not_near_equal<LeftType, RightType>(lhs[0], rhs[0], e);
		ret[1] = is_not_near_equal<LeftType, RightType>(lhs[1], rhs[1], e);
		return ret;
	}

	template<typename LeftType, typename RightType, typename EpsT>
	bools<4, 0> is_not_near_equal_each(const quat_base<LeftType>& lhs, const vec4_t<RightType>& rhs, const EpsT& e) {
		bools<4, 0> ret;
		ret[0] = is_not_near_equal<LeftType, RightType>(lhs[0], rhs[0], e);
		ret[1] = is_not_near_equal<LeftType, RightType>(lhs[1], rhs[1], e);
		ret[2] = is_not_near_equal<LeftType, RightType>(lhs[2], rhs[2], e);
		ret[3] = is_not_near_equal<LeftType, RightType>(lhs[3], rhs[3], e);
		return ret;
	}

	template<typename LeftType, typename RightType, typename EpsT>
	bools<4, 0> is_not_near_equal_each(const vec4_t<LeftType>& lhs, const quat_base<RightType>& rhs, const EpsT& e) {
		bools<4, 0> ret;
		ret[0] = is_not_near_equal<LeftType, RightType>(lhs[0], rhs[0], e);
		ret[1] = is_not_near_equal<LeftType, RightType>(lhs[1], rhs[1], e);
		ret[2] = is_not_near_equal<LeftType, RightType>(lhs[2], rhs[2], e);
		ret[3] = is_not_near_equal<LeftType, RightType>(lhs[3], rhs[3], e);
		return ret;
	}

	template<typename LeftType, typename RightType, typename EpsT>
	bools<2, 0> is_not_near_equal_each(const complex_base<LeftType>& lhs, const vec2_t<RightType>& rhs, const EpsT& e) {
		bools<2, 0> ret;
		ret[0] = is_not_near_equal<LeftType, RightType>(lhs[0], rhs[0], e);
		ret[1] = is_not_near_equal<LeftType, RightType>(lhs[1], rhs[1], e);
		return ret;
	}

	template<typename LeftType, typename RightType, typename EpsT>
	bools<2, 0> is_not_near_equal_each(const vec2_t<LeftType>& lhs, const complex_base<RightType>& rhs, const EpsT& e) {
		bools<2, 0> ret;
		ret[0] = is_not_near_equal<LeftType, RightType>(lhs[0], rhs[0], e);
		ret[1] = is_not_near_equal<LeftType, RightType>(lhs[1], rhs[1], e);
		return ret;
	}
#pragma endregion

//unary
#pragma region(is_near_zero_each)
	template<typename LeftType, int Row, int Col, typename EpsT>
	bools<Row, Col> is_near_zero_each(const mat_base<LeftType, Row, Col>& v, const EpsT& e) {
		bools<Row, Col> ret;
		for (int row = 0; row != Row; ++row) {
			ret[row] = is_near_zero_each<LeftType, Col>(v[row], e);
		}
		return ret;
	}

	template<typename LeftType, int Dim, typename EpsT>
	bools<Dim, 0> is_near_zero_each(const vec_base<LeftType, Dim>& v, const EpsT& e) {
		bools<Dim, 0> ret;
		for (int d = 0; d != Dim; ++d) {
			ret[d] = is_near_zero<LeftType>(v[d], e);
		}
		return ret;
	}

	template<typename LeftType, typename EpsT>
	bools<4, 0> is_near_zero_each(const quat_base<LeftType>& v, const EpsT& e) {
		bools<4, 0> ret;
		ret[0] = is_near_zero<LeftType>(v[0], e);
		ret[1] = is_near_zero<LeftType>(v[1], e);
		ret[2] = is_near_zero<LeftType>(v[2], e);
		ret[3] = is_near_zero<LeftType>(v[3], e);
		return ret;
	}

	template<typename LeftType, typename EpsT>
	bools<2, 0> is_near_zero_each(const complex_base<LeftType>& v, const EpsT& e) {
		bools<2, 0> ret;
		ret[0] = is_near_zero<LeftType>(v[0], e);
		ret[1] = is_near_zero<LeftType>(v[1], e);
		return ret;
	}		
#pragma endregion

//unary
#pragma region(is_zero_each)
	template<typename Type, int Row, int Col>
	bools<Row, Col> is_zero_each(const mat_base<Type, Row, Col>& v) {
		bools<Row, Col> ret;
		for (int row = 0; row != Row; ++row) {
			ret[row] = is_zero_each<Type, Col>(v[row]);
		}
		return ret;
	}

	template<typename Type, int Dim>
	bools<Dim, 0> is_zero_each(const vec_base<Type, Dim>& v) {
		bools<Dim, 0> ret;
		for (int d = 0; d != Dim; ++d) {
			ret[d] = is_zero<Type>(v[d]);
		}
		return ret;
	}

	template<typename Type>
	bools<4, 0> is_zero_each(const quat_base<Type>& v) {
		bools<4, 0> ret;
		ret[0] = is_zero<Type>(v[0]);
		ret[1] = is_zero<Type>(v[1]);
		ret[2] = is_zero<Type>(v[2]);
		ret[3] = is_zero<Type>(v[3]);
		return ret;
	}

	template<typename Type>
	bools<2, 0> is_zero_each(const complex_base<Type>& v) {
		bools<2, 0> ret;
		ret[0] = is_zero<Type>(v[0]);
		ret[1] = is_zero<Type>(v[1]);
		return ret;
	}
#pragma endregion

//binary
#pragma region(is_less_each)
	template<typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol>
	bool is_less_each(const mat_base<LeftType, LeftRow, LeftCol>& lhs, const mat_base<RightType, RightRow, RightCol>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int LeftDim, int RightDim>
	bool is_less_each(const vec_base<LeftType, LeftDim>& lhs, const vec_base<RightType, RightDim>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_less_each(const complex_base<LeftType>& c, const vec_base<RightType, Dim>& v) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_less_each(const vec_base<LeftType, Dim>& v, const complex_base<RightType>& c) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_less_each(const quat_base<LeftType>& q, const vec_base<RightType, Dim>& v) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Row, int Col>
	bools<Row, Col> is_less_each(const mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		bools<Row, Col> ret;
		for (int row = 0; row != Row; ++row) {
			ret[row] = is_less_each<LeftType, RightType, Col>(lhs[row], rhs[row]);
		}
		return ret;
	}

	template<typename LeftType, typename RightType, int Dim>
	bools<Dim, 0> is_less_each(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		bools<Dim, 0> ret;
		for (int d = 0; d != Dim; ++d) {
			ret[d] = is_less<LeftType, RightType>(lhs[d], rhs[d]);
		}
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<4, 0> is_less_each(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		bools<4, 0> ret;
		ret[0] = is_less<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_less<LeftType, RightType>(lhs[1], rhs[1]);
		ret[2] = is_less<LeftType, RightType>(lhs[2], rhs[2]);
		ret[3] = is_less<LeftType, RightType>(lhs[3], rhs[3]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<2, 0> is_less_each(const complex_base<LeftType>& lhs, const complex_base<RightType>& rhs) {
		bools<2, 0> ret;
		ret[0] = is_less<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_less<LeftType, RightType>(lhs[1], rhs[1]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<4, 0> is_less_each(const quat_base<LeftType>& lhs, const vec4_t<RightType>& rhs) {
		bools<4, 0> ret;
		ret[0] = is_less<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_less<LeftType, RightType>(lhs[1], rhs[1]);
		ret[2] = is_less<LeftType, RightType>(lhs[2], rhs[2]);
		ret[3] = is_less<LeftType, RightType>(lhs[3], rhs[3]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<4, 0> is_less_each(const vec4_t<LeftType>& lhs, const quat_base<RightType>& rhs) {
		bools<4, 0> ret;
		ret[0] = is_less<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_less<LeftType, RightType>(lhs[1], rhs[1]);
		ret[2] = is_less<LeftType, RightType>(lhs[2], rhs[2]);
		ret[3] = is_less<LeftType, RightType>(lhs[3], rhs[3]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<2, 0> is_less_each(const complex_base<LeftType>& lhs, const vec2_t<RightType>& rhs) {
		bools<2, 0> ret;
		ret[0] = is_less<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_less<LeftType, RightType>(lhs[1], rhs[1]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<2, 0> is_less_each(const vec2_t<LeftType>& lhs, const complex_base<RightType>& rhs) {
		bools<2, 0> ret;
		ret[0] = is_less<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_less<LeftType, RightType>(lhs[1], rhs[1]);
		return ret;
	}
#pragma endregion

//binary
#pragma region(is_less_equal_each)
	template<typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol>
	bool is_less_equal_each(const mat_base<LeftType, LeftRow, LeftCol>& lhs, const mat_base<RightType, RightRow, RightCol>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int LeftDim, int RightDim>
	bool is_less_equal_each(const vec_base<LeftType, LeftDim>& lhs, const vec_base<RightType, RightDim>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_less_equal_each(const complex_base<LeftType>& c, const vec_base<RightType, Dim>& v) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_less_equal_each(const vec_base<LeftType, Dim>& v, const complex_base<RightType>& c) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_less_equal_each(const quat_base<LeftType>& q, const vec_base<RightType, Dim>& v) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Row, int Col>
	bools<Row, Col> is_less_equal_each(const mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		bools<Row, Col> ret;
		for (int row = 0; row != Row; ++row) {
			ret[row] = is_less_equal_each<LeftType, RightType, Col>(lhs[row], rhs[row]);
		}
		return ret;
	}

	template<typename LeftType, typename RightType, int Dim>
	bools<Dim, 0> is_less_equal_each(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		bools<Dim, 0> ret;
		for (int d = 0; d != Dim; ++d) {
			ret[d] = is_less_equal<LeftType, RightType>(lhs[d], rhs[d]);
		}
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<4, 0> is_less_equal_each(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		bools<4, 0> ret;
		ret[0] = is_less_equal<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_less_equal<LeftType, RightType>(lhs[1], rhs[1]);
		ret[2] = is_less_equal<LeftType, RightType>(lhs[2], rhs[2]);
		ret[3] = is_less_equal<LeftType, RightType>(lhs[3], rhs[3]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<2, 0> is_less_equal_each(const complex_base<LeftType>& lhs, const complex_base<RightType>& rhs) {
		bools<2, 0> ret;
		ret[0] = is_less_equal<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_less_equal<LeftType, RightType>(lhs[1], rhs[1]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<4, 0> is_less_equal_each(const quat_base<LeftType>& lhs, const vec4_t<RightType>& rhs) {
		bools<4, 0> ret;
		ret[0] = is_less_equal<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_less_equal<LeftType, RightType>(lhs[1], rhs[1]);
		ret[2] = is_less_equal<LeftType, RightType>(lhs[2], rhs[2]);
		ret[3] = is_less_equal<LeftType, RightType>(lhs[3], rhs[3]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<4, 0> is_less_equal_each(const vec4_t<LeftType>& lhs, const quat_base<RightType>& rhs) {
		bools<4, 0> ret;
		ret[0] = is_less_equal<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_less_equal<LeftType, RightType>(lhs[1], rhs[1]);
		ret[2] = is_less_equal<LeftType, RightType>(lhs[2], rhs[2]);
		ret[3] = is_less_equal<LeftType, RightType>(lhs[3], rhs[3]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<2, 0> is_less_equal_each(const complex_base<LeftType>& lhs, const vec2_t<RightType>& rhs) {
		bools<2, 0> ret;
		ret[0] = is_less_equal<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_less_equal<LeftType, RightType>(lhs[1], rhs[1]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<2, 0> is_less_equal_each(const vec2_t<LeftType>& lhs, const complex_base<RightType>& rhs) {
		bools<2, 0> ret;
		ret[0] = is_less_equal<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_less_equal<LeftType, RightType>(lhs[1], rhs[1]);
		return ret;
	}
#pragma endregion

//binary 
#pragma region(is_greater_each)
	template<typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol>
	bool is_greater_each(const mat_base<LeftType, LeftRow, LeftCol>& lhs, const mat_base<RightType, RightRow, RightCol>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int LeftDim, int RightDim>
	bool is_greater_each(const vec_base<LeftType, LeftDim>& lhs, const vec_base<RightType, RightDim>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_greater_each(const complex_base<LeftType>& c, const vec_base<RightType, Dim>& v) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_greater_each(const vec_base<LeftType, Dim>& v, const complex_base<RightType>& c) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_greater_each(const quat_base<LeftType>& q, const vec_base<RightType, Dim>& v) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Row, int Col>
	bools<Row, Col> is_greater_each(const mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		bools<Row, Col> ret;
		for (int row = 0; row != Row; ++row) {
			ret[row] = is_greater_each<LeftType, RightType, Col>(lhs[row], rhs[row]);
		}
		return ret;
	}

	template<typename LeftType, typename RightType, int Dim>
	bools<Dim, 0> is_greater_each(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		bools<Dim, 0> ret;
		for (int d = 0; d != Dim; ++d) {
			ret[d] = is_greater<LeftType, RightType>(lhs[d], rhs[d]);
		}
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<4, 0> is_greater_each(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		bools<4, 0> ret;
		ret[0] = is_greater<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_greater<LeftType, RightType>(lhs[1], rhs[1]);
		ret[2] = is_greater<LeftType, RightType>(lhs[2], rhs[2]);
		ret[3] = is_greater<LeftType, RightType>(lhs[3], rhs[3]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<2, 0> is_greater_each(const complex_base<LeftType>& lhs, const complex_base<RightType>& rhs) {
		bools<2, 0> ret;
		ret[0] = is_greater<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_greater<LeftType, RightType>(lhs[1], rhs[1]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<4, 0> is_greater_each(const quat_base<LeftType>& lhs, const vec4_t<RightType>& rhs) {
		bools<4, 0> ret;
		ret[0] = is_greater<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_greater<LeftType, RightType>(lhs[1], rhs[1]);
		ret[2] = is_greater<LeftType, RightType>(lhs[2], rhs[2]);
		ret[3] = is_greater<LeftType, RightType>(lhs[3], rhs[3]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<4, 0> is_greater_each(const vec4_t<LeftType>& lhs, const quat_base<RightType>& rhs) {
		bools<4, 0> ret;
		ret[0] = is_greater<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_greater<LeftType, RightType>(lhs[1], rhs[1]);
		ret[2] = is_greater<LeftType, RightType>(lhs[2], rhs[2]);
		ret[3] = is_greater<LeftType, RightType>(lhs[3], rhs[3]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<2, 0> is_greater_each(const complex_base<LeftType>& lhs, const vec2_t<RightType>& rhs) {
		bools<2, 0> ret;
		ret[0] = is_greater<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_greater<LeftType, RightType>(lhs[1], rhs[1]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<2, 0> is_greater_each(const vec2_t<LeftType>& lhs, const complex_base<RightType>& rhs) {
		bools<2, 0> ret;
		ret[0] = is_greater<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_greater<LeftType, RightType>(lhs[1], rhs[1]);
		return ret;
	}
#pragma endregion

//binary
#pragma region(is_greater_equal_each)
	template<typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol>
	bool is_greater_equal_each(const mat_base<LeftType, LeftRow, LeftCol>& lhs, const mat_base<RightType, RightRow, RightCol>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int LeftDim, int RightDim>
	bool is_greater_equal_each(const vec_base<LeftType, LeftDim>& lhs, const vec_base<RightType, RightDim>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_greater_equal_each(const complex_base<LeftType>& c, const vec_base<RightType, Dim>& v) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_greater_equal_each(const vec_base<LeftType, Dim>& v, const complex_base<RightType>& c) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Dim>
	bool is_greater_equal_each(const quat_base<LeftType>& q, const vec_base<RightType, Dim>& v) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	template<typename LeftType, typename RightType, int Row, int Col>
	bools<Row, Col> is_greater_equal_each(const mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		bools<Row, Col> ret;
		for (int row = 0; row != Row; ++row) {
			ret[row] = is_greater_equal_each<LeftType, RightType, Col>(lhs[row], rhs[row]);
		}
		return ret;
	}

	template<typename LeftType, typename RightType, int Dim>
	bools<Dim, 0> is_greater_equal_each(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		bools<Dim, 0> ret;
		for (int d = 0; d != Dim; ++d) {
			ret[d] = is_greater_equal<LeftType, RightType>(lhs[d], rhs[d]);
		}
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<4, 0> is_greater_equal_each(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		bools<4, 0> ret;
		ret[0] = is_greater_equal<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_greater_equal<LeftType, RightType>(lhs[1], rhs[1]);
		ret[2] = is_greater_equal<LeftType, RightType>(lhs[2], rhs[2]);
		ret[3] = is_greater_equal<LeftType, RightType>(lhs[3], rhs[3]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<2, 0> is_greater_equal_each(const complex_base<LeftType>& lhs, const complex_base<RightType>& rhs) {
		bools<2, 0> ret;
		ret[0] = is_greater_equal<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_greater_equal<LeftType, RightType>(lhs[1], rhs[1]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<4, 0> is_greater_equal_each(const quat_base<LeftType>& lhs, const vec4_t<RightType>& rhs) {
		bools<4, 0> ret;
		ret[0] = is_greater_equal<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_greater_equal<LeftType, RightType>(lhs[1], rhs[1]);
		ret[2] = is_greater_equal<LeftType, RightType>(lhs[2], rhs[2]);
		ret[3] = is_greater_equal<LeftType, RightType>(lhs[3], rhs[3]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<4, 0> is_greater_equal_each(const vec4_t<LeftType>& lhs, const quat_base<RightType>& rhs) {
		bools<4, 0> ret;
		ret[0] = is_greater_equal<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_greater_equal<LeftType, RightType>(lhs[1], rhs[1]);
		ret[2] = is_greater_equal<LeftType, RightType>(lhs[2], rhs[2]);
		ret[3] = is_greater_equal<LeftType, RightType>(lhs[3], rhs[3]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<2, 0> is_greater_equal_each(const complex_base<LeftType>& lhs, const vec2_t<RightType>& rhs) {
		bools<2, 0> ret;
		ret[0] = is_greater_equal<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_greater_equal<LeftType, RightType>(lhs[1], rhs[1]);
		return ret;
	}

	template<typename LeftType, typename RightType>
	bools<2, 0> is_greater_equal_each(const vec2_t<LeftType>& lhs, const complex_base<RightType>& rhs) {
		bools<2, 0> ret;
		ret[0] = is_greater_equal<LeftType, RightType>(lhs[0], rhs[0]);
		ret[1] = is_greater_equal<LeftType, RightType>(lhs[1], rhs[1]);
		return ret;
	}
#pragma endregion

//unary 
#pragma region(abs)
	template<typename RetType, typename Type, int Row, int Col>
	auto abs(const mat_base<Type, Row, Col>& v) {
		using RetT = mat_base<RetType, Row, Col>;
		RetT ret;
		for (int r = 0; r != Row; ++r) {
			ret[r] = abs<RetType, Type, Col>(v[r]);
		}
		return ret;
	}

	template<typename RetType, typename Type, int Dim>
	auto abs(const vec_base<Type, Dim>& v) {
		using RetT = vec_base<RetType, Dim>;
		RetT ret;
		for (int d = 0; d != Dim; ++d) {
			ret[d] = abs<RetType, Type>(v[d]);
		}
		return ret;
	}

	template<typename RetType, typename Type>
	auto abs(const mat4x4_t<Type>& v) {
		using RetT = mat4x4_t<RetType>;
		return RetT{
			abs<RetType, Type>(v[0]),
			abs<RetType, Type>(v[1]),
			abs<RetType, Type>(v[2]),
			abs<RetType, Type>(v[3])
		};
	}

	template<typename RetType, typename Type>
	auto abs(const mat3x3_t<Type>& v) {
		using RetT = mat3x3_t<RetType>;
		return RetT{
			abs<RetType, Type>(v[0]),
			abs<RetType, Type>(v[1]),
			abs<RetType, Type>(v[2])
		};
	}

	template<typename RetType, typename Type>
	auto abs(const mat2x2_t<Type>& v) {
		using RetT = mat2x2_t<RetType>;
		return RetT{
			abs<RetType, Type>(v[0]),
			abs<RetType, Type>(v[1])
		};
	}

	template<typename RetType, typename Type>
	auto abs(const vec4_t<Type>& v) {
		using RetT = vec4_t<RetType>;
		return RetT{
			abs<RetType, Type>(v[0]),
			abs<RetType, Type>(v[1]),
			abs<RetType, Type>(v[2]),
			abs<RetType, Type>(v[3])
		};
	}

	template<typename RetType, typename Type>
	auto abs(const vec3_t<Type>& v) {
		using RetT = vec3_t<RetType>;
		return RetT{
			abs<RetType, Type>(v[0]),
			abs<RetType, Type>(v[1]),
			abs<RetType, Type>(v[2])
		};
	}

	template<typename RetType, typename Type>
	auto abs(const vec2_t<Type>& v) {
		using RetT = vec2_t<RetType>;
		return RetT{
			abs<RetType, Type>(v[0]),
			abs<RetType, Type>(v[1])
		};
	}

	template<typename RetType, typename Type>
	auto abs(const quat_base<Type>& q) {
		using RetT = quat_base<RetType>;
		return RetT{
			abs<RetType, Type>(q[0]),
			abs<RetType, Type>(q[1]),
			abs<RetType, Type>(q[2]),
			abs<RetType, Type>(q[3])
		};
	}

	template<typename RetType, typename Type>
	auto abs(const complex_base<Type>& c) {
		using RetT = complex_base<RetType>;
		return RetT{
			abs<RetType, Type>(c[0]),
			abs<RetType, Type>(c[1])
		};
	}

#pragma endregion

//unary 
#pragma region(minus)
	template<typename RetType, typename Type, int Row, int Col>
	auto minus(const mat_base<Type, Row, Col>& v) {
		using RetT = mat_base<RetType, Row, Col>;
		RetT ret;
		for (int r = 0; r != Row; ++r) {
			ret[r] = minus<RetType, Type, Col>(v[r]);
		}
		return ret;
	}

	template<typename RetType, typename Type, int Dim>
	auto minus(const vec_base<Type, Dim>& v) {
		using RetT = vec_base<RetType, Dim>;
		RetT ret;
		for (int d = 0; d != Dim; ++d) {
			ret[d] = minus<RetType, Type>(v[d]);
		}
		return ret;
	}

	template<typename RetType, typename Type>
	auto minus(const mat4x4_t<Type>& v) {
		using RetT = mat4x4_t<RetType>;
		return RetT{
			minus<RetType, Type>(v[0]),
			minus<RetType, Type>(v[1]),
			minus<RetType, Type>(v[2]),
			minus<RetType, Type>(v[3])
		};
	}

	template<typename RetType, typename Type>
	auto minus(const mat3x3_t<Type>& v) {
		using RetT = mat3x3_t<RetType>;
		return RetT{
			minus<RetType, Type>(v[0]),
			minus<RetType, Type>(v[1]),
			minus<RetType, Type>(v[2])
		};
	}

	template<typename RetType, typename Type>
	auto minus(const mat2x2_t<Type>& v) {
		using RetT = mat2x2_t<RetType>;
		return RetT{
			minus<RetType, Type>(v[0]),
			minus<RetType, Type>(v[1])
		};
	}
	
	template<typename RetType, typename Type>
	auto minus(const vec4_t<Type>& v) {
		using RetT = vec4_t<RetType>;
		return RetT{
			minus<RetType, Type>(v[0]),
			minus<RetType, Type>(v[1]),
			minus<RetType, Type>(v[2]),
			minus<RetType, Type>(v[3])
		};
	}

	template<typename RetType, typename Type>
	auto minus(const vec3_t<Type>& v) {
		using RetT = vec3_t<RetType>;
		return RetT{
			minus<RetType, Type>(v[0]),
			minus<RetType, Type>(v[1]),
			minus<RetType, Type>(v[2])
		};
	}

	template<typename RetType, typename Type>
	auto minus(const vec2_t<Type>& v) {
		using RetT = vec2_t<RetType>;
		return RetT{
			minus<RetType, Type>(v[0]),
			minus<RetType, Type>(v[1])
		};
	}

	template<typename RetType, typename Type>
	auto minus(const quat_base<Type>& q) {
		using RetT = quat_base<RetType>;
		return RetT{
			minus<RetType, Type>(q[0]),
			minus<RetType, Type>(q[1]),
			minus<RetType, Type>(q[2]),
			minus<RetType, Type>(q[3])
		};
	}

	template<typename RetType, typename Type>
	auto minus(const complex_base<Type>& c) {
		using RetT = complex_base<RetType>;
		return RetT{
			minus<RetType, Type>(c[0]),
			minus<RetType, Type>(c[1])
		};
	}

#pragma endregion
	
//unary 
#pragma region(plus)
	template<typename RetType, typename Type, int Row, int Col>
	auto plus(const mat_base<Type, Row, Col>& v) {
		using RetT = mat_base<RetType, Row, Col>;
		RetT ret;
		for (int r = 0; r != Row; ++r) {
			ret[r] = plus<RetType, Type, Col>(v[r]);
		}
		return ret;
	}

	template<typename RetType, typename Type, int Dim>
	auto plus(const vec_base<Type, Dim>& v) {
		using RetT = vec_base<RetType, Dim>;
		RetT ret;
		for (int d = 0; d != Dim; ++d) {
			ret[d] = plus<RetType, Type>(v[d]);
		}
		return ret;
	}

	template<typename RetType, typename Type>
	auto plus(const mat4x4_t<Type>& v) {
		using RetT = mat4x4_t<RetType>;
		return RetT{
			plus<RetType, Type>(v[0]),
			plus<RetType, Type>(v[1]),
			plus<RetType, Type>(v[2]),
			plus<RetType, Type>(v[3])
		};
	}

	template<typename RetType, typename Type>
	auto plus(const mat3x3_t<Type>& v) {
		using RetT = mat3x3_t<RetType>;
		return RetT{
			plus<RetType, Type>(v[0]),
			plus<RetType, Type>(v[1]),
			plus<RetType, Type>(v[2])
		};
	}

	template<typename RetType, typename Type>
	auto plus(const mat2x2_t<Type>& v) {
		using RetT = mat2x2_t<RetType>;
		return RetT{
			plus<RetType, Type>(v[0]),
			plus<RetType, Type>(v[1])
		};
	}
	
	template<typename RetType, typename Type>
	auto plus(const vec4_t<Type>& v) {
		using RetT = vec4_t<RetType>;
		return RetT{
			plus<RetType, Type>(v[0]),
			plus<RetType, Type>(v[1]),
			plus<RetType, Type>(v[2]),
			plus<RetType, Type>(v[3])
		};
	}

	template<typename RetType, typename Type>
	auto plus(const vec3_t<Type>& v) {
		using RetT = vec3_t<RetType>;
		return RetT{
			plus<RetType, Type>(v[0]),
			plus<RetType, Type>(v[1]),
			plus<RetType, Type>(v[2])
		};
	}

	template<typename RetType, typename Type>
	auto plus(const vec2_t<Type>& v) {
		using RetT = vec2_t<RetType>;
		return RetT{
			plus<RetType, Type>(v[0]),
			plus<RetType, Type>(v[1])
		};
	}

	template<typename RetType, typename Type>
	auto plus(const quat_base<Type>& q) {
		using RetT = quat_base<RetType>;
		return RetT{
			plus<RetType, Type>(q[0]),
			plus<RetType, Type>(q[1]),
			plus<RetType, Type>(q[2]),
			plus<RetType, Type>(q[3])
		};
	}

	template<typename RetType, typename Type>
	auto plus(const complex_base<Type>& c) {
		using RetT = complex_base<RetType>;
		return RetT{
			plus<RetType, Type>(c[0]),
			plus<RetType, Type>(c[1])
		};
	}

#pragma endregion

//binary
#pragma region(add)
	template<typename RetType, typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol>
	auto add(const mat_base<LeftType, LeftRow, LeftCol>& lhs, const mat_base<RightType, RightRow, RightCol>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return mat_base<LeftType, LeftRow, LeftCol>();
	}

	template<typename RetType, typename LeftType, typename RightType, int LeftDim, int RightDim>
	auto add(const vec_base<LeftType, LeftDim>& lhs, const vec_base<RightType, RightDim>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return vec_base<LeftType, LeftDim>();
	}

	//n x m matrix
	template<typename RetType, typename LeftType, typename RightType, int Row, int Col>
	auto add(const mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		using RetT = mat_base<RetType, Row, Col>;
		RetT ret;
		for (int i = 0; i != Row; ++i) {
			ret[i] = add<RetType, LeftType, RightType, Col>(lhs[i], rhs[i]);
		}
		return ret;
	}

	//n dim vector 
	template<typename RetType, typename LeftType, typename RightType, int Dim>
	auto add(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		using RetT = vec_base<RetType, Dim>;
		RetT ret;
		for (int i = 0; i != Dim; ++i) {
			ret[i] = add<RetType, LeftType, RightType>(lhs[i], rhs[i]);
		}
		return ret;
	}

	//4 x 4 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto add(const mat4x4_t<LeftType>& lhs, const mat4x4_t<RightType>& rhs) {
		using RetT = mat4x4_t<RetType>;
		return RetT{
			add<RetType, LeftType, RightType>(lhs[0], rhs[0]),
			add<RetType, LeftType, RightType>(lhs[1], rhs[1]),
			add<RetType, LeftType, RightType>(lhs[2], rhs[2]),
			add<RetType, LeftType, RightType>(lhs[3], rhs[3])
		};
	}

	//3 x 3 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto add(const mat3x3_t<LeftType>& lhs, const mat3x3_t<RightType>& rhs) {
		using RetT = mat3x3_t<RetType>;
		return RetT{
			add<RetType, LeftType, RightType>(lhs[0] , rhs[0]),
			add<RetType, LeftType, RightType>(lhs[1] , rhs[1]),
			add<RetType, LeftType, RightType>(lhs[2] , rhs[2])
		};
	}

	//2 x 2 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto add(const mat2x2_t<LeftType>& lhs, const mat2x2_t<RightType>& rhs) {
		using RetT = mat2x2_t<RetType>;
		return RetT{
			add<RetType, LeftType, RightType>(lhs[0] , rhs[0]),
			add<RetType, LeftType, RightType>(lhs[1] , rhs[1])
		};
	}

	//4 dim vector 
	template<typename RetType, typename LeftType, typename RightType>
	auto add(const vec4_t<LeftType>& lhs, const vec4_t<RightType>& rhs) {
		using RetT = vec4_t<RetType>;
		return RetT{
			add<RetType, LeftType, RightType>(lhs[0], rhs[0]),
			add<RetType, LeftType, RightType>(lhs[1], rhs[1]),
			add<RetType, LeftType, RightType>(lhs[2], rhs[2]),
			add<RetType, LeftType, RightType>(lhs[3], rhs[3])
		};
	}

	//3 dim vector 
	template<typename RetType, typename LeftType, typename RightType>
	auto add(const vec3_t<LeftType>& lhs, const vec3_t<RightType>& rhs) {
		using RetT = vec3_t<RetType>;
		return RetT{
			add<RetType, LeftType, RightType>(lhs[0], rhs[0]),
			add<RetType, LeftType, RightType>(lhs[1], rhs[1]),
			add<RetType, LeftType, RightType>(lhs[2], rhs[2])
		};
	}

	//2 dim vector 
	template<typename RetType, typename LeftType, typename RightType>
	auto add(const vec2_t<LeftType>& lhs, const vec2_t<RightType>& rhs) {
		using RetT = vec2_t<RetType>;
		return RetT{
			add<RetType, LeftType, RightType>(lhs[0], rhs[0]),
			add<RetType, LeftType, RightType>(lhs[1], rhs[1])
		};
	}


	//quat
	template<typename RetType, typename LeftType, typename RightType>
	auto add(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		using RetT = quat_base<RetType>;
		return RetT{
			add<RetType, LeftType, RightType>(lhs[0], rhs[0]),
			add<RetType, LeftType, RightType>(lhs[1], rhs[1]),
			add<RetType, LeftType, RightType>(lhs[2], rhs[2]),
			add<RetType, LeftType, RightType>(lhs[3], rhs[3])
		};
	}

	//complex
	template<typename RetType, typename LeftType, typename RightType>
	auto add(const complex_base<LeftType>& lhs, const complex_base<RightType>& rhs) {
		using RetT = complex_base<RetType>;
		return RetT{
			add<RetType, LeftType, RightType>(lhs[0], rhs[0]),
			add<RetType, LeftType, RightType>(lhs[1], rhs[1])
		};
	}
#pragma endregion
	
//binary
#pragma region(sub)
	template<typename RetType, typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol>
	auto sub(const mat_base<LeftType, LeftRow, LeftCol>& lhs, const mat_base<RightType, RightRow, RightCol>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return mat_base<LeftType, LeftRow, LeftCol>();
	}

	template<typename RetType, typename LeftType, typename RightType, int LeftDim, int RightDim>
	auto sub(const vec_base<LeftType, LeftDim>& lhs, const vec_base<RightType, RightDim>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return vec_base<LeftType, LeftDim>();
	}

	//n x m matrix
	template<typename RetType, typename LeftType, typename RightType, int Row, int Col>
	auto sub(const mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		using RetT = mat_base<RetType, Row, Col>;
		RetT ret;
		for (int i = 0; i != Row; ++i) {
			ret[i] = sub<RetType, LeftType, RightType, Col>(lhs[i], rhs[i]);
		}
		return ret;
	}

	//n dim vector 
	template<typename RetType, typename LeftType, typename RightType, int Dim>
	auto sub(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		using RetT = vec_base<RetType, Dim>;
		RetT ret;
		for (int i = 0; i != Dim; ++i) {
			ret[i] = sub<RetType, LeftType, RightType>(lhs[i], rhs[i]);
		}
		return ret;
	}

	//4 x 4 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto sub(const mat4x4_t<LeftType>& lhs, const mat4x4_t<RightType>& rhs) {
		using RetT = mat4x4_t<RetType>;
		return RetT{
			sub<RetType, LeftType, RightType>(lhs[0], rhs[0]),
			sub<RetType, LeftType, RightType>(lhs[1], rhs[1]),
			sub<RetType, LeftType, RightType>(lhs[2], rhs[2]),
			sub<RetType, LeftType, RightType>(lhs[3], rhs[3])
		};
	}

	//3 x 3 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto sub(const mat3x3_t<LeftType>& lhs, const mat3x3_t<RightType>& rhs) {
		using RetT = mat3x3_t<RetType>;
		return RetT{
			sub<RetType, LeftType, RightType>(lhs[0] , rhs[0]),
			sub<RetType, LeftType, RightType>(lhs[1] , rhs[1]),
			sub<RetType, LeftType, RightType>(lhs[2] , rhs[2])
		};
	}

	//2 x 2 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto sub(const mat2x2_t<LeftType>& lhs, const mat2x2_t<RightType>& rhs) {
		using RetT = mat2x2_t<RetType>;
		return RetT{
			sub<RetType, LeftType, RightType>(lhs[0] , rhs[0]),
			sub<RetType, LeftType, RightType>(lhs[1] , rhs[1])
		};
	}

	//4 dim vector 
	template<typename RetType, typename LeftType, typename RightType>
	auto sub(const vec4_t<LeftType>& lhs, const vec4_t<RightType>& rhs) {
		using RetT = vec4_t<RetType>;
		return RetT{
			sub<RetType, LeftType, RightType>(lhs[0], rhs[0]),
			sub<RetType, LeftType, RightType>(lhs[1], rhs[1]),
			sub<RetType, LeftType, RightType>(lhs[2], rhs[2]),
			sub<RetType, LeftType, RightType>(lhs[3], rhs[3])
		};
	}

	//3 dim vector 
	template<typename RetType, typename LeftType, typename RightType>
	auto sub(const vec3_t<LeftType>& lhs, const vec3_t<RightType>& rhs) {
		using RetT = vec3_t<RetType>;
		return RetT{
			sub<RetType, LeftType, RightType>(lhs[0], rhs[0]),
			sub<RetType, LeftType, RightType>(lhs[1], rhs[1]),
			sub<RetType, LeftType, RightType>(lhs[2], rhs[2])
		};
	}

	//2 dim vector 
	template<typename RetType, typename LeftType, typename RightType>
	auto sub(const vec2_t<LeftType>& lhs, const vec2_t<RightType>& rhs) {
		using RetT = vec2_t<RetType>;
		return RetT{
			sub<RetType, LeftType, RightType>(lhs[0], rhs[0]),
			sub<RetType, LeftType, RightType>(lhs[1], rhs[1])
		};
	}


	//quat
	template<typename RetType, typename LeftType, typename RightType>
	auto sub(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		using RetT = quat_base<RetType>;
		return RetT{
			sub<RetType, LeftType, RightType>(lhs[0], rhs[0]),
			sub<RetType, LeftType, RightType>(lhs[1], rhs[1]),
			sub<RetType, LeftType, RightType>(lhs[2], rhs[2]),
			sub<RetType, LeftType, RightType>(lhs[3], rhs[3])
		};
	}

	//complex
	template<typename RetType, typename LeftType, typename RightType>
	auto sub(const complex_base<LeftType>& lhs, const complex_base<RightType>& rhs) {
		using RetT = complex_base<RetType>;
		return RetT{
			sub<RetType, LeftType, RightType>(lhs[0], rhs[0]),
			sub<RetType, LeftType, RightType>(lhs[1], rhs[1])
		};
	}
#pragma endregion
	
//binary
#pragma region(div)
	template<typename RetType, typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol>
	auto div(const mat_base<LeftType, LeftRow, LeftCol>& lhs, const mat_base<RightType, RightRow, RightCol>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return mat_base<LeftType, LeftRow, LeftCol>();
	}

	template<typename RetType, typename LeftType, typename RightType, int LeftDim, int RightDim>
	auto div(const vec_base<LeftType, LeftDim>& lhs, const vec_base<RightType, RightDim>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return vec_base<LeftType, LeftDim>();
	}

	//n x m matrix
	template<typename RetType, typename LeftType, typename RightType, int Row, int Col>
	auto div(const mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		using RetT = mat_base<RetType, Row, Col>;
		RetT ret;
		for (int r = 0; r != Row; ++r) {
			ret[r] = div<RetType, LeftType, RightType, Col>(lhs[r], rhs[r]);
		}
		return ret;
	}

	//n dim vector
	template<typename RetType, typename LeftType, typename RightType, int Dim>
	auto div(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		using RetT = vec_base<RetType, Dim>;
		RetT ret;
		for (int d = 0; d != Dim; ++d) {
			ret[d] = div<RetType>(lhs[d], rhs[d]);
		}
		return ret;
	}

	//4 x 4 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto div(const mat4x4_t<LeftType>& lhs, const mat4x4_t<RightType>& rhs) {
		using RetT = mat4x4_t<RetType>;
		return RetT{
			div<RetType>(lhs[0] , rhs[0]),
			div<RetType>(lhs[1] , rhs[1]),
			div<RetType>(lhs[2] , rhs[2]),
			div<RetType>(lhs[3] , rhs[3])
		};
	}

	//3 x 3 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto div(const mat3x3_t<LeftType>& lhs, const mat3x3_t<RightType>& rhs) {
		using RetT = mat3x3_t<RetType>;
		return RetT{
			div<RetType>(lhs[0] , rhs[0]),
			div<RetType>(lhs[1] , rhs[1]),
			div<RetType>(lhs[2] , rhs[2])
		};
	}

	//2 x 2 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto div(const mat2x2_t<LeftType>& lhs, const mat2x2_t<RightType>& rhs) {
		using RetT = mat2x2_t<RetType>;
		return RetT{
			div<RetType>(lhs[0] , rhs[0]),
			div<RetType>(lhs[1] , rhs[1])
		};
	}

	//4 dim vector
	template<typename RetType, typename LeftType, typename RightType>
	auto div(const vec4_t<LeftType>& lhs,  const vec4_t<RightType>& rhs) {
		using RetT = vec4_t<RetType>;
		return RetT{
			div<RetType>(lhs[0] , rhs[0]),
			div<RetType>(lhs[1] , rhs[1]),
			div<RetType>(lhs[2] , rhs[2]),
			div<RetType>(lhs[3] , rhs[3])
		};
	}

	//3 dim vector
	template<typename RetType, typename LeftType, typename RightType>
	auto div(const vec3_t<LeftType>& lhs, const vec3_t<RightType>& rhs) {
		using RetT = vec3_t<RetType>;
		return RetT{
			div<RetType>(lhs[0] , rhs[0]),
			div<RetType>(lhs[1] , rhs[1]),
			div<RetType>(lhs[2] , rhs[2])
		};
	}

	//2 dim vector
	template<typename RetType, typename LeftType, typename RightType>
	auto div(const vec2_t<LeftType>& lhs, const vec2_t<RightType>& rhs) {
		using RetT = vec2_t<RetType>;
		return RetT{
			div<RetType>(lhs[0] , rhs[0]),
			div<RetType>(lhs[1] , rhs[1])
		};
	}

	//quat
	template<typename RetType, typename LeftType, typename RightType>
	auto div(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		using RetT = quat_base<RetType>;
		return RetT{
			div<RetType>(lhs[0] , rhs[0]),
			div<RetType>(lhs[1] , rhs[1]),
			div<RetType>(lhs[2] , rhs[2]),
			div<RetType>(lhs[3] , rhs[3])
		};
	}

	//complex
	template<typename RetType, typename LeftType, typename RightType>
	auto div(const complex_base<LeftType>& lhs, const complex_base<RightType>& rhs) {
		using RetT = complex_base<RetType>;
		return RetT{
			div<RetType>(lhs[0] , rhs[0]),
			div<RetType>(lhs[1] , rhs[1])
		};
	}

	//left scalar
	//n x m matrix
	template<typename RetType, typename LeftType, typename RightType, int Row, int Col>
	auto div(const LeftType& lhs, const mat_base<RightType, Row, Col>& rhs) {
		using RetT = mat_base<RetType, Row, Col>;
		RetT ret;
		for (int r = 0; r != Row; ++r) {
			ret[r] = div<RetType, LeftType, RightType, Col>(lhs, rhs[r]);
		}
		return ret;
	}

	//left scalar
	//n dim vector
	template<typename RetType, typename LeftType, typename RightType, int Dim>
	auto div(const LeftType& lhs, const vec_base<RightType, Dim>& rhs) {
		using RetT = vec_base<RetType, Dim>;
		RetT ret;
		for (int d = 0; d != Dim; ++d) {
			ret[d] = div<RetType>(lhs, rhs[d]);
		}
		return ret;
	}

	//left scalar
	//4 x 4 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto div(const LeftType& lhs, const mat4x4_t<RightType>& rhs) {
		using RetT = mat4x4_t<RetType>;
		return RetT{
			div<RetType>(lhs, rhs[0]),
			div<RetType>(lhs, rhs[1]),
			div<RetType>(lhs, rhs[2]),
			div<RetType>(lhs, rhs[3])
		};
	}

	//left scalar
	//3 x 3 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto div(const LeftType& lhs, const mat3x3_t<RightType>& rhs) {
		using RetT = mat3x3_t<RetType>;
		return RetT{
			div<RetType>(lhs, rhs[0]),
			div<RetType>(lhs, rhs[1]),
			div<RetType>(lhs, rhs[2])
		};
	}

	//left scalar
	//2 x 2 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto div(const LeftType& lhs, const mat2x2_t<RightType>& rhs) {
		using RetT = mat2x2_t<RetType>;
		return RetT{
			div<RetType>(lhs, rhs[0]),
			div<RetType>(lhs, rhs[1])
		};
	}

	//left scalar
	//4 dim vector
	template<typename RetType, typename LeftType, typename RightType>
	auto div(const LeftType& lhs, const vec4_t<RightType>& rhs) {
		using RetT = vec4_t<RetType>;
		return RetT{
			div<RetType>(lhs, rhs[0]),
			div<RetType>(lhs, rhs[1]),
			div<RetType>(lhs, rhs[2]),
			div<RetType>(lhs, rhs[3])
		};
	}

	//left scalar
	//3 dim vector
	template<typename RetType, typename LeftType, typename RightType>
	auto div(const LeftType& lhs, const vec3_t<RightType>& rhs) {
		using RetT = vec3_t<RetType>;
		return RetT{
			div<RetType>(lhs, rhs[0]),
			div<RetType>(lhs, rhs[1]),
			div<RetType>(lhs, rhs[2])
		};
	}

	//left scalar
	//2 dim vector
	template<typename RetType, typename LeftType, typename RightType>
	auto div(const LeftType& lhs, const vec2_t<RightType>& rhs) {
		using RetT = vec2_t<RetType>;
		return RetT{
			div<RetType>(lhs, rhs[0]),
			div<RetType>(lhs, rhs[1])
		};
	}

	//left scalar
	//quat
	template<typename RetType, typename LeftType, typename RightType>
	auto div(const LeftType& lhs, const quat_base<RightType>& rhs) {
		using RetT = quat_base<RetType>;
		return RetT{
			div<RetType>(lhs, rhs[0]),
			div<RetType>(lhs, rhs[1]),
			div<RetType>(lhs, rhs[2]),
			div<RetType>(lhs, rhs[3])
		};
	}

	//left scalar
	//complex
	template<typename RetType, typename LeftType, typename RightType>
	auto div(const LeftType& lhs, const complex_base<RightType>& rhs) {
		using RetT = complex_base<RetType>;
		return RetT{
			div<RetType>(lhs, rhs[0]),
			div<RetType>(lhs, rhs[1])
		};
	}

	//right scalar
	//n x m matrix
	template<typename RetType, typename LeftType, typename RightType, int Row, int Col>
	auto div(const mat_base<LeftType, Row, Col>& lhs, const RightType& rhs) {
		using RetT = mat_base<RetType, Row, Col>;
		RetT ret;
		for (int r = 0; r != Row; ++r) {
			ret[r] = div<RetType, LeftType, RightType, Col>(lhs[r], rhs);
		}
		return ret;
	}

	//right scalar
	//n dim vector
	template<typename RetType, typename LeftType, typename RightType, int Dim>
	auto div(const vec_base<LeftType, Dim>& lhs, const RightType& rhs) {
		using RetT = vec_base<RetType, Dim>;
		RetT ret;
		for (int d = 0; d != Dim; ++d) {
			ret[d] = div<RetType>(lhs[d], rhs);
		}
		return ret;
	}

	//right scalar
	//4 x 4 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto div(const mat4x4_t<LeftType>& lhs, const RightType& rhs) {
		using RetT = mat4x4_t<RetType>;
		return RetT{
			div<RetType>(lhs[0], rhs),
			div<RetType>(lhs[1], rhs),
			div<RetType>(lhs[2], rhs),
			div<RetType>(lhs[3], rhs)
		};
	}

	//right scalar
	//3 x 3 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto div(const mat3x3_t<LeftType>& lhs, const RightType& rhs) {
		using RetT = mat3x3_t<RetType>;
		return RetT{
			div<RetType>(lhs[0], rhs),
			div<RetType>(lhs[1], rhs),
			div<RetType>(lhs[2], rhs)
		};
	}

	//right scalar
	//2 x 2 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto div(const mat2x2_t<LeftType>& lhs, const RightType& rhs) {
		using RetT = mat2x2_t<RetType>;
		return RetT{
			div<RetType>(lhs[0], rhs),
			div<RetType>(lhs[1], rhs)
		};
	}

	//right scalar
	//4 dim vector
	template<typename RetType, typename LeftType, typename RightType>
	auto div(const vec4_t<LeftType>& lhs, const RightType& rhs) {
		using RetT = vec4_t<RetType>;
		return RetT{
			div<RetType>(lhs[0], rhs),
			div<RetType>(lhs[1], rhs),
			div<RetType>(lhs[2], rhs),
			div<RetType>(lhs[3], rhs)
		};
	}

	//right scalar
	//3 dim vector
	template<typename RetType, typename LeftType, typename RightType>
	auto div(const vec3_t<LeftType>& lhs, const RightType& rhs) {
		using RetT = vec3_t<RetType>;
		return RetT{
			div<RetType>(lhs[0], rhs),
			div<RetType>(lhs[1], rhs),
			div<RetType>(lhs[2], rhs)
		};
	}

	//right scalar
	//2 dim vector
	template<typename RetType, typename LeftType, typename RightType>
	auto div(const vec2_t<LeftType>& lhs, const RightType& rhs) {
		using RetT = vec2_t<RetType>;
		return RetT{
			div<RetType>(lhs[0], rhs),
			div<RetType>(lhs[1], rhs)
		};
	}

	//right scalar
	//quat
	template<typename RetType, typename LeftType, typename RightType>
	auto div(const quat_base<LeftType>& lhs, const RightType& rhs) {
		using RetT = quat_base<LeftType>;
		return RetT{
			div<RetType>(lhs[0], rhs),
			div<RetType>(lhs[1], rhs),
			div<RetType>(lhs[2], rhs),
			div<RetType>(lhs[3], rhs)
		};
	}

	//right scalar
	//complex
	template<typename RetType, typename LeftType, typename RightType>
	auto div(const complex_base<LeftType>& lhs, const RightType& rhs) {
		using RetT = complex_base<LeftType>;
		return RetT{
			div<RetType>(lhs[0], rhs),
			div<RetType>(lhs[1], rhs)
		};
	}
#pragma endregion

//binary
#pragma region(mul)
namespace Details {
	template<bool b> struct mat_mul_mat_able {};

	//n x m matrix
	template<typename RetType, typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol>
	auto mul(const mat_base<LeftType, LeftRow, LeftCol>& lhs, const mat_base<RightType, RightRow, RightCol>& rhs, mat_mul_mat_able<true> unused) {
		using RetT = mat_base<RetType, RightCol, LeftRow>;
		RetT ret = RetT();

		for (int row = 0; row < LeftRow; ++row) {
			for (int i = 0; i < LeftCol; ++i) {
				for (int col = 0; col < LeftRow; ++col) {
					ret[row][col] += Math::mul<RetType, LeftType, RightType>(lhs[row][i], rhs[i][col]);
				}
			}
		}

		return ret;
	}

	//n x m matrix 
	template<typename RetType, typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol>
	auto mul(const mat_base<LeftType, LeftRow, LeftCol>& lhs, const mat_base<RightType, RightRow, RightCol>& rhs, mat_mul_mat_able<false> unused) {
		Z_MATH_ASSERT(0, __FUNCTION__ " (mat mul mat) : invalid mutiply able type.");
		return 0;
	}

}//namespace Details 

	//n x m matrix
	template<typename RetType, typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol>
	auto mul(const mat_base<LeftType, LeftRow, LeftCol>& lhs, const mat_base<RightType, RightRow, RightCol>& rhs) {
		return Details::mul<RetType, LeftType, RightType, LeftRow, LeftCol, RightRow, RightCol>(lhs, rhs, Details::mat_mul_mat_able<LeftCol == RightRow>());
	}

	//n dim vector error 
	template<typename LeftType, typename RightType, int LeftDim, int RightDim>
	auto mul(const vec_base<LeftType, LeftDim>& lhs, const vec_base<RightType, RightDim>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return vec_base<LeftType, LeftDim>();
	}

	//n dim vector 
	template<typename RetType, typename LeftType, typename RightType, int Dim>
	auto mul(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		using RetT = vec_base<RetType, Dim>;
		RetT ret;
		for (int i = 0; i != Dim; ++i) {
			ret[i] = mul<RetType, LeftType, RightType>(lhs[i], rhs[i]);
		}
		return ret;
	}

	//4 x 4 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto mul(const mat4x4_t<LeftType>& lhs, const mat4x4_t<RightType>& rhs) {
		using RetT = mat4x4_t<RetType>;
		using VecT = typename RetT::element_type;

		return RetT{ 
			VecT{
			mul<RetType, LeftType, RightType>(lhs[0][0] , rhs[0][0]) + mul<RetType, LeftType, RightType>(lhs[0][1] , rhs[1][0]) + mul<RetType, LeftType, RightType>(lhs[0][2] , rhs[2][0]) + mul<RetType, LeftType, RightType>(lhs[0][3] , rhs[3][0]) ,
			mul<RetType, LeftType, RightType>(lhs[0][0] , rhs[0][1]) + mul<RetType, LeftType, RightType>(lhs[0][1] , rhs[1][1]) + mul<RetType, LeftType, RightType>(lhs[0][2] , rhs[2][1]) + mul<RetType, LeftType, RightType>(lhs[0][3] , rhs[3][1]) ,
			mul<RetType, LeftType, RightType>(lhs[0][0] , rhs[0][2]) + mul<RetType, LeftType, RightType>(lhs[0][1] , rhs[1][2]) + mul<RetType, LeftType, RightType>(lhs[0][2] , rhs[2][2]) + mul<RetType, LeftType, RightType>(lhs[0][3] , rhs[3][2]) ,
			mul<RetType, LeftType, RightType>(lhs[0][0] , rhs[0][3]) + mul<RetType, LeftType, RightType>(lhs[0][1] , rhs[1][3]) + mul<RetType, LeftType, RightType>(lhs[0][2] , rhs[2][3]) + mul<RetType, LeftType, RightType>(lhs[0][3] , rhs[3][3])},					
			VecT{
			mul<RetType, LeftType, RightType>(lhs[1][0] , rhs[0][0]) + mul<RetType, LeftType, RightType>(lhs[1][1] , rhs[1][0]) + mul<RetType, LeftType, RightType>(lhs[1][2] , rhs[2][0]) + mul<RetType, LeftType, RightType>(lhs[1][3] , rhs[3][0]) ,
			mul<RetType, LeftType, RightType>(lhs[1][0] , rhs[0][1]) + mul<RetType, LeftType, RightType>(lhs[1][1] , rhs[1][1]) + mul<RetType, LeftType, RightType>(lhs[1][2] , rhs[2][1]) + mul<RetType, LeftType, RightType>(lhs[1][3] , rhs[3][1]) ,
			mul<RetType, LeftType, RightType>(lhs[1][0] , rhs[0][2]) + mul<RetType, LeftType, RightType>(lhs[1][1] , rhs[1][2]) + mul<RetType, LeftType, RightType>(lhs[1][2] , rhs[2][2]) + mul<RetType, LeftType, RightType>(lhs[1][3] , rhs[3][2]) ,
			mul<RetType, LeftType, RightType>(lhs[1][0] , rhs[0][3]) + mul<RetType, LeftType, RightType>(lhs[1][1] , rhs[1][3]) + mul<RetType, LeftType, RightType>(lhs[1][2] , rhs[2][3]) + mul<RetType, LeftType, RightType>(lhs[1][3] , rhs[3][3])},					
			VecT{
			mul<RetType, LeftType, RightType>(lhs[2][0] , rhs[0][0]) + mul<RetType, LeftType, RightType>(lhs[2][1] , rhs[1][0]) + mul<RetType, LeftType, RightType>(lhs[2][2] , rhs[2][0]) + mul<RetType, LeftType, RightType>(lhs[2][3] , rhs[3][0]) ,
			mul<RetType, LeftType, RightType>(lhs[2][0] , rhs[0][1]) + mul<RetType, LeftType, RightType>(lhs[2][1] , rhs[1][1]) + mul<RetType, LeftType, RightType>(lhs[2][2] , rhs[2][1]) + mul<RetType, LeftType, RightType>(lhs[2][3] , rhs[3][1]) ,
			mul<RetType, LeftType, RightType>(lhs[2][0] , rhs[0][2]) + mul<RetType, LeftType, RightType>(lhs[2][1] , rhs[1][2]) + mul<RetType, LeftType, RightType>(lhs[2][2] , rhs[2][2]) + mul<RetType, LeftType, RightType>(lhs[2][3] , rhs[3][2]) ,
			mul<RetType, LeftType, RightType>(lhs[2][0] , rhs[0][3]) + mul<RetType, LeftType, RightType>(lhs[2][1] , rhs[1][3]) + mul<RetType, LeftType, RightType>(lhs[2][2] , rhs[2][3]) + mul<RetType, LeftType, RightType>(lhs[2][3] , rhs[3][3])},					
			VecT{
			mul<RetType, LeftType, RightType>(lhs[3][0] , rhs[0][0]) + mul<RetType, LeftType, RightType>(lhs[3][1] , rhs[1][0]) + mul<RetType, LeftType, RightType>(lhs[3][2] , rhs[2][0]) + mul<RetType, LeftType, RightType>(lhs[3][3] , rhs[3][0]) ,
			mul<RetType, LeftType, RightType>(lhs[3][0] , rhs[0][1]) + mul<RetType, LeftType, RightType>(lhs[3][1] , rhs[1][1]) + mul<RetType, LeftType, RightType>(lhs[3][2] , rhs[2][1]) + mul<RetType, LeftType, RightType>(lhs[3][3] , rhs[3][1]) ,
			mul<RetType, LeftType, RightType>(lhs[3][0] , rhs[0][2]) + mul<RetType, LeftType, RightType>(lhs[3][1] , rhs[1][2]) + mul<RetType, LeftType, RightType>(lhs[3][2] , rhs[2][2]) + mul<RetType, LeftType, RightType>(lhs[3][3] , rhs[3][2]) ,
			mul<RetType, LeftType, RightType>(lhs[3][0] , rhs[0][3]) + mul<RetType, LeftType, RightType>(lhs[3][1] , rhs[1][3]) + mul<RetType, LeftType, RightType>(lhs[3][2] , rhs[2][3]) + mul<RetType, LeftType, RightType>(lhs[3][3] , rhs[3][3])}
		};
	}

	//3 x 3 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto mul(const mat3x3_t<LeftType>& lhs, const mat3x3_t<RightType>& rhs) {
		using RetT = mat3x3_t<RetType>;
		using VecT = typename RetT::element_type;

		return RetT{ 
			VecT{
			mul<RetType, LeftType, RightType>(lhs[0][0], rhs[0][0]) + mul<RetType, LeftType, RightType>(lhs[0][1], rhs[1][0]) + mul<RetType, LeftType, RightType>(lhs[0][2], rhs[2][0]) ,
			mul<RetType, LeftType, RightType>(lhs[0][0], rhs[0][1]) + mul<RetType, LeftType, RightType>(lhs[0][1], rhs[1][1]) + mul<RetType, LeftType, RightType>(lhs[0][2], rhs[2][1]) ,
			mul<RetType, LeftType, RightType>(lhs[0][0], rhs[0][2]) + mul<RetType, LeftType, RightType>(lhs[0][1], rhs[1][2]) + mul<RetType, LeftType, RightType>(lhs[0][2], rhs[2][2])},				VecT{
			mul<RetType, LeftType, RightType>(lhs[1][0], rhs[0][0]) + mul<RetType, LeftType, RightType>(lhs[1][1], rhs[1][0]) + mul<RetType, LeftType, RightType>(lhs[1][2], rhs[2][0]) ,
			mul<RetType, LeftType, RightType>(lhs[1][0], rhs[0][1]) + mul<RetType, LeftType, RightType>(lhs[1][1], rhs[1][1]) + mul<RetType, LeftType, RightType>(lhs[1][2], rhs[2][1]) ,
			mul<RetType, LeftType, RightType>(lhs[1][0], rhs[0][2]) + mul<RetType, LeftType, RightType>(lhs[1][1], rhs[1][2]) + mul<RetType, LeftType, RightType>(lhs[1][2], rhs[2][2])},				VecT{
			mul<RetType, LeftType, RightType>(lhs[2][0], rhs[0][0]) + mul<RetType, LeftType, RightType>(lhs[2][1], rhs[1][0]) + mul<RetType, LeftType, RightType>(lhs[2][2], rhs[2][0]) ,
			mul<RetType, LeftType, RightType>(lhs[2][0], rhs[0][1]) + mul<RetType, LeftType, RightType>(lhs[2][1], rhs[1][1]) + mul<RetType, LeftType, RightType>(lhs[2][2], rhs[2][1]) ,
			mul<RetType, LeftType, RightType>(lhs[2][0], rhs[0][2]) + mul<RetType, LeftType, RightType>(lhs[2][1], rhs[1][2]) + mul<RetType, LeftType, RightType>(lhs[2][2], rhs[2][2])}
		};
	}

	//2 x 2 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto mul(const mat2x2_t<LeftType>& lhs, const mat2x2_t<RightType>& rhs) {
		using RetT = mat2x2_t<RetType>;
		using VecT = typename RetT::element_type;

		return RetT{
			VecT{
			mul<RetType, LeftType, RightType>(lhs[0][0],rhs[0][0]) + mul<RetType, LeftType, RightType>(lhs[0][1], rhs[1][0]),
			mul<RetType, LeftType, RightType>(lhs[0][0],rhs[0][1]) + mul<RetType, LeftType, RightType>(lhs[0][1], rhs[1][1])},

			VecT{
			mul<RetType, LeftType, RightType>(lhs[1][0],rhs[0][0]) + mul<RetType, LeftType, RightType>(lhs[1][1], rhs[1][0]) ,
			mul<RetType, LeftType, RightType>(lhs[1][0],rhs[0][1]) + mul<RetType, LeftType, RightType>(lhs[1][1], rhs[1][1])}
		};
	}

	//4 dim vector 
	template<typename RetType, typename LeftType, typename RightType>
	auto mul(const vec4_t<LeftType>& lhs, const vec4_t<RightType>& rhs) {
		using RetT = vec4_t<RetType>;
		return RetT{ mul<RetType, LeftType, RightType>(lhs[0], rhs[0]), mul<RetType, LeftType, RightType>(lhs[1], rhs[1]), mul<RetType, LeftType, RightType>(lhs[2], rhs[2]), mul<RetType, LeftType, RightType>(lhs[3], rhs[3]) };
	}

	//3 dim vector 
	template<typename RetType, typename LeftType, typename RightType>
	auto mul(const vec3_t<LeftType>& lhs, const vec3_t<RightType>& rhs) {
		using RetT = vec3_t<RetType>;
		return RetT{ mul<RetType, LeftType, RightType>(lhs[0], rhs[0]), mul<RetType, LeftType, RightType>(lhs[1], rhs[1]), mul<RetType, LeftType, RightType>(lhs[2], rhs[2]) };
	}

	//2 dim vector 
	template<typename RetType, typename LeftType, typename RightType>
	auto mul(const vec2_t<LeftType>& lhs, const vec2_t<RightType>& rhs) {
		using RetT = vec2_t<RetType>;
		return RetT{ mul<RetType, LeftType, RightType>(lhs[0], rhs[0]), mul<RetType, LeftType, RightType>(lhs[1], rhs[1]) };
	}

	//quat
	// Returns the product rhs*lhs (which is the concatenation of a rotation lhs followed by the rotation rhs)
	template<typename RetType, typename LeftType, typename RightType>
	auto mul(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		using RetT = quat_base<RetType>;
		return RetT{
			mul<RetType, RightType, LeftType>(rhs[3],lhs[0]) + mul<RetType, RightType, LeftType>(rhs[0],lhs[3]) + mul<RetType, RightType, LeftType>(rhs[1],lhs[2]) - mul<RetType, RightType, LeftType>(rhs[2],lhs[1]),
			mul<RetType, RightType, LeftType>(rhs[3],lhs[1]) - mul<RetType, RightType, LeftType>(rhs[0],lhs[2]) + mul<RetType, RightType, LeftType>(rhs[1],lhs[3]) + mul<RetType, RightType, LeftType>(rhs[2],lhs[0]),
			mul<RetType, RightType, LeftType>(rhs[3],lhs[2]) + mul<RetType, RightType, LeftType>(rhs[0],lhs[1]) - mul<RetType, RightType, LeftType>(rhs[1],lhs[0]) + mul<RetType, RightType, LeftType>(rhs[2],lhs[3]),
			mul<RetType, RightType, LeftType>(rhs[3],lhs[3]) - mul<RetType, RightType, LeftType>(rhs[0],lhs[0]) - mul<RetType, RightType, LeftType>(rhs[1],lhs[1]) - mul<RetType, RightType, LeftType>(rhs[2],lhs[2])
		};
	}

	//complex
	template<typename RetType, typename LeftType, typename RightType>
	auto mul(const complex_base<LeftType>& lhs, const complex_base<RightType>& rhs) {
		using RetT = complex_base<RetType>;
		return RetT{
			mul<RetType, LeftType, RightType>(lhs[0] * rhs[0]) - mul<RetType, LeftType, RightType>(lhs[1] * rhs[1]),
			mul<RetType, LeftType, RightType>(lhs[0] * rhs[1]) + mul<RetType, LeftType, RightType>(lhs[1] * rhs[0])
		};
	}

namespace Details {
	template<bool b> struct mat_mul_vec_able {};

	//n x m matrix * vec n t
	template<typename RetType, typename MatType, typename VecType, int Row, int Col, int Dim>
	auto mul(const mat_base<MatType, Row, Col>& m, const vec_base<VecType, Dim>& v, mat_mul_vec_able<true> unused) {
		using RetT = vec_base<RetType, Row>;
		RetT ret = RetT();

		for (int d = 0; d != Row; ++d) {
			for (int r = 0; r != Col; ++r) {
				ret[d] += Math::mul<RetType, VecType, MatType>(m[d][r], v[r]);
			}
		}

		return ret;
	}

	//n x m matrix * vec n t
	template<typename RetType, typename MatType, typename VecType, int Row, int Col, int Dim>
	auto mul(const mat_base<MatType, Row, Col>& lhs, const vec_base<VecType, Dim>& rhs, mat_mul_vec_able<false> unused) {
		Z_MATH_ASSERT(0, __FUNCTION__ " (mat mul vec) : invalid mutiply able type.");
		return 0;
	}
}//namespace Details

	//n x m matrix * vec n t
	template<typename RetType, typename MatType, typename VecType, int Row, int Col, int Dim>
	auto mul(const mat_base<MatType, Row, Col>& m, const vec_base<VecType, Dim>& v) {
		return Details::mul<RetType, MatType, VecType, Row, Col, Dim>(m, v, Details::mat_mul_vec_able<Col == Dim>());
	}

	//4 x 4 matrix * vec 4t
	template<typename RetType, typename MatType, typename VecType>
	auto mul(const mat4x4_t<MatType>& m, const vec4_t<VecType>& v) {
		using RetT = vec4_t<RetType>;

		return RetT{
			mul<RetType, VecType, MatType>(v[0] , m[0][0]) + mul<RetType, VecType, MatType>(v[1] , m[0][1]) + mul<RetType, VecType, MatType>(v[2] , m[0][2]) + mul<RetType, VecType, MatType>(v[3] , m[0][3]) ,
			mul<RetType, VecType, MatType>(v[0] , m[1][0]) + mul<RetType, VecType, MatType>(v[1] , m[1][1]) + mul<RetType, VecType, MatType>(v[2] , m[1][2]) + mul<RetType, VecType, MatType>(v[3] , m[1][3]) ,
			mul<RetType, VecType, MatType>(v[0] , m[2][0]) + mul<RetType, VecType, MatType>(v[1] , m[2][1]) + mul<RetType, VecType, MatType>(v[2] , m[2][2]) + mul<RetType, VecType, MatType>(v[3] , m[2][3]) ,
			mul<RetType, VecType, MatType>(v[0] , m[3][0]) + mul<RetType, VecType, MatType>(v[1] , m[3][1]) + mul<RetType, VecType, MatType>(v[2] , m[3][2]) + mul<RetType, VecType, MatType>(v[3] , m[3][3])
		};
	}

	//3 x 3 matrix * vec 3t
	template<typename RetType, typename MatType, typename VecType>
	auto mul(const mat3x3_t<MatType>& m, const vec3_t<VecType>& v) {
		using RetT = vec3_t<RetType>;

		return RetT{
			mul<RetType, VecType, MatType>(v[0] , m[0][0]) + mul<RetType, VecType, MatType>(v[1] , m[0][1]) + mul<RetType, VecType, MatType>(v[2] , m[0][2]) ,
			mul<RetType, VecType, MatType>(v[0] , m[1][0]) + mul<RetType, VecType, MatType>(v[1] , m[1][1]) + mul<RetType, VecType, MatType>(v[2] , m[1][2]) ,
			mul<RetType, VecType, MatType>(v[0] , m[2][0]) + mul<RetType, VecType, MatType>(v[1] , m[2][1]) + mul<RetType, VecType, MatType>(v[2] , m[2][2]) 
		};
	}

	//2 x 2 matrix * vec 2t
	template<typename RetType, typename MatType, typename VecType>
	auto mul(const mat2x2_t<MatType>& m, const vec2_t<VecType>& v) {
		using RetT = vec2_t<RetType>;

		return RetT{
			mul<RetType, VecType, MatType>(v[0] , m[0][0]) + mul<RetType, VecType, MatType>(v[1] , m[0][1]) ,
			mul<RetType, VecType, MatType>(v[0] , m[1][0]) + mul<RetType, VecType, MatType>(v[1] , m[1][1]) 
		};
	}

namespace Details {
	template<bool b> struct vec_mul_mat_able {};

	// vec n t * m x n matrix
	template<typename RetType, typename MatType, typename VecType, int Row, int Col, int Dim>
	auto mul(const vec_base<VecType, Dim>& lhs, const mat_base<MatType, Row, Col>& rhs, mat_mul_vec_able<true> unused) {
		using RetT = vec_base<RetType, Col>;
		RetT ret = RetT();

		for (int d = 0; d != Col; ++d) {
			for (int r = 0; r != Row; ++r) {
				ret[d] += Math::mul<RetType, VecType, MatType>(lhs[r], rhs[r][d]);
			}
		}

		return ret;
	}

	// vec n t * m x n matrix
	template<typename RetType, typename MatType, typename VecType, int Row, int Col, int Dim>
	auto mul(const vec_base<VecType, Dim>& lhs, const mat_base<MatType, Row, Col>& rhs, vec_mul_mat_able<false> unused) {
		Z_MATH_ASSERT(0, __FUNCTION__ " (vec mul mat) : invalid mutiply able type.");
		return 0;
	}
}//namespace Details 

	// vec n t * m x b matrix
	template<typename MatType, typename VecType, int Row, int Col, int Dim>
	auto mul(const vec_base<VecType, Dim>& lhs, const mat_base<MatType, Row, Col>& rhs) {
		return Details::mul(lhs, rhs, Details::vec_mul_mat_able<Row == Dim>());
	}

	//vec 4t * 4 x 4 matrix
	template<typename RetType, typename VecType, typename MatType>
	auto mul(const vec4_t<VecType>& v, const mat4x4_t<MatType>& m) {
		using RetT = vec4_t<RetType>;

		return RetT{
			mul<RetType, VecType, MatType>(v[0], m[0][0]) + mul<RetType, VecType, MatType>(v[1], m[1][0]) + mul<RetType, VecType, MatType>(v[2], m[2][0]) + mul<RetType, VecType, MatType>(v[3], m[3][0]) ,
			mul<RetType, VecType, MatType>(v[0], m[0][1]) + mul<RetType, VecType, MatType>(v[1], m[1][1]) + mul<RetType, VecType, MatType>(v[2], m[2][1]) + mul<RetType, VecType, MatType>(v[3], m[3][1]) ,
			mul<RetType, VecType, MatType>(v[0], m[0][2]) + mul<RetType, VecType, MatType>(v[1], m[1][2]) + mul<RetType, VecType, MatType>(v[2], m[2][2]) + mul<RetType, VecType, MatType>(v[3], m[3][2]) ,
			mul<RetType, VecType, MatType>(v[0], m[0][3]) + mul<RetType, VecType, MatType>(v[1], m[1][3]) + mul<RetType, VecType, MatType>(v[2], m[2][3]) + mul<RetType, VecType, MatType>(v[3], m[3][3])
		};
	}

	//vec 3t * 3 x 3 matrix
	template<typename RetType, typename VecType, typename MatType>
	auto mul(const vec3_t<VecType>& v, const mat3x3_t<MatType>& m) {
		using RetT = vec3_t<RetType>;

		return RetT{
			mul<RetType, VecType, MatType>(v[0], m[0][0]) + mul<RetType, VecType, MatType>(v[1], m[1][0]) + mul<RetType, VecType, MatType>(v[2], m[2][0]) ,
			mul<RetType, VecType, MatType>(v[0], m[0][1]) + mul<RetType, VecType, MatType>(v[1], m[1][1]) + mul<RetType, VecType, MatType>(v[2], m[2][1]) ,
			mul<RetType, VecType, MatType>(v[0], m[0][2]) + mul<RetType, VecType, MatType>(v[1], m[1][2]) + mul<RetType, VecType, MatType>(v[2], m[2][2]) 
		};
	}

	//vec 2t * 2 x 2 matrix
	template<typename RetType, typename VecType, typename MatType>
	auto mul(const vec2_t<VecType>& v, const mat2x2_t<MatType>& m) {
		using RetT = vec2_t<RetType>;

		return RetT{
			mul<RetType, VecType, MatType>(v[0], m[0][0]) + mul<RetType, VecType, MatType>(v[1], m[1][0]) ,
			mul<RetType, VecType, MatType>(v[0], m[0][1]) + mul<RetType, VecType, MatType>(v[1], m[1][1]) 
		};
	}

	//n x m matrix * scalar
	template<typename RetType, typename LeftType, typename RightType , int Row, int Col>
	auto mul(const mat_base<LeftType, Row, Col>& m , const RightType& s) {
		static_assert(is_arithmetic<RightType>::value, __FUNCTION__ " : must arithmetic type.");
		using RetT = mat_base<RetType, Row, Col>;
		RetT ret;
		for (int row = 0; row != Row; ++row) {
			ret[row] = mul<RetType, LeftType, RightType, Col>(m[row], s);
		}
		return ret;
	}

	//scalar * n x m matrix 
	template<typename RetType, typename LeftType, typename RightType, int Row, int Col>
	auto mul(const LeftType& s, const mat_base<RightType, Row, Col>& m) {
		static_assert(is_arithmetic<LeftType>::value, __FUNCTION__ " : must arithmetic type.");
		return mul<RetType, RightType, LeftType, Row, Col>(m, s);
	}

	//vec nt * scalar
	template<typename RetType, typename LeftType, typename RightType, int Dim>
	auto mul(const vec_base<LeftType, Dim>& v, const RightType& s) {
		static_assert(is_arithmetic<RightType>::value, __FUNCTION__ " : must arithmetic type.");
		using RetT = vec_base<RetType, Dim>;
		RetT ret;
		for (int d = 0; d != Dim; ++d) {
			ret[d] = mul<RetType, LeftType, RightType>(v[d], s);
		}
		return ret;
	}

	//scalar * vec nt 
	template<typename RetType, typename LeftType, typename RightType, int Dim>
	auto mul(const LeftType& s, const vec_base<RightType, Dim>& v) {
		static_assert(is_arithmetic<LeftType>::value, __FUNCTION__ " : must arithmetic type.");
		return mul<RetType, RightType, LeftType, Dim>(v, s);
	}

	//vec 4t * scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto mul(const vec4_t<LeftType>& v, const RightType& s) {
		static_assert(is_arithmetic<RightType>::value, __FUNCTION__ " : must arithmetic type.");
		using RetT = vec4_t<RetType>;
		return RetT{
			mul<RetType, LeftType, RightType>(v[0], s), 
			mul<RetType, LeftType, RightType>(v[1], s), 
			mul<RetType, LeftType, RightType>(v[2], s), 
			mul<RetType, LeftType, RightType>(v[3], s)
		};
	}
	
	//vec 3t * scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto mul(const vec3_t<LeftType>& v, const RightType& s) {
		static_assert(is_arithmetic<RightType>::value, __FUNCTION__ " : must arithmetic type.");
		using RetT = vec3_t<RetType>;
		return RetT{
			mul<RetType, LeftType, RightType>(v[0], s),
			mul<RetType, LeftType, RightType>(v[1], s),
			mul<RetType, LeftType, RightType>(v[2], s)
		};
	}
	
	//vec 2t * scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto mul(const vec2_t<LeftType>& v, const RightType& s) {
		static_assert(is_arithmetic<RightType>::value, __FUNCTION__ " : must arithmetic type.");
		using RetT = vec2_t<RetType>;
		return RetT{
			mul<RetType, LeftType, RightType>(v[0], s),
			mul<RetType, LeftType, RightType>(v[1], s)
		};
	}

	//scalar * vec 4t
	template<typename RetType, typename LeftType, typename RightType>
	auto mul(const LeftType& s, const vec4_t<RightType>& v) {
		static_assert(is_arithmetic<LeftType>::value, __FUNCTION__ " : must arithmetic type.");
		using RetT = vec4_t<RetType>;
		return RetT{
			mul<RetType, LeftType, RightType>(s, v[0]),
			mul<RetType, LeftType, RightType>(s, v[1]),
			mul<RetType, LeftType, RightType>(s, v[2]),
			mul<RetType, LeftType, RightType>(s, v[3])
		};
	}

	//scalar * vec 3t
	template<typename RetType, typename LeftType, typename RightType>
	auto mul(const LeftType& s, const vec3_t<RightType>& v) {
		static_assert(is_arithmetic<LeftType>::value, __FUNCTION__ " : must arithmetic type.");
		using RetT = vec3_t<RetType>;
		return RetT{
			mul<RetType, LeftType, RightType>(s, v[0]),
			mul<RetType, LeftType, RightType>(s, v[1]),
			mul<RetType, LeftType, RightType>(s, v[2])
		};
	}

	//scalar * vec 2t
	template<typename RetType, typename LeftType, typename RightType>
	auto mul(const LeftType& s, const vec2_t<RightType>& v) {
		static_assert(is_arithmetic<LeftType>::value, __FUNCTION__ " : must arithmetic type.");
		using RetT = vec2_t<RetType>;
		return RetT{
			mul<RetType, LeftType, RightType>(s, v[0]),
			mul<RetType, LeftType, RightType>(s, v[1])
		};
	}

	//quat * scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto mul(const quat_base<LeftType>& q, const RightType& s) {
		static_assert(is_arithmetic<RightType>::value, __FUNCTION__ " : must arithmetic type.");
		using RetT = quat_base<RetType>;
		return RetT{
			mul<RetType, LeftType, RightType>(q[0], s),
			mul<RetType, LeftType, RightType>(q[1], s),
			mul<RetType, LeftType, RightType>(q[2], s),
			mul<RetType, LeftType, RightType>(q[3], s)
		};
	}

	//scalar * quat 
	template<typename RetType, typename LeftType, typename RightType>
	auto mul(const LeftType& s, const quat_base<RightType>& q) {
		static_assert(is_arithmetic<LeftType>::value, __FUNCTION__ " : must arithmetic type.");
		using RetT = quat_base<RetType>;
		return RetT{
			mul<RetType, LeftType, RightType>(s, q[0]),
			mul<RetType, LeftType, RightType>(s, q[1]),
			mul<RetType, LeftType, RightType>(s, q[2]),
			mul<RetType, LeftType, RightType>(s, q[3])
		};
	}

	//complex * scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto mul(const complex_base<LeftType>& c, const RightType& s) {
		static_assert(is_arithmetic<RightType>::value, __FUNCTION__ " : must arithmetic type.");
		using RetT = complex_base<RetType>;
		return RetT{
			mul<RetType, LeftType, RightType>(c[0], s),
			mul<RetType, LeftType, RightType>(c[1], s)
		};
	}

	//scalar * complex 
	template<typename RetType, typename LeftType, typename RightType>
	auto mul(const LeftType& s, const complex_base<RightType>& c) {
		static_assert(is_arithmetic<LeftType>::value, __FUNCTION__ " : must arithmetic type.");
		using RetT = complex_base<RetType>;
		return RetT{
			mul<RetType, LeftType, RightType>(s, c[0]),
			mul<RetType, LeftType, RightType>(s, c[1])
		};
	}


#pragma endregion

//binary
#pragma region(dot)
	//vec nt dot vec nt
	template<typename RetType, typename LeftType, typename RightType, int Dim>
	RetType dot(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		RetType ret = RetType();
		for (int d = 0; d != Dim; ++d) {
			ret += mul<RetType, LeftType, RightType>(lhs[d], rhs[d]);
		}
		return ret;
	}

	//vec 4t dot vec 4t
	template<typename RetType, typename LeftType, typename RightType>
	RetType dot(const vec4_t<LeftType>& lhs, const vec4_t<RightType>& rhs) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return RetType{
			mul<RetType, LeftType, RightType>(lhs[0], rhs[0]) +
			mul<RetType, LeftType, RightType>(lhs[1], rhs[1]) +
			mul<RetType, LeftType, RightType>(lhs[2], rhs[2]) +
			mul<RetType, LeftType, RightType>(lhs[3], rhs[3])
		};
	}

	//vec 4t dot vec 4t
	template<typename RetType, typename Type>
	RetType dot(const vec4_t<Type>& lhs, const vec4_t<Type>& rhs) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return RetType{
			mul<RetType, Type, Type>(lhs[0], rhs[0]) +
			mul<RetType, Type, Type>(lhs[1], rhs[1]) +
			mul<RetType, Type, Type>(lhs[2], rhs[2]) +
			mul<RetType, Type, Type>(lhs[3], rhs[3])
		};
	}

	//vec 3t dot vec 3t
	template<typename RetType, typename LeftType, typename RightType>
	RetType dot(const vec3_t<LeftType>& lhs, const vec3_t<RightType>& rhs) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return RetType{
			mul<RetType, LeftType, RightType>(lhs[0], rhs[0]) +
			mul<RetType, LeftType, RightType>(lhs[1], rhs[1]) +
			mul<RetType, LeftType, RightType>(lhs[2], rhs[2])
		};
	}

	//vec 2t dot vec 2t
	template<typename RetType, typename LeftType, typename RightType>
	RetType dot(const vec2_t<LeftType>& lhs, const vec2_t<RightType>& rhs) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return RetType{
			mul<RetType, LeftType, RightType>(lhs[0], rhs[0]) +
			mul<RetType, LeftType, RightType>(lhs[1], rhs[1])
		};
	}

	//quat dot quat
	template<typename RetType, typename LeftType, typename RightType>
	RetType dot(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return RetType{
			mul<RetType, LeftType, RightType>(lhs[0], rhs[0]) +
			mul<RetType, LeftType, RightType>(lhs[1], rhs[1]) +
			mul<RetType, LeftType, RightType>(lhs[2], rhs[2]) +
			mul<RetType, LeftType, RightType>(lhs[3], rhs[3])
		};
	}

	//complex dot complex
	template<typename RetType, typename LeftType, typename RightType>
	RetType dot(const complex_base<LeftType>& lhs, const complex_base<RightType>& rhs) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return RetType{
			mul<RetType, LeftType, RightType>(lhs[0], rhs[0])+
			mul<RetType, LeftType, RightType>(lhs[1], rhs[1])
		};
	}
#pragma endregion

#pragma region(lengthSq)
	template<typename RetType, typename Type, int Dim>
	RetType lengthSq(const vec_base<Type, Dim>& v) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		RetType ret = RetType();
		for (int d = 0; d != Dim; ++d) {
			ret += mul<RetType, Type, Type>(v[d], v[d]);
		}
		return ret;
	}

	template<typename RetType, typename Type>
	RetType lengthSq(const vec4_t<Type>& v) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return RetType(
			mul<RetType, Type, Type>(v[0], v[0]) +
			mul<RetType, Type, Type>(v[1], v[1]) +
			mul<RetType, Type, Type>(v[2], v[2]) +
			mul<RetType, Type, Type>(v[3], v[3])
		);
	}

	template<typename RetType, typename Type>
	RetType lengthSq(const vec3_t<Type>& v) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return RetType(
			mul<RetType, Type, Type>(v[0], v[0]) +
			mul<RetType, Type, Type>(v[1], v[1]) +
			mul<RetType, Type, Type>(v[2], v[2])
		);
	}

	template<typename RetType, typename Type>
	RetType lengthSq(const vec2_t<Type>& v) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return RetType(
			mul<RetType, Type, Type>(v[0], v[0]) +
			mul<RetType, Type, Type>(v[1], v[1])
		);
	}

	template<typename RetType, typename Type>
	RetType lengthSq(const quat_base<Type>& q) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return RetType(
			mul<RetType,Type,Type>(q[0],q[0]) +
			mul<RetType,Type,Type>(q[1],q[1]) +
			mul<RetType,Type,Type>(q[2],q[2]) +
			mul<RetType,Type,Type>(q[3],q[3]) 
		);
	}

	template<typename RetType, typename Type>
	RetType lengthSq(const complex_base<Type>& c) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return RetType(
			mul<RetType, Type, Type>(c[0], c[0]) +
			mul<RetType, Type, Type>(c[0], c[0])
		);
	}
#pragma endregion

#pragma region(length)
	template<typename RetType, typename Type, int Dim>
	RetType length(const vec_base<Type, Dim>& v) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return static_cast<RetType>(std::sqrt(lengthSq<RetType, Type, Dim>(v)));
	}

	template<typename RetType, typename Type>
	RetType length(const vec4_t<Type>& v) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return static_cast<RetType>(std::sqrt(lengthSq<RetType, Type>(v)));
	}

	template<typename RetType, typename Type>
	RetType length(const vec3_t<Type>& v) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return static_cast<RetType>(std::sqrt(lengthSq<RetType, Type>(v)));
	}

	template<typename RetType, typename Type>
	RetType length(const vec2_t<Type>& v) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return static_cast<RetType>(std::sqrt(lengthSq<RetType, Type>(v)));
	}

	template<typename RetType, typename Type>
	RetType length(const quat_base<Type>& q) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return static_cast<RetType>(std::sqrt(lengthSq<RetType, Type>(q)));
	}

	template<typename RetType, typename Type>
	RetType length(const complex_base<Type>& c) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return static_cast<RetType>(std::sqrt(lengthSq<RetType, Type>(c)));
	}
#pragma endregion

#pragma region(distanceSq)
	template<typename RetType, typename LeftType, typename RightType, int LeftDim, int RightDim>
	RetType distanceSq(const vec_base<LeftType, LeftDim>& lhs, const vec_base<RightType, RightDim>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return 0;
	}

	//vec nt distanceSq vec nt
	template<typename RetType, typename LeftType, typename RightType, int Dim>
	RetType distanceSq(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return lengthSq<RetType, RetType, RetType>(sub<RetType, LeftType, RightType, Dim>(lhs, rhs));
	}

	//vec 4t distanceSq vec 4t
	template<typename RetType, typename LeftType, typename RightType>
	RetType distanceSq(const vec4_t<LeftType>& lhs, const vec4_t<RightType>& rhs) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return lengthSq<RetType, RetType, RetType>(sub<RetType, LeftType, RightType>(lhs, rhs));
	}

	//vec 3t distanceSq vec 3t
	template<typename RetType, typename LeftType, typename RightType>
	RetType distanceSq(const vec3_t<LeftType>& lhs, const vec3_t<RightType>& rhs) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return lengthSq<RetType, RetType, RetType>(sub<RetType, LeftType, RightType>(lhs, rhs));
	}

	//vec 2t distanceSq vec 2t
	template<typename RetType, typename LeftType, typename RightType>
	RetType distanceSq(const vec2_t<LeftType>& lhs, const vec2_t<RightType>& rhs) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return lengthSq<RetType, RetType, RetType>(sub<RetType, LeftType, RightType>(lhs, rhs));
	}

	//quat distanceSq quat
	template<typename RetType, typename LeftType, typename RightType>
	RetType distanceSq(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return lengthSq<RetType, RetType, RetType>(sub<RetType, LeftType, RightType>(lhs, rhs));
	}

	//complex distanceSq complex
	template<typename RetType, typename LeftType, typename RightType>
	RetType distanceSq(const complex_base<LeftType>& lhs, const complex_base<RightType>& rhs) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return lengthSq<RetType, RetType, RetType>(sub<RetType, LeftType, RightType>(lhs, rhs));
	}
#pragma endregion

#pragma region(distance)
	template<typename RetType, typename LeftType, typename RightType, int LeftDim, int RightDim>
	RetType distance(const vec_base<LeftType, LeftDim>& lhs, const vec_base<RightType, RightDim>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return 0;
	}

	//vec nt distance vec nt
	template<typename RetType, typename LeftType, typename RightType, int Dim>
	RetType distance(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return static_cast<RetType>(std::sqrt(distanceSq<RetType, LeftType, RightType, Dim>(lhs, rhs)));
	}

	//vec 4t distance vec 4t
	template<typename RetType, typename LeftType, typename RightType>
	RetType distance(const vec4_t<LeftType>& lhs, const vec4_t<RightType>& rhs) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return static_cast<RetType>(std::sqrt(distanceSq<RetType, LeftType, RightType>(lhs, rhs)));
	}

	//vec 3t distance vec 3t
	template<typename RetType, typename LeftType, typename RightType>
	RetType distance(const vec3_t<LeftType>& lhs, const vec3_t<RightType>& rhs) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return static_cast<RetType>(std::sqrt(distanceSq<RetType, LeftType, RightType>(lhs, rhs)));
	}

	//vec 2t distance vec 2t
	template<typename RetType, typename LeftType, typename RightType>
	RetType distance(const vec2_t<LeftType>& lhs, const vec2_t<RightType>& rhs) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return static_cast<RetType>(std::sqrt(distanceSq<RetType, LeftType, RightType>(lhs, rhs)));
	}

	//quat distance quat
	template<typename RetType, typename LeftType, typename RightType>
	RetType distance(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return static_cast<RetType>(std::sqrt(distanceSq<RetType, LeftType, RightType>(lhs, rhs)));
	}

	//complex distance complex
	template<typename RetType, typename LeftType, typename RightType>
	RetType distance(const complex_base<LeftType>& lhs, const complex_base<RightType>& rhs) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		return static_cast<RetType>(std::sqrt(distanceSq<RetType, LeftType, RightType>(lhs, rhs)));
	}
#pragma endregion

#pragma region(min)
	template<typename RetType, typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol>
	auto min(const mat_base<LeftType, LeftRow, LeftCol>& value, const mat_base<RightType, RightRow, RightCol>& range) {
		Z_MATH_ASSERT_EQUIVALENT;
		return mat_base<LeftType, LeftRow, LeftCol>();
	}
	
	template<typename RetType, typename LeftType, typename RightType, int LeftDim, int RightDim>
	auto min(const vec_base<LeftType, LeftDim>& value, const vec_base<RightType, RightDim>& range) {
		Z_MATH_ASSERT_EQUIVALENT;
		return vec_base<LeftType, LeftDim>();
	}

	//n x m matrix
	template<typename RetType, typename LeftType, typename RightType, int Row, int Col>
	auto min(const mat_base<LeftType, Row, Col>& value, const mat_base<RightType, Row, Col>& range) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		using RetT = mat_base<RetType, Row, Col>;
		RetT ret;
		for (int i = 0; i != Row; ++i) {
			ret[i] = min<RetType, LeftType, RightType, Col>(value[i], range[i]);
		}
		return ret;
	}

	//n x m matrix min scalar
	template<typename RetType, typename LeftType, typename RightType, int Row, int Col>
	auto min(const mat_base<LeftType, Row, Col>& value, const RightType& range) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		using RetT = mat_base<RetType, Row, Col>;
		RetT ret;
		for (int i = 0; i != Row; ++i) {
			ret[i] = min<RetType, LeftType, RightType, Col>(value[i], range);
		}
		return ret;
	}

	//n x m matrix min scalar
	template<typename RetType, typename LeftType, typename RightType, int Row, int Col>
	auto min(const LeftType& value, const mat_base<RightType, Row, Col>& range) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		using RetT = mat_base<RetType, Row, Col>;
		RetT ret;
		for (int i = 0; i != Row; ++i) {
			ret[i] = min<RetType, LeftType, RightType, Col>(value, range[i]);
		}
		return ret;
	}

	//n dim vector 
	template<typename RetType, typename LeftType, typename RightType, int Dim>
	auto min(const vec_base<LeftType, Dim>& value, const vec_base<RightType, Dim>& range) {
		using RetT = vec_base<RetType, Dim>;
		RetT ret;
		for (int i = 0; i != Dim; ++i) {
			ret[i] = min<RetType, LeftType, RightType>(value[i], range[i]);
		}
		return ret;
	}

	//n dim vector min scalar 
	template<typename RetType, typename LeftType, typename RightType, int Dim>
	auto min(const vec_base<LeftType, Dim>& value, const RightType& range) {
		using RetT = vec_base<RetType, Dim>;
		RetT ret;
		for (int i = 0; i != Dim; ++i) {
			ret[i] = min<RetType, LeftType, RightType>(value[i], range);
		}
		return ret;
	}

	//4 x 4 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto min(const mat4x4_t<LeftType>& value, const mat4x4_t<RightType>& range) {
		using RetT = mat4x4_t<RetType>;
		return RetT{
			min<RetType, LeftType, RightType>(value[0], range[0]),
			min<RetType, LeftType, RightType>(value[1], range[1]),
			min<RetType, LeftType, RightType>(value[2], range[2]),
			min<RetType, LeftType, RightType>(value[3], range[3])
		};
	}

	//3 x 3 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto min(const mat3x3_t<LeftType>& value, const mat3x3_t<RightType>& range) {
		using RetT = mat3x3_t<RetType>;
		return RetT{
			min<RetType, LeftType, RightType>(value[0] , range[0]),
			min<RetType, LeftType, RightType>(value[1] , range[1]),
			min<RetType, LeftType, RightType>(value[2] , range[2])
		};
	}

	//2 x 2 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto min(const mat2x2_t<LeftType>& value, const mat2x2_t<RightType>& range) {
		using RetT = mat2x2_t<RetType>;
		return RetT{
			min<RetType, LeftType, RightType>(value[0] , range[0]),
			min<RetType, LeftType, RightType>(value[1] , range[1])
		};
	}

	//4 dim vector 
	template<typename RetType, typename LeftType, typename RightType>
	auto min(const vec4_t<LeftType>& value, const vec4_t<RightType>& range) {
		using RetT = vec4_t<RetType>;
		return RetT{
			min<RetType, LeftType, RightType>(value[0], range[0]),
			min<RetType, LeftType, RightType>(value[1], range[1]),
			min<RetType, LeftType, RightType>(value[2], range[2]),
			min<RetType, LeftType, RightType>(value[3], range[3])
		};
	}

	//3 dim vector 
	template<typename RetType, typename LeftType, typename RightType>
	auto min(const vec3_t<LeftType>& value, const vec3_t<RightType>& range) {
		using RetT = vec3_t<RetType>;
		return RetT{
			min<RetType, LeftType, RightType>(value[0], range[0]),
			min<RetType, LeftType, RightType>(value[1], range[1]),
			min<RetType, LeftType, RightType>(value[2], range[2])
		};
	}

	//2 dim vector 
	template<typename RetType, typename LeftType, typename RightType>
	auto min(const vec2_t<LeftType>& value, const vec2_t<RightType>& range) {
		using RetT = vec2_t<RetType>;
		return RetT{
			min<RetType, LeftType, RightType>(value[0], range[0]),
			min<RetType, LeftType, RightType>(value[1], range[1])
		};
	}

	//quat
	template<typename RetType, typename LeftType, typename RightType>
	auto min(const quat_base<LeftType>& value, const quat_base<RightType>& range) {
		using RetT = quat_base<RetType>;
		return RetT{
			min<RetType, LeftType, RightType>(value[0], range[0]),
			min<RetType, LeftType, RightType>(value[1], range[1]),
			min<RetType, LeftType, RightType>(value[2], range[2]),
			min<RetType, LeftType, RightType>(value[3], range[3])
		};
	}

	//complex
	template<typename RetType, typename LeftType, typename RightType>
	auto min(const complex_base<LeftType>& value, const complex_base<RightType>& range) {
		using RetT = complex_base<RetType>;
		return RetT{
			min<RetType, LeftType, RightType>(value[0], range[0]),
			min<RetType, LeftType, RightType>(value[1], range[1])
		};
	}

	//2 x 2 matrix min scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto min(const mat2x2_t<LeftType>& value, const RightType& range) {
		using RetT = mat2x2_t<RetType>;
		return RetT{
			min<RetType, LeftType, RightType>(value[0], range),
			min<RetType, LeftType, RightType>(value[1], range)
		};
	}

	//3 x 3 matrix min scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto min(const mat3x3_t<LeftType>& value, const RightType& range) {
		using RetT = mat3x3_t<RetType>;
		return RetT{
			min<RetType, LeftType, RightType>(value[0], range),
			min<RetType, LeftType, RightType>(value[1], range),
			min<RetType, LeftType, RightType>(value[2], range)
		};
	}

	//4 x 4 matrix min scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto min(const mat4x4_t<LeftType>& value, const RightType& range) {
		using RetT = mat4x4_t<RetType>;
		return RetT{
			min<RetType, LeftType, RightType>(value[0], range),
			min<RetType, LeftType, RightType>(value[1], range),
			min<RetType, LeftType, RightType>(value[2], range),
			min<RetType, LeftType, RightType>(value[3], range)
		};
	}

	//vec 2t min scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto min(const vec2_t<LeftType>& value, const RightType& range) {
		using RetT = vec2_t<RetType>;
		return RetT{
			min<RetType, LeftType, RightType>(value[0], range),
			min<RetType, LeftType, RightType>(value[1], range)
		};
	}

	//vec 3t min scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto min(const vec3_t<LeftType>& value, const RightType& range) {
		using RetT = vec3_t<RetType>;
		return RetT{
			min<RetType, LeftType, RightType>(value[0], range),
			min<RetType, LeftType, RightType>(value[1], range),
			min<RetType, LeftType, RightType>(value[2], range)
		};
	}

	//vec 4t min scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto min(const vec4_t<LeftType>& value, const RightType& range) {
		using RetT = vec4_t<RetType>;
		return RetT{
			min<RetType, LeftType, RightType>(value[0], range),
			min<RetType, LeftType, RightType>(value[1], range),
			min<RetType, LeftType, RightType>(value[2], range),
			min<RetType, LeftType, RightType>(value[3], range)
		};
	}

	//quat min scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto min(const quat_base<LeftType>& value, const RightType& range) {
		using RetT = quat_base<RetType>;
		return RetT{
			min<RetType, LeftType, RightType>(value[0], range),
			min<RetType, LeftType, RightType>(value[1], range),
			min<RetType, LeftType, RightType>(value[2], range),
			min<RetType, LeftType, RightType>(value[3], range)
		};
	}

	//complex min scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto min(const complex_base<LeftType>& value, const RightType& range) {
		using RetT = complex_base<RetType>;
		return RetT{
			min<RetType, LeftType, RightType>(value[0], range),
			min<RetType, LeftType, RightType>(value[1], range)
		};
	}
#pragma endregion

#pragma region(max)
	template<typename RetType, typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol>
	auto max(const mat_base<LeftType, LeftRow, LeftCol>& value, const mat_base<RightType, RightRow, RightCol>& range) {
		Z_MATH_ASSERT_EQUIVALENT;
		return mat_base<LeftType, LeftRow, LeftCol>();
	}

	template<typename RetType, typename LeftType, typename RightType, int LeftDim, int RightDim>
	auto max(const vec_base<LeftType, LeftDim>& value, const vec_base<RightType, RightDim>& range) {
		Z_MATH_ASSERT_EQUIVALENT;
		return vec_base<LeftType, LeftDim>();
	}

	//n x m matrix
	template<typename RetType, typename LeftType, typename RightType, int Row, int Col>
	auto max(const mat_base<LeftType, Row, Col>& value, const mat_base<RightType, Row, Col>& range) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		using RetT = mat_base<RetType, Row, Col>;
		RetT ret;
		for (int i = 0; i != Row; ++i) {
			ret[i] = max<RetType, LeftType, RightType, Col>(value[i], range[i]);
		}
		return ret;
	}

	//n x m matrix max scalar
	template<typename RetType, typename LeftType, typename RightType, int Row, int Col>
	auto max(const mat_base<LeftType, Row, Col>& value, const RightType& range) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		using RetT = mat_base<RetType, Row, Col>;
		RetT ret;
		for (int i = 0; i != Row; ++i) {
			ret[i] = max<RetType, LeftType, RightType, Col>(value[i], range);
		}
		return ret;
	}

	//n x m matrix max scalar
	template<typename RetType, typename LeftType, typename RightType, int Row, int Col>
	auto max(const LeftType& value, const mat_base<RightType, Row, Col>& range) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		using RetT = mat_base<RetType, Row, Col>;
		RetT ret;
		for (int i = 0; i != Row; ++i) {
			ret[i] = max<RetType, LeftType, RightType, Col>(value, range[i]);
		}
		return ret;
	}

	//n dim vector 
	template<typename RetType, typename LeftType, typename RightType, int Dim>
	auto max(const vec_base<LeftType, Dim>& value, const vec_base<RightType, Dim>& range) {
		using RetT = vec_base<RetType, Dim>;
		RetT ret;
		for (int i = 0; i != Dim; ++i) {
			ret[i] = max<RetType, LeftType, RightType>(value[i], range[i]);
		}
		return ret;
	}

	//n dim vector max scalar 
	template<typename RetType, typename LeftType, typename RightType, int Dim>
	auto max(const vec_base<LeftType, Dim>& value, const RightType& range) {
		using RetT = vec_base<RetType, Dim>;
		RetT ret;
		for (int i = 0; i != Dim; ++i) {
			ret[i] = max<RetType, LeftType, RightType>(value[i], range);
		}
		return ret;
	}

	//4 x 4 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto max(const mat4x4_t<LeftType>& value, const mat4x4_t<RightType>& range) {
		using RetT = mat4x4_t<RetType>;
		return RetT{
			max<RetType, LeftType, RightType>(value[0], range[0]),
			max<RetType, LeftType, RightType>(value[1], range[1]),
			max<RetType, LeftType, RightType>(value[2], range[2]),
			max<RetType, LeftType, RightType>(value[3], range[3])
		};
	}

	//3 x 3 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto max(const mat3x3_t<LeftType>& value, const mat3x3_t<RightType>& range) {
		using RetT = mat3x3_t<RetType>;
		return RetT{
			max<RetType, LeftType, RightType>(value[0] , range[0]),
			max<RetType, LeftType, RightType>(value[1] , range[1]),
			max<RetType, LeftType, RightType>(value[2] , range[2])
		};
	}

	//2 x 2 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto max(const mat2x2_t<LeftType>& value, const mat2x2_t<RightType>& range) {
		using RetT = mat2x2_t<RetType>;
		return RetT{
			max<RetType, LeftType, RightType>(value[0] , range[0]),
			max<RetType, LeftType, RightType>(value[1] , range[1])
		};
	}

	//4 dim vector 
	template<typename RetType, typename LeftType, typename RightType>
	auto max(const vec4_t<LeftType>& value, const vec4_t<RightType>& range) {
		using RetT = vec4_t<RetType>;
		return RetT{
			max<RetType, LeftType, RightType>(value[0], range[0]),
			max<RetType, LeftType, RightType>(value[1], range[1]),
			max<RetType, LeftType, RightType>(value[2], range[2]),
			max<RetType, LeftType, RightType>(value[3], range[3])
		};
	}

	//3 dim vector 
	template<typename RetType, typename LeftType, typename RightType>
	auto max(const vec3_t<LeftType>& value, const vec3_t<RightType>& range) {
		using RetT = vec3_t<RetType>;
		return RetT{
			max<RetType, LeftType, RightType>(value[0], range[0]),
			max<RetType, LeftType, RightType>(value[1], range[1]),
			max<RetType, LeftType, RightType>(value[2], range[2])
		};
	}

	//2 dim vector 
	template<typename RetType, typename LeftType, typename RightType>
	auto max(const vec2_t<LeftType>& value, const vec2_t<RightType>& range) {
		using RetT = vec2_t<RetType>;
		return RetT{
			max<RetType, LeftType, RightType>(value[0], range[0]),
			max<RetType, LeftType, RightType>(value[1], range[1])
		};
	}

	//quat
	template<typename RetType, typename LeftType, typename RightType>
	auto max(const quat_base<LeftType>& value, const quat_base<RightType>& range) {
		using RetT = quat_base<RetType>;
		return RetT{
			max<RetType, LeftType, RightType>(value[0], range[0]),
			max<RetType, LeftType, RightType>(value[1], range[1]),
			max<RetType, LeftType, RightType>(value[2], range[2]),
			max<RetType, LeftType, RightType>(value[3], range[3])
		};
	}

	//complex
	template<typename RetType, typename LeftType, typename RightType>
	auto max(const complex_base<LeftType>& value, const complex_base<RightType>& range) {
		using RetT = complex_base<RetType>;
		return RetT{
			max<RetType, LeftType, RightType>(value[0], range[0]),
			max<RetType, LeftType, RightType>(value[1], range[1])
		};
	}

	//2 x 2 matrix max scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto max(const mat2x2_t<LeftType>& value, const RightType& range) {
		using RetT = mat2x2_t<RetType>;
		return RetT{
			max<RetType, LeftType, RightType>(value[0], range),
			max<RetType, LeftType, RightType>(value[1], range)
		};
	}

	//3 x 3 matrix max scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto max(const mat3x3_t<LeftType>& value, const RightType& range) {
		using RetT = mat3x3_t<RetType>;
		return RetT{
			max<RetType, LeftType, RightType>(value[0], range),
			max<RetType, LeftType, RightType>(value[1], range),
			max<RetType, LeftType, RightType>(value[2], range)
		};
	}

	//4 x 4 matrix max scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto max(const mat4x4_t<LeftType>& value, const RightType& range) {
		using RetT = mat4x4_t<RetType>;
		return RetT{
			max<RetType, LeftType, RightType>(value[0], range),
			max<RetType, LeftType, RightType>(value[1], range),
			max<RetType, LeftType, RightType>(value[2], range),
			max<RetType, LeftType, RightType>(value[3], range)
		};
	}

	//vec 2t max scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto max(const vec2_t<LeftType>& value, const RightType& range) {
		using RetT = vec2_t<RetType>;
		return RetT{
			max<RetType, LeftType, RightType>(value[0], range),
			max<RetType, LeftType, RightType>(value[1], range)
		};
	}

	//vec 3t max scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto max(const vec3_t<LeftType>& value, const RightType& range) {
		using RetT = vec3_t<RetType>;
		return RetT{
			max<RetType, LeftType, RightType>(value[0], range),
			max<RetType, LeftType, RightType>(value[1], range),
			max<RetType, LeftType, RightType>(value[2], range)
		};
	}

	//vec 4t max scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto max(const vec4_t<LeftType>& value, const RightType& range) {
		using RetT = vec4_t<RetType>;
		return RetT{
			max<RetType, LeftType, RightType>(value[0], range),
			max<RetType, LeftType, RightType>(value[1], range),
			max<RetType, LeftType, RightType>(value[2], range),
			max<RetType, LeftType, RightType>(value[3], range)
		};
	}

	//quat max scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto max(const quat_base<LeftType>& value, const RightType& range) {
		using RetT = quat_base<RetType>;
		return RetT{
			max<RetType, LeftType, RightType>(value[0], range),
			max<RetType, LeftType, RightType>(value[1], range),
			max<RetType, LeftType, RightType>(value[2], range),
			max<RetType, LeftType, RightType>(value[3], range)
		};
	}

	//complex max scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto max(const complex_base<LeftType>& value, const RightType& range) {
		using RetT = complex_base<RetType>;
		return RetT{
			max<RetType, LeftType, RightType>(value[0], range),
			max<RetType, LeftType, RightType>(value[1], range)
		};
	}
#pragma endregion

#pragma region(clamp)
	template<typename RetType, typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol>
	auto clamp(const mat_base<LeftType, LeftRow, LeftCol>& value, const mat_base<RightType, RightRow, RightCol>& begin, const mat_base<RightType, RightRow, RightCol>& end) {
		Z_MATH_ASSERT_EQUIVALENT;
		return mat_base<LeftType, LeftRow, LeftCol>();
	}

	template<typename RetType, typename LeftType, typename RightType, int LeftDim, int RightDim>
	auto clamp(const vec_base<LeftType, LeftDim>& value, const vec_base<RightType, RightDim>& begin, const vec_base<RightType, RightDim>& end) {
		Z_MATH_ASSERT_EQUIVALENT;
		return vec_base<LeftType, LeftDim>();
	}

	//n x m matrix
	template<typename RetType, typename LeftType, typename RightType, int Row, int Col>
	auto clamp(const mat_base<LeftType, Row, Col>& value, const mat_base<RightType, Row, Col>& begin, const mat_base<RightType, Row, Col>& end) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		using RetT = mat_base<RetType, Row, Col>;
		RetT ret;
		for (int i = 0; i != Row; ++i) {
			ret[i] = clamp<RetType, LeftType, RightType, Col>(value[i], begin[i], end[i]);
		}
		return ret;
	}

	//n x m matrix clamp scalar
	template<typename RetType, typename LeftType, typename RightType, int Row, int Col>
	auto clamp(const mat_base<LeftType, Row, Col>& value, const RightType& begin, const RightType& end) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		using RetT = mat_base<RetType, Row, Col>;
		RetT ret;
		for (int i = 0; i != Row; ++i) {
			ret[i] = clamp<RetType, LeftType, RightType, Col>(value[i], begin, end);
		}
		return ret;
	}

	//n x m matrix clamp scalar
	template<typename RetType, typename LeftType, typename RightType, int Row, int Col>
	auto clamp(const LeftType& value, const mat_base<RightType, Row, Col>& begin, const mat_base<RightType, Row, Col>& end) {
		static_assert(is_arithmetic<RetType>::value, __FUNCTION__ " : Return type must arithmetic type.");
		using RetT = mat_base<RetType, Row, Col>;
		RetT ret;
		for (int i = 0; i != Row; ++i) {
			ret[i] = clamp<RetType, LeftType, RightType, Col>(value, begin[i], end[i]);
		}
		return ret;
	}

	//n dim vector 
	template<typename RetType, typename LeftType, typename RightType, int Dim>
	auto clamp(const vec_base<LeftType, Dim>& value, const vec_base<RightType, Dim>& begin, const vec_base<RightType, Dim>& end) {
		using RetT = vec_base<RetType, Dim>;
		RetT ret;
		for (int i = 0; i != Dim; ++i) {
			ret[i] = clamp<RetType, LeftType, RightType>(value[i], begin[i], end[i]);
		}
		return ret;
	}

	//n dim vector clamp scalar 
	template<typename RetType, typename LeftType, typename RightType, int Dim>
	auto clamp(const vec_base<LeftType, Dim>& value, const RightType& begin, const RightType& end) {
		using RetT = vec_base<RetType, Dim>;
		RetT ret;
		for (int i = 0; i != Dim; ++i) {
			ret[i] = clamp<RetType, LeftType, RightType>(value[i], begin, end);
		}
		return ret;
	}

	//4 x 4 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto clamp(const mat4x4_t<LeftType>& value, const mat4x4_t<RightType>& begin, const mat4x4_t<RightType>& end) {
		using RetT = mat4x4_t<RetType>;
		return RetT{
			clamp<RetType, LeftType, RightType>(value[0], begin[0], end[0]),
			clamp<RetType, LeftType, RightType>(value[1], begin[1], end[1]),
			clamp<RetType, LeftType, RightType>(value[2], begin[2], end[2]),
			clamp<RetType, LeftType, RightType>(value[3], begin[3], end[3])
		};
	}

	//3 x 3 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto clamp(const mat3x3_t<LeftType>& value, const mat3x3_t<RightType>& begin, const mat3x3_t<RightType>& end) {
		using RetT = mat3x3_t<RetType>;
		return RetT{
			clamp<RetType, LeftType, RightType>(value[0] , begin[0], end[0]),
			clamp<RetType, LeftType, RightType>(value[1] , begin[1], end[1]),
			clamp<RetType, LeftType, RightType>(value[2] , begin[2], end[2])
		};
	}

	//2 x 2 matrix
	template<typename RetType, typename LeftType, typename RightType>
	auto clamp(const mat2x2_t<LeftType>& value, const mat2x2_t<RightType>& begin, const mat2x2_t<RightType>& end) {
		using RetT = mat2x2_t<RetType>;
		return RetT{
			clamp<RetType, LeftType, RightType>(value[0] , begin[0], end[0]),
			clamp<RetType, LeftType, RightType>(value[1] , begin[1], end[1])
		};
	}

	//4 dim vector 
	template<typename RetType, typename LeftType, typename RightType>
	auto clamp(const vec4_t<LeftType>& value, const vec4_t<RightType>& begin, const vec4_t<RightType>& end) {
		using RetT = vec4_t<RetType>;
		return RetT{
			clamp<RetType, LeftType, RightType>(value[0], begin[0], end[0]),
			clamp<RetType, LeftType, RightType>(value[1], begin[1], end[1]),
			clamp<RetType, LeftType, RightType>(value[2], begin[2], end[2]),
			clamp<RetType, LeftType, RightType>(value[3], begin[3], end[3])
		};
	}

	//3 dim vector 
	template<typename RetType, typename LeftType, typename RightType>
	auto clamp(const vec3_t<LeftType>& value, const vec3_t<RightType>& begin, const vec3_t<RightType>& end) {
		using RetT = vec3_t<RetType>;
		return RetT{
			clamp<RetType, LeftType, RightType>(value[0], begin[0], end[0]),
			clamp<RetType, LeftType, RightType>(value[1], begin[1], end[1]),
			clamp<RetType, LeftType, RightType>(value[2], begin[2], end[2])
		};
	}

	//2 dim vector 
	template<typename RetType, typename LeftType, typename RightType>
	auto clamp(const vec2_t<LeftType>& value, const vec2_t<RightType>& begin, const vec2_t<RightType>& end) {
		using RetT = vec2_t<RetType>;
		return RetT{
			clamp<RetType, LeftType, RightType>(value[0], begin[0], end[0]),
			clamp<RetType, LeftType, RightType>(value[1], begin[1], end[1])
		};
	}

	//quat
	template<typename RetType, typename LeftType, typename RightType>
	auto clamp(const quat_base<LeftType>& value, const quat_base<RightType>& begin, const quat_base<RightType>& end) {
		using RetT = quat_base<RetType>;
		return RetT{
			clamp<RetType, LeftType, RightType>(value[0], begin[0], end[0]),
			clamp<RetType, LeftType, RightType>(value[1], begin[1], end[1]),
			clamp<RetType, LeftType, RightType>(value[2], begin[2], end[2]),
			clamp<RetType, LeftType, RightType>(value[3], begin[3], end[3])
		};
	}

	//complex
	template<typename RetType, typename LeftType, typename RightType>
	auto clamp(const complex_base<LeftType>& value, const complex_base<RightType>& begin, const complex_base<RightType>& end) {
		using RetT = complex_base<RetType>;
		return RetT{
			clamp<RetType, LeftType, RightType>(value[0], begin[0], end[0]),
			clamp<RetType, LeftType, RightType>(value[1], begin[1], end[1])
		};
	}

	//2 x 2 matrix clamp scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto clamp(const mat2x2_t<LeftType>& value, const RightType& begin, const RightType& end) {
		using RetT = mat2x2_t<RetType>;
		return RetT{
			clamp<RetType, LeftType, RightType>(value[0], begin, end),
			clamp<RetType, LeftType, RightType>(value[1], begin, end)
		};
	}

	//3 x 3 matrix clamp scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto clamp(const mat3x3_t<LeftType>& value, const RightType& begin, const RightType& end) {
		using RetT = mat3x3_t<RetType>;
		return RetT{
			clamp<RetType, LeftType, RightType>(value[0], begin, end),
			clamp<RetType, LeftType, RightType>(value[1], begin, end),
			clamp<RetType, LeftType, RightType>(value[2], begin, end)
		};
	}

	//4 x 4 matrix clamp scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto clamp(const mat4x4_t<LeftType>& value, const RightType& begin, const RightType& end) {
		using RetT = mat4x4_t<RetType>;
		return RetT{
			clamp<RetType, LeftType, RightType>(value[0], begin, end),
			clamp<RetType, LeftType, RightType>(value[1], begin, end),
			clamp<RetType, LeftType, RightType>(value[2], begin, end),
			clamp<RetType, LeftType, RightType>(value[3], begin, end)
		};
	}

	//vec 2t clamp scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto clamp(const vec2_t<LeftType>& value, const RightType& begin, const RightType& end) {
		using RetT = vec2_t<RetType>;
		return RetT{
			clamp<RetType, LeftType, RightType>(value[0], begin, end),
			clamp<RetType, LeftType, RightType>(value[1], begin, end)
		};
	}

	//vec 3t clamp scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto clamp(const vec3_t<LeftType>& value, const RightType& begin, const RightType& end) {
		using RetT = vec3_t<RetType>;
		return RetT{
			clamp<RetType, LeftType, RightType>(value[0], begin, end),
			clamp<RetType, LeftType, RightType>(value[1], begin, end),
			clamp<RetType, LeftType, RightType>(value[2], begin, end)
		};
	}

	//vec 4t clamp scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto clamp(const vec4_t<LeftType>& value, const RightType& begin, const RightType& end) {
		using RetT = vec4_t<RetType>;
		return RetT{
			clamp<RetType, LeftType, RightType>(value[0], begin, end),
			clamp<RetType, LeftType, RightType>(value[1], begin, end),
			clamp<RetType, LeftType, RightType>(value[2], begin, end),
			clamp<RetType, LeftType, RightType>(value[3], begin, end)
		};
	}

	//quat clamp scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto clamp(const quat_base<LeftType>& value, const RightType& begin, const RightType& end) {
		using RetT = quat_base<RetType>;
		return RetT{
			clamp<RetType, LeftType, RightType>(value[0], begin, end),
			clamp<RetType, LeftType, RightType>(value[1], begin, end),
			clamp<RetType, LeftType, RightType>(value[2], begin, end),
			clamp<RetType, LeftType, RightType>(value[3], begin, end)
		};
	}

	//complex clamp scalar
	template<typename RetType, typename LeftType, typename RightType>
	auto clamp(const complex_base<LeftType>& value, const RightType& begin, const RightType& end) {
		using RetT = complex_base<RetType>;
		return RetT{
			clamp<RetType, LeftType, RightType>(value[0], begin, end),
			clamp<RetType, LeftType, RightType>(value[1], begin, end)
		};
	}
#pragma endregion

#pragma region(saturate)
	//n x m matrix
	template<typename RetType, typename Type, int Row, int Col>
	auto saturate(const mat_base<Type, Row, Col>& v) {
		using RetT = mat_base<RetType, Row, Col>;
		RetT ret;
		for (int i = 0; i != Row; ++i) {
			ret[i] = saturate<RetType, Type, Col>(v[i]);
		}
		return ret;
	}

	//n dim vector 
	template<typename RetType, typename Type, int Dim>
	auto saturate(const vec_base<Type, Dim>& v) {
		using RetT = vec_base<RetType, Dim>;
		RetT ret;
		for (int d = 0; d != Dim; ++d) {
			ret[d] = saturate<RetType, Type>(v[d]);
		}
		return ret;
	}

	//4 x 4 matrix
	template<typename RetType, typename Type>
	auto saturate(const mat4x4_t<Type>& v) {
		using RetT = mat4x4_t<RetType>;
		return RetT{
			saturate<RetType, Type>(v[0]),
			saturate<RetType, Type>(v[1]),
			saturate<RetType, Type>(v[2]),
			saturate<RetType, Type>(v[3])
		};
	}

	//3 x 3 matrix
	template<typename RetType, typename Type>
	auto saturate(const mat3x3_t<Type>& v) {
		using RetT = mat3x3_t<RetType>;
		return RetT{
			saturate<RetType, Type>(v[0]),
			saturate<RetType, Type>(v[1]),
			saturate<RetType, Type>(v[2])
		};
	}

	//2 x 2 matrix
	template<typename RetType, typename Type>
	auto saturate(const mat2x2_t<Type>& v) {
		using RetT = mat4x4_t<RetType>;
		return RetT{
			saturate<RetType, Type>(v[0]),
			saturate<RetType, Type>(v[1])
		};
	}

	//4 dim vector 
	template<typename RetType, typename Type>
	auto saturate(const vec4_t<Type>& v) {
		using RetT = vec4_t<RetType>;
		return RetT{
			saturate<RetType, Type>(v[0]),
			saturate<RetType, Type>(v[1]),
			saturate<RetType, Type>(v[2]),
			saturate<RetType, Type>(v[3])
		};
	}

	//3 dim vector 
	template<typename RetType, typename Type>
	auto saturate(const vec3_t<Type>& v) {
		using RetT = vec3_t<RetType>;
		return RetT{
			saturate<RetType, Type>(v[0]),
			saturate<RetType, Type>(v[1]),
			saturate<RetType, Type>(v[2])
		};
	}

	//2 dim vector 
	template<typename RetType, typename Type>
	auto saturate(const vec2_t<Type>& v) {
		using RetT = vec2_t<RetType>;
		return RetT{
			saturate<RetType, Type>(v[0]),
			saturate<RetType, Type>(v[1])
		};
	}

	//quat
	template<typename RetType, typename Type>
	auto saturate(const quat_base<Type>& v) {
		using RetT = quat_base<RetType>;
		return RetT{
			saturate<RetType, Type>(v[0]),
			saturate<RetType, Type>(v[1]),
			saturate<RetType, Type>(v[2]),
			saturate<RetType, Type>(v[3])
		};
	}

	//complex
	template<typename RetType, typename Type>
	auto saturate(const complex_base<Type>& v) {
		using RetT = complex_base<RetType>;
		return RetT{
			saturate<RetType, Type>(v[0]),
			saturate<RetType, Type>(v[1])
		};
	}
#pragma endregion

#pragma region(lerp)
	template<typename RetType, typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol, typename FloatT>
	auto lerp(const mat_base<LeftType, LeftRow, LeftCol>& value, const mat_base<RightType, RightRow, RightCol>& begin, const mat_base<RightType, RightRow, RightCol>& end) {
		Z_MATH_ASSERT_EQUIVALENT;
		return mat_base<LeftType, LeftRow, LeftCol>();
	}

	template<typename RetType, typename LeftType, typename RightType, int LeftDim, int RightDim>
	auto lerp(const vec_base<LeftType, LeftDim>& value, const vec_base<RightType, RightDim>& begin, const vec_base<RightType, RightDim>& end) {
		Z_MATH_ASSERT_EQUIVALENT;
		return vec_base<LeftType, LeftDim>();
	}

	//n x m matrix
	template<typename RetType, typename LeftType, typename RightType, int Row, int Col, typename FloatT>
	auto lerp(const mat_base<LeftType, Row, Col>& v0, const mat_base<RightType, Row, Col>& v1, const FloatT& t) {
		using RetT = mat_base<RetType, Row, Col>;
		RetT ret;
		for (int i = 0; i != Row; ++i) {
			ret[i] = lerp<RetType, LeftType, RightType, Col, FloatT>(v0[i], v1[i], t);
		}
		return ret;
	}

	//vec nt
	template<typename RetType, typename LeftType, typename RightType, int Dim, typename FloatT>
	auto lerp(const vec_base<LeftType, Dim>& v0, const vec_base<RightType, Dim>& v1, const FloatT& t) {
		using RetT = vec_base<RetType, Dim>;
		RetT ret;
		for (int d = 0; d != Dim; ++d) {
			ret[d] = lerp<RetType, LeftType, RightType, FloatT>(v0[d], v1[d], t);
		}
		return ret;
	}

	//4 x 4 matrix
	template<typename RetType, typename LeftType, typename RightType, typename FloatT>
	auto lerp(const mat4x4_t<LeftType>& v0, const mat4x4_t<RightType>& v1, const FloatT& t) {
		using RetT = mat4x4_t<RetType>;
		return RetT{
			lerp<RetType, LeftType, RightType, FloatT>(v0[0],v1[0],t),
			lerp<RetType, LeftType, RightType, FloatT>(v0[1],v1[1],t),
			lerp<RetType, LeftType, RightType, FloatT>(v0[2],v1[2],t),
			lerp<RetType, LeftType, RightType, FloatT>(v0[3],v1[3],t)
		};
	}

	//3 x 3 matrix
	template<typename RetType, typename LeftType, typename RightType, typename FloatT>
	auto lerp(const mat3x3_t<LeftType>& v0, const mat3x3_t<RightType>& v1, const FloatT& t) {
		using RetT = mat3x3_t<RetType>;
		return RetT{
			lerp<RetType, LeftType, RightType, FloatT>(v0[0],v1[0],t),
			lerp<RetType, LeftType, RightType, FloatT>(v0[1],v1[1],t),
			lerp<RetType, LeftType, RightType, FloatT>(v0[2],v1[2],t)
		};
	}

	//2 x 2 matrix
	template<typename RetType, typename LeftType, typename RightType, typename FloatT>
	auto lerp(const mat2x2_t<LeftType>& v0, const mat2x2_t<RightType>& v1, const FloatT& t) {
		using RetT = mat2x2_t<RetType>;
		return RetT{
			lerp<RetType, LeftType, RightType, FloatT>(v0[0],v1[0],t),
			lerp<RetType, LeftType, RightType, FloatT>(v0[1],v1[1],t)
		};
	}

	//vec 4t
	template<typename RetType, typename LeftType, typename RightType, typename FloatT>
	auto lerp(const vec4_t<LeftType>& v0, const vec4_t<RightType>& v1, const FloatT& t) {
		using RetT = vec4_t<RetType>;
		return RetT{
			lerp<RetType,LeftType,RightType,FloatT>(v0[0],v1[0],t),
			lerp<RetType,LeftType,RightType,FloatT>(v0[1],v1[1],t),
			lerp<RetType,LeftType,RightType,FloatT>(v0[2],v1[2],t),
			lerp<RetType,LeftType,RightType,FloatT>(v0[3],v1[3],t)
		};
	}

	//vec 3t
	template<typename RetType, typename LeftType, typename RightType, typename FloatT>
	auto lerp(const vec3_t<LeftType>& v0, const vec3_t<RightType>& v1, const FloatT& t) {
		using RetT = vec3_t<RetType>;
		return RetT{
			lerp<RetType,LeftType,RightType,FloatT>(v0[0],v1[0],t),
			lerp<RetType,LeftType,RightType,FloatT>(v0[1],v1[1],t),
			lerp<RetType,LeftType,RightType,FloatT>(v0[2],v1[2],t)
		};
	}

	//vec 2t
	template<typename RetType, typename LeftType, typename RightType, typename FloatT>
	auto lerp(const vec2_t<LeftType>& v0, const vec2_t<RightType>& v1, const FloatT& t) {
		using RetT = vec2_t<RetType>;
		return RetT{
			lerp<RetType,LeftType,RightType,FloatT>(v0[0],v1[0],t),
			lerp<RetType,LeftType,RightType,FloatT>(v0[1],v1[1],t)
		};
	}

	//quat
	template<typename RetType, typename LeftType, typename RightType, typename FloatT>
	auto lerp(const quat_base<LeftType>& v0, const quat_base<RightType>& v1, const FloatT& t) {
		using RetT = quat_base<RetType>;
		return RetT{
			lerp<RetType,LeftType,RightType,FloatT>(v0[0],v1[0],t),
			lerp<RetType,LeftType,RightType,FloatT>(v0[1],v1[1],t),
			lerp<RetType,LeftType,RightType,FloatT>(v0[2],v1[2],t),
			lerp<RetType,LeftType,RightType,FloatT>(v0[3],v1[3],t)
		};
	}

	//complex
	template<typename RetType, typename LeftType, typename RightType, typename FloatT>
	auto lerp(const complex_base<LeftType>& v0, const complex_base<RightType>& v1, const FloatT& t) {
		using RetT = complex_base<RetType>;
		return RetT{
			lerp<RetType,LeftType,RightType,FloatT>(v0[0],v1[0],t),
			lerp<RetType,LeftType,RightType,FloatT>(v0[1],v1[1],t)
		};
	}
#pragma endregion

#pragma region(hermite1)

	//n x m matrix
	template<typename RetType, typename Type, int Row, int Col, typename FloatT>
	auto hermite(
		const mat_base<Type, Row, Col>& v0, 
		const mat_base<Type, Row, Col>& tan0, 
		const mat_base<Type, Row, Col>& v1, 
		const mat_base<Type, Row, Col>& tan1, const FloatT& t ) {
		using RetT = mat_base<RetType, Row, Col>;

		using PromotionType = promotion_t<RetType, Type, FloatT>;
		PromotionType exp[4];
		Details::hermite_support(static_cast<PromotionType>(t), exp);

		RetT ret;

		for (int r = 0; r != Row; ++r) {
			for (int c = 0; c != Col; ++c) {
				ret[r][c] = Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[r][c], tan0[r][c], v1[r][c], tan1[r][c]);
			}
		}

		return ret;
	}

	// vec nt
	template<typename RetType, typename Type, int Dim, typename FloatT>
	auto hermite(
		const vec_base<Type, Dim>& v0,
		const vec_base<Type, Dim>& tan0,
		const vec_base<Type, Dim>& v1,
		const vec_base<Type, Dim>& tan1, const FloatT& t) {

		using RetT = vec_base<RetType, Dim>;

		using PromotionType = promotion_t<RetType, Type, FloatT>;
		PromotionType exp[4];
		Details::hermite_support(static_cast<PromotionType>(t), exp);

		RetT ret;

		for (int d = 0; d != Dim; ++d) {
			ret[d] = Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[d], tan0[d], v1[d], tan1[d]);
		}

		return ret;
	}

	// vec4t
	template<typename RetType, typename Type, typename FloatT>
	auto hermite(
		const vec4_t<Type>& v0,
		const vec4_t<Type>& tan0,
		const vec4_t<Type>& v1,
		const vec4_t<Type>& tan1, const FloatT& t) {

		using RetT = vec4_t<RetType>;

		using PromotionType = promotion_t<RetType, Type, FloatT>;
		PromotionType exp[4];
		Details::hermite_support(static_cast<PromotionType>(t), exp);

		return RetT{
			Details::hermite_calculate<RetType,PromotionType,Type>(exp, v0[0], tan0[0], v1[0], tan1[0]),
			Details::hermite_calculate<RetType,PromotionType,Type>(exp, v0[1], tan0[1], v1[1], tan1[1]),
			Details::hermite_calculate<RetType,PromotionType,Type>(exp, v0[2], tan0[2], v1[2], tan1[2]),
			Details::hermite_calculate<RetType,PromotionType,Type>(exp, v0[3], tan0[3], v1[3], tan1[3])
		};
	}

	// vec3t
	template<typename RetType, typename Type, typename FloatT>
	auto hermite(
		const vec3_t<Type>& v0,
		const vec3_t<Type>& tan0,
		const vec3_t<Type>& v1,
		const vec3_t<Type>& tan1, const FloatT& t) {

		using RetT = vec3_t<RetType>;

		using PromotionType = promotion_t<RetType, Type, FloatT>;
		PromotionType exp[4];
		Details::hermite_support(static_cast<PromotionType>(t), exp);

		return RetT{
			Details::hermite_calculate<RetType,PromotionType,Type>(exp, v0[0], tan0[0], v1[0], tan1[0]),
			Details::hermite_calculate<RetType,PromotionType,Type>(exp, v0[1], tan0[1], v1[1], tan1[1]),
			Details::hermite_calculate<RetType,PromotionType,Type>(exp, v0[2], tan0[2], v1[2], tan1[2])
		};
	}

	// vec2t
	template<typename RetType, typename Type, typename FloatT>
	auto hermite(
		const vec2_t<Type>& v0,
		const vec2_t<Type>& tan0,
		const vec2_t<Type>& v1,
		const vec2_t<Type>& tan1, const FloatT& t) {

		using RetT = vec2_t<RetType>;

		using PromotionType = promotion_t<RetType, Type, FloatT>;
		PromotionType exp[4];
		Details::hermite_support(static_cast<PromotionType>(t), exp);

		return RetT{
			Details::hermite_calculate<RetType,PromotionType,Type>(exp, v0[0], tan0[0], v1[0], tan1[0]),
			Details::hermite_calculate<RetType,PromotionType,Type>(exp, v0[1], tan0[1], v1[1], tan1[1])
		};
	}

	//4 x 4 matrix
	template<typename RetType, typename Type, typename FloatT>
	auto hermite(
		const mat4x4_t<Type>& v0,
		const mat4x4_t<Type>& tan0,
		const mat4x4_t<Type>& v1,
		const mat4x4_t<Type>& tan1, const FloatT& t) {

		using RetT = mat4x4_t<Type>;
		using VecT = typename RetT::element_type;

		using PromotionType = promotion_t<RetType, Type, FloatT>;
		PromotionType exp[4];
		Details::hermite_support(static_cast<PromotionType>(t), exp);

		return RetT{
			VecT{
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[0][0], tan0[0][0], v1[0][0], tan1[0][0]) ,
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[0][1], tan0[0][1], v1[0][1], tan1[0][1]) ,
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[0][2], tan0[0][2], v1[0][2], tan1[0][2]) ,
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[0][3], tan0[0][3], v1[0][3], tan1[0][3])},

			VecT{
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[1][0], tan0[1][0], v1[1][0], tan1[1][0]) ,
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[1][1], tan0[1][1], v1[1][1], tan1[1][1]) ,
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[1][2], tan0[1][2], v1[1][2], tan1[1][2]) ,
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[1][3], tan0[1][3], v1[1][3], tan1[1][3])},

			VecT{
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[2][0], tan0[2][0], v1[2][0], tan1[2][0]) ,
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[2][1], tan0[2][1], v1[2][1], tan1[2][1]) ,
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[2][2], tan0[2][2], v1[2][2], tan1[2][2]) ,
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[2][3], tan0[2][3], v1[2][3], tan1[2][3])},

			VecT{
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[3][0], tan0[3][0], v1[3][0], tan1[3][0]) ,
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[3][1], tan0[3][1], v1[3][1], tan1[3][1]) ,
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[3][2], tan0[3][2], v1[3][2], tan1[3][2]) ,
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[3][3], tan0[3][3], v1[3][3], tan1[3][3])}
		};
	}

	//3 x 3 matrix
	template<typename RetType, typename Type, typename FloatT>
	auto hermite(
		const mat3x3_t<Type>& v0,
		const mat3x3_t<Type>& tan0,
		const mat3x3_t<Type>& v1,
		const mat3x3_t<Type>& tan1, const FloatT& t) {

		using RetT = mat3x3_t<Type>;
		using VecT = typename RetT::element_type;

		using PromotionType = promotion_t<RetType, Type, FloatT>;
		PromotionType exp[4];
		Details::hermite_support(static_cast<PromotionType>(t), exp);

		return RetT{
			VecT{
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[0][0], tan0[0][0], v1[0][0], tan1[0][0]) ,
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[0][1], tan0[0][1], v1[0][1], tan1[0][1]) ,
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[0][2], tan0[0][2], v1[0][2], tan1[0][2])},

			VecT{
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[1][0], tan0[1][0], v1[1][0], tan1[1][0]) ,
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[1][1], tan0[1][1], v1[1][1], tan1[1][1]) ,
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[1][2], tan0[1][2], v1[1][2], tan1[1][2])},

			VecT{
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[2][0], tan0[2][0], v1[2][0], tan1[2][0]) ,
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[2][1], tan0[2][1], v1[2][1], tan1[2][1]) ,
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[2][2], tan0[2][2], v1[2][2], tan1[2][2])}
		};
	}

	//2 x 2 matrix
	template<typename RetType, typename Type, typename FloatT>
	auto hermite(
		const mat2x2_t<Type>& v0,
		const mat2x2_t<Type>& tan0,
		const mat2x2_t<Type>& v1,
		const mat2x2_t<Type>& tan1, const FloatT& t) {

		using RetT = mat2x2_t<Type>;
		using VecT = typename RetT::element_type;

		using PromotionType = promotion_t<RetType, Type, FloatT>;
		PromotionType exp[4];
		Details::hermite_support(static_cast<PromotionType>(t), exp);

		return RetT{
			VecT{
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[0][0], tan0[0][0], v1[0][0], tan1[0][0]) ,
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[0][1], tan0[0][1], v1[0][1], tan1[0][1])},

			VecT{
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[1][0], tan0[1][0], v1[1][0], tan1[1][0]) ,
			Details::hermite_calculate<RetType, PromotionType, Type>(exp, v0[1][1], tan0[1][1], v1[1][1], tan1[1][1])}
		};
	}

	// quat
	template<typename RetType, typename Type, typename FloatT>
	auto hermite(
		const quat_base<Type>& v0,
		const quat_base<Type>& tan0,
		const quat_base<Type>& v1,
		const quat_base<Type>& tan1, const FloatT& t) {

		using RetT = quat_base<RetType>;

		using PromotionType = promotion_t<RetType, Type, FloatT>;
		PromotionType exp[4];
		Details::hermite_support(static_cast<PromotionType>(t), exp);

		return RetT{
			Details::hermite_calculate<RetType,PromotionType,Type>(exp, v0[0], tan0[0], v1[0], tan1[0]),
			Details::hermite_calculate<RetType,PromotionType,Type>(exp, v0[1], tan0[1], v1[1], tan1[1]),
			Details::hermite_calculate<RetType,PromotionType,Type>(exp, v0[2], tan0[2], v1[2], tan1[2]),
			Details::hermite_calculate<RetType,PromotionType,Type>(exp, v0[3], tan0[3], v1[3], tan1[3])
		};
	}

	// complex
	template<typename RetType, typename Type, typename FloatT>
	auto hermite(
		const complex_base<Type>& v0,
		const complex_base<Type>& tan0,
		const complex_base<Type>& v1,
		const complex_base<Type>& tan1, const FloatT& t) {

		using RetT = complex_base<RetType>;

		using PromotionType = promotion_t<RetType, Type, FloatT>;
		PromotionType exp[4];
		Details::hermite_support(static_cast<PromotionType>(t), exp);

		return RetT{
			Details::hermite_calculate<RetType,PromotionType,Type>(exp, v0[0], tan0[0], v1[0], tan1[0]),
			Details::hermite_calculate<RetType,PromotionType,Type>(exp, v0[1], tan0[1], v1[1], tan1[1])
		};
	}
#pragma endregion

#pragma region(hermite2)

	//n x m matrix
	template<typename RetType, typename Type, int Row, int Col, typename FloatT>
	auto hermite(
		const mat_base<Type, Row, Col>& v0,
		const mat_base<Type, Row, Col>& v1,
		const mat_base<Type, Row, Col>& v2, const FloatT& t) {
		return hermite<RetType, Type, Row, Col, FloatT>(
			v0, sub<Type>(v1, v0) , v1, sub<Type>(v2, v1) , t
		);
	}

	// vec nt
	template<typename RetType, typename Type, int Dim, typename FloatT>
	auto hermite(
		const vec_base<Type, Dim>& v0,
		const vec_base<Type, Dim>& v1,
		const vec_base<Type, Dim>& v2, const FloatT& t) {
		return hermite<RetType, Type, Dim, FloatT>(
			v0, sub<Type>(v1, v0), v1, sub<Type>(v2, v1), t
		);
	}

	// vec4t
	template<typename RetType, typename Type, typename FloatT>
	auto hermite(
		const vec4_t<Type>& v0,
		const vec4_t<Type>& v1,
		const vec4_t<Type>& v2, const FloatT& t) {
		return hermite<RetType, Type, FloatT>(v0, sub<Type>(v1, v0), v1, sub<Type>(v2, v1), t);
	}

	// vec3t
	template<typename RetType, typename Type, typename FloatT>
	auto hermite(
		const vec3_t<Type>& v0,
		const vec3_t<Type>& v1,
		const vec3_t<Type>& v2,const FloatT& t) {
		return hermite<RetType, Type, FloatT>(v0, sub<Type>(v1, v0), v1, sub<Type>(v2, v1), t);
	}

	// vec2t
	template<typename RetType, typename Type, typename FloatT>
	auto hermite(
		const vec2_t<Type>& v0,
		const vec2_t<Type>& v1,
		const vec2_t<Type>& v2, const FloatT& t) {
		return hermite<RetType, Type, FloatT>(v0, sub<Type>(v1, v0), v1, sub<Type>(v2, v1), t);
	}

	//4 x 4 matrix
	template<typename RetType, typename Type, typename FloatT>
	auto hermite(
		const mat4x4_t<Type>& v0,
		const mat4x4_t<Type>& v1,
		const mat4x4_t<Type>& v2, const FloatT& t) {
		return hermite<RetType, Type, FloatT>(v0, sub<Type>(v1, v0), v1, sub<Type>(v2, v1), t);
	}

	//3 x 3 matrix
	template<typename RetType, typename Type, typename FloatT>
	auto hermite(
		const mat3x3_t<Type>& v0,
		const mat3x3_t<Type>& v1,
		const mat3x3_t<Type>& v2, const FloatT& t) {
		return hermite<RetType, Type, FloatT>(v0, sub<Type>(v1, v0), v1, sub<Type>(v2, v1), t);
	}

	//2 x 2 matrix
	template<typename RetType, typename Type, typename FloatT>
	auto hermite(
		const mat2x2_t<Type>& v0,
		const mat2x2_t<Type>& v1,
		const mat2x2_t<Type>& v2, const FloatT& t) {
		return hermite<RetType, Type, FloatT>(v0, sub<Type>(v1, v0), v1, sub<Type>(v2, v1), t);
	}

	// quat
	template<typename RetType, typename Type, typename FloatT>
	auto hermite(
		const quat_base<Type>& v0,
		const quat_base<Type>& v1,
		const quat_base<Type>& v2, const FloatT& t) {
		return hermite<RetType, Type, FloatT>(v0, sub<Type>(v1, v0), v1, sub<Type>(v2, v1), t);
	}

	// complex
	template<typename RetType, typename Type, typename FloatT>
	auto hermite(
		const complex_base<Type>& v0,
		const complex_base<Type>& v1,
		const complex_base<Type>& v2, const FloatT& t) {
		return hermite<RetType, Type, FloatT>(v0, sub<Type>(v1, v0), v1, sub<Type>(v2, v1), t);
	}
#pragma endregion

#pragma region(catmullRom)
	//n x m matrix
	template<typename RetType, typename Type, int Row, int Col, typename FloatT>
	auto catmullRom(
		const mat_base<Type, Row, Col>& v0,
		const mat_base<Type, Row, Col>& v1,
		const mat_base<Type, Row, Col>& v2,
		const mat_base<Type, Row, Col>& v3, const FloatT& t) {
		using RetT = mat_base<RetType, Row, Col>;

		using PromotionType = promotion_t<RetType, Type, FloatT>;
		PromotionType exp[4];
		Details::catmullRom_support(static_cast<PromotionType>(t), exp);

		RetT ret;

		for (int r = 0; r != Row; ++r) {
			for (int c = 0; c != Col; ++c) {
				ret[r][c] = Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[r][c], v1[r][c], v2[r][c], v3[r][c]);
			}
		}

		return ret;
	}

	// vec nt
	template<typename RetType, typename Type, int Dim, typename FloatT>
	auto catmullRom(
		const vec_base<Type, Dim>& v0,
		const vec_base<Type, Dim>& v1,
		const vec_base<Type, Dim>& v2,
		const vec_base<Type, Dim>& v3, const FloatT& t) {

		using RetT = vec_base<RetType, Dim>;

		using PromotionType = promotion_t<RetType, Type, FloatT>;
		PromotionType exp[4];
		Details::catmullRom_support(static_cast<PromotionType>(t), exp);

		RetT ret;

		for (int d = 0; d != Dim; ++d) {
			ret[d] = Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[d], v1[d], v2[d], v3[d]);
		}

		return ret;
	}

	// vec4t
	template<typename RetType, typename Type, typename FloatT>
	auto catmullRom(
		const vec4_t<Type>& v0,
		const vec4_t<Type>& v1,
		const vec4_t<Type>& v2,
		const vec4_t<Type>& v3, const FloatT& t) {

		using RetT = vec4_t<RetType>;

		using PromotionType = promotion_t<RetType, Type, FloatT>;
		PromotionType exp[4];
		Details::catmullRom_support(static_cast<PromotionType>(t), exp);

		return RetT{
			Details::catmullRom_calculate<RetType,PromotionType,Type>(exp, v0[0], v1[0], v2[0], v3[0]),
			Details::catmullRom_calculate<RetType,PromotionType,Type>(exp, v0[1], v1[1], v2[1], v3[1]),
			Details::catmullRom_calculate<RetType,PromotionType,Type>(exp, v0[2], v1[2], v2[2], v3[2]),
			Details::catmullRom_calculate<RetType,PromotionType,Type>(exp, v0[3], v1[3], v2[3], v3[3])
		};
	}

	// vec3t
	template<typename RetType, typename Type, typename FloatT>
	auto catmullRom(
		const vec3_t<Type>& v0,
		const vec3_t<Type>& v1,
		const vec3_t<Type>& v2,
		const vec3_t<Type>& v3, const FloatT& t) {

		using RetT = vec3_t<RetType>;

		using PromotionType = promotion_t<RetType, Type, FloatT>;
		PromotionType exp[4];
		Details::catmullRom_support(static_cast<PromotionType>(t), exp);

		return RetT{
			Details::catmullRom_calculate<RetType,PromotionType,Type>(exp, v0[0], v1[0], v2[0], v3[0]),
			Details::catmullRom_calculate<RetType,PromotionType,Type>(exp, v0[1], v1[1], v2[1], v3[1]),
			Details::catmullRom_calculate<RetType,PromotionType,Type>(exp, v0[2], v1[2], v2[2], v3[2])
		};
	}

	// vec2t
	template<typename RetType, typename Type, typename FloatT>
	auto catmullRom(
		const vec2_t<Type>& v0,
		const vec2_t<Type>& v1,
		const vec2_t<Type>& v2,
		const vec2_t<Type>& v3, const FloatT& t) {

		using RetT = vec2_t<RetType>;

		using PromotionType = promotion_t<RetType, Type, FloatT>;
		PromotionType exp[4];
		Details::catmullRom_support(static_cast<PromotionType>(t), exp);

		return RetT{
			Details::catmullRom_calculate<RetType,PromotionType,Type>(exp, v0[0], v1[0], v2[0], v3[0]),
			Details::catmullRom_calculate<RetType,PromotionType,Type>(exp, v0[1], v1[1], v2[1], v3[1])
		};
	}

	//4 x 4 matrix
	template<typename RetType, typename Type, typename FloatT>
	auto catmullRom(
		const mat4x4_t<Type>& v0,
		const mat4x4_t<Type>& v1,
		const mat4x4_t<Type>& v2,
		const mat4x4_t<Type>& v3, const FloatT& t) {

		using RetT = mat4x4_t<Type>;
		using VecT = typename RetT::element_type;

		using PromotionType = promotion_t<RetType, Type, FloatT>;
		PromotionType exp[4];
		Details::catmullRom_support(static_cast<PromotionType>(t), exp);

		return RetT{
			VecT{
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[0][0], v1[0][0], v2[0][0], v3[0][0]) ,
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[0][1], v1[0][1], v2[0][1], v3[0][1]) ,
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[0][2], v1[0][2], v2[0][2], v3[0][2]) ,
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[0][3], v1[0][3], v2[0][3], v3[0][3])},

			VecT{
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[1][0], v1[1][0], v2[1][0], v3[1][0]) ,
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[1][1], v1[1][1], v2[1][1], v3[1][1]) ,
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[1][2], v1[1][2], v2[1][2], v3[1][2]) ,
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[1][3], v1[1][3], v2[1][3], v3[1][3])},

			VecT{
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[2][0], v1[2][0], v2[2][0], v3[2][0]) ,
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[2][1], v1[2][1], v2[2][1], v3[2][1]) ,
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[2][2], v1[2][2], v2[2][2], v3[2][2]) ,
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[2][3], v1[2][3], v2[2][3], v3[2][3])},

			VecT{
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[3][0], v1[3][0], v2[3][0], v3[3][0]) ,
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[3][1], v1[3][1], v2[3][1], v3[3][1]) ,
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[3][2], v1[3][2], v2[3][2], v3[3][2]) ,
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[3][3], v1[3][3], v2[3][3], v3[3][3])}
		};
	}

	//3 x 3 matrix
	template<typename RetType, typename Type, typename FloatT>
	auto catmullRom(
		const mat3x3_t<Type>& v0,
		const mat3x3_t<Type>& v1,
		const mat3x3_t<Type>& v2,
		const mat3x3_t<Type>& v3, const FloatT& t) {

		using RetT = mat3x3_t<Type>;
		using VecT = typename RetT::element_type;

		using PromotionType = promotion_t<RetType, Type, FloatT>;
		PromotionType exp[4];
		Details::catmullRom_support(static_cast<PromotionType>(t), exp);

		return RetT{
			VecT{
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[0][0], v1[0][0], v2[0][0], v3[0][0]) ,
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[0][1], v1[0][1], v2[0][1], v3[0][1]) ,
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[0][2], v1[0][2], v2[0][2], v3[0][2])},

			VecT{
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[1][0], v1[1][0], v2[1][0], v3[1][0]) ,
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[1][1], v1[1][1], v2[1][1], v3[1][1]) ,
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[1][2], v1[1][2], v2[1][2], v3[1][2])},

			VecT{
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[2][0], v1[2][0], v2[2][0], v3[2][0]) ,
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[2][1], v1[2][1], v2[2][1], v3[2][1]) ,
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[2][2], v1[2][2], v2[2][2], v3[2][2])}
		};
	}

	//2 x 2 matrix
	template<typename RetType, typename Type, typename FloatT>
	auto catmullRom(
		const mat2x2_t<Type>& v0,
		const mat2x2_t<Type>& v1,
		const mat2x2_t<Type>& v2,
		const mat2x2_t<Type>& v3, const FloatT& t) {

		using RetT = mat2x2_t<Type>;
		using VecT = typename RetT::element_type;

		using PromotionType = promotion_t<RetType, Type, FloatT>;
		PromotionType exp[4];
		Details::catmullRom_support(static_cast<PromotionType>(t), exp);

		return RetT{
			VecT{
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[0][0], v1[0][0], v2[0][0], v3[0][0]) ,
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[0][1], v1[0][1], v2[0][1], v3[0][1])},

			VecT{
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[1][0], v1[1][0], v2[1][0], v3[1][0]) ,
			Details::catmullRom_calculate<RetType, PromotionType, Type>(exp, v0[1][1], v1[1][1], v2[1][1], v3[1][1])}
		};
	}

	// quat
	template<typename RetType, typename Type, typename FloatT>
	auto catmullRom(
		const quat_base<Type>& v0,
		const quat_base<Type>& v1,
		const quat_base<Type>& v2,
		const quat_base<Type>& v3, const FloatT& t) {

		using RetT = quat_base<RetType>;

		using PromotionType = promotion_t<RetType, Type, FloatT>;
		PromotionType exp[4];
		Details::catmullRom_support(static_cast<PromotionType>(t), exp);

		return RetT{
			Details::catmullRom_calculate<RetType,PromotionType,Type>(exp, v0[0], v1[0], v2[0], v3[0]),
			Details::catmullRom_calculate<RetType,PromotionType,Type>(exp, v0[1], v1[1], v2[1], v3[1]),
			Details::catmullRom_calculate<RetType,PromotionType,Type>(exp, v0[2], v1[2], v2[2], v3[2]),
			Details::catmullRom_calculate<RetType,PromotionType,Type>(exp, v0[3], v1[3], v2[3], v3[3])
		};
	}

	// complex
	template<typename RetType, typename Type, typename FloatT>
	auto catmullRom(
		const complex_base<Type>& v0,
		const complex_base<Type>& v1,
		const complex_base<Type>& v2,
		const complex_base<Type>& v3, const FloatT& t) {

		using RetT = complex_base<RetType>;

		using PromotionType = promotion_t<RetType, Type, FloatT>;
		PromotionType exp[4];
		Details::catmullRom_support(static_cast<PromotionType>(t), exp);

		return RetT{
			Details::catmullRom_calculate<RetType,PromotionType,Type>(exp, v0[0], v1[0], v2[0], v3[0]),
			Details::catmullRom_calculate<RetType,PromotionType,Type>(exp, v0[1], v1[1], v2[1], v3[1])
		};
	}

#pragma endregion

#pragma region(smoothStep)
	template<typename RetType, typename LeftType, typename RightType, int LeftRow, int LeftCol, int RightRow, int RightCol, typename FloatT>
	auto smoothStep(const mat_base<LeftType, LeftRow, LeftCol>& value, const mat_base<RightType, RightRow, RightCol>& begin, const mat_base<RightType, RightRow, RightCol>& end) {
		Z_MATH_ASSERT_EQUIVALENT;
		return mat_base<LeftType, LeftRow, LeftCol>();
	}

	template<typename RetType, typename LeftType, typename RightType, int LeftDim, int RightDim>
	auto smoothStep(const vec_base<LeftType, LeftDim>& value, const vec_base<RightType, RightDim>& begin, const vec_base<RightType, RightDim>& end) {
		Z_MATH_ASSERT_EQUIVALENT;
		return vec_base<LeftType, LeftDim>();
	}

	//n x m matrix
	template<typename RetType, typename LeftType, typename RightType, int Row, int Col, typename FloatT>
	auto smoothStep(const mat_base<LeftType, Row, Col>& v0, const mat_base<RightType, Row, Col>& v1, const FloatT& t) {
		using RetT = mat_base<RetType, Row, Col>;
		RetT ret;
		for (int i = 0; i != Row; ++i) {
			ret[i] = smoothStep<RetType, LeftType, RightType, Col, FloatT>(v0[i], v1[i], t);
		}
		return ret;
	}

	//vec nt
	template<typename RetType, typename LeftType, typename RightType, int Dim, typename FloatT>
	auto smoothStep(const vec_base<LeftType, Dim>& v0, const vec_base<RightType, Dim>& v1, const FloatT& t) {
		using RetT = vec_base<RetType, Dim>;
		RetT ret;
		for (int d = 0; d != Dim; ++d) {
			ret[d] = smoothStep<RetType, LeftType, RightType, FloatT>(v0[d], v1[d], t);
		}
		return ret;
	}

	//4 x 4 matrix
	template<typename RetType, typename LeftType, typename RightType, typename FloatT>
	auto smoothStep(const mat4x4_t<LeftType>& v0, const mat4x4_t<RightType>& v1, const FloatT& t) {
		using RetT = mat4x4_t<RetType>;
		return RetT{
			smoothStep<RetType, LeftType, RightType, FloatT>(v0[0],v1[0],t),
			smoothStep<RetType, LeftType, RightType, FloatT>(v0[1],v1[1],t),
			smoothStep<RetType, LeftType, RightType, FloatT>(v0[2],v1[2],t),
			smoothStep<RetType, LeftType, RightType, FloatT>(v0[3],v1[3],t)
		};
	}

	//3 x 3 matrix
	template<typename RetType, typename LeftType, typename RightType, typename FloatT>
	auto smoothStep(const mat3x3_t<LeftType>& v0, const mat3x3_t<RightType>& v1, const FloatT& t) {
		using RetT = mat3x3_t<RetType>;
		return RetT{
			smoothStep<RetType, LeftType, RightType, FloatT>(v0[0],v1[0],t),
			smoothStep<RetType, LeftType, RightType, FloatT>(v0[1],v1[1],t),
			smoothStep<RetType, LeftType, RightType, FloatT>(v0[2],v1[2],t)
		};
	}

	//2 x 2 matrix
	template<typename RetType, typename LeftType, typename RightType, typename FloatT>
	auto smoothStep(const mat2x2_t<LeftType>& v0, const mat2x2_t<RightType>& v1, const FloatT& t) {
		using RetT = mat2x2_t<RetType>;
		return RetT{
			smoothStep<RetType, LeftType, RightType, FloatT>(v0[0],v1[0],t),
			smoothStep<RetType, LeftType, RightType, FloatT>(v0[1],v1[1],t)
		};
	}

	//vec 4t
	template<typename RetType, typename LeftType, typename RightType, typename FloatT>
	auto smoothStep(const vec4_t<LeftType>& v0, const vec4_t<RightType>& v1, const FloatT& t) {
		using RetT = vec4_t<RetType>;
		return RetT{
			smoothStep<RetType,LeftType,RightType,FloatT>(v0[0],v1[0],t),
			smoothStep<RetType,LeftType,RightType,FloatT>(v0[1],v1[1],t),
			smoothStep<RetType,LeftType,RightType,FloatT>(v0[2],v1[2],t),
			smoothStep<RetType,LeftType,RightType,FloatT>(v0[3],v1[3],t)
		};
	}

	//vec 3t
	template<typename RetType, typename LeftType, typename RightType, typename FloatT>
	auto smoothStep(const vec3_t<LeftType>& v0, const vec3_t<RightType>& v1, const FloatT& t) {
		using RetT = vec3_t<RetType>;
		return RetT{
			smoothStep<RetType,LeftType,RightType,FloatT>(v0[0],v1[0],t),
			smoothStep<RetType,LeftType,RightType,FloatT>(v0[1],v1[1],t),
			smoothStep<RetType,LeftType,RightType,FloatT>(v0[2],v1[2],t)
		};
	}

	//vec 2t
	template<typename RetType, typename LeftType, typename RightType, typename FloatT>
	auto smoothStep(const vec2_t<LeftType>& v0, const vec2_t<RightType>& v1, const FloatT& t) {
		using RetT = vec2_t<RetType>;
		return RetT{
			smoothStep<RetType,LeftType,RightType,FloatT>(v0[0],v1[0],t),
			smoothStep<RetType,LeftType,RightType,FloatT>(v0[1],v1[1],t)
		};
	}

	//quat
	template<typename RetType, typename LeftType, typename RightType, typename FloatT>
	auto smoothStep(const quat_base<LeftType>& v0, const quat_base<RightType>& v1, const FloatT& t) {
		using RetT = quat_base<RetType>;
		return RetT{
			smoothStep<RetType,LeftType,RightType,FloatT>(v0[0],v1[0],t),
			smoothStep<RetType,LeftType,RightType,FloatT>(v0[1],v1[1],t),
			smoothStep<RetType,LeftType,RightType,FloatT>(v0[2],v1[2],t),
			smoothStep<RetType,LeftType,RightType,FloatT>(v0[3],v1[3],t)
		};
	}

	//complex
	template<typename RetType, typename LeftType, typename RightType, typename FloatT>
	auto smoothStep(const complex_base<LeftType>& v0, const complex_base<RightType>& v1, const FloatT& t) {
		using RetT = complex_base<RetType>;
		return RetT{
			smoothStep<RetType,LeftType,RightType,FloatT>(v0[0],v1[0],t),
			smoothStep<RetType,LeftType,RightType,FloatT>(v0[1],v1[1],t)
		};
	}
#pragma endregion

#pragma region(make_normalize)
	template<typename RetType, typename Type, int Dim>
	auto make_normalize(const vec_base<Type, Dim>& v) {
		using RetT = vec_base<RetType, Dim>;
		RetType invLen = div(static_cast<RetType>(1), length<RetType>(v));
		RetT ret;
		for (int d = 0; d != Dim; ++d) {
			ret[d] = invLen * v[d];
		}
		return ret;
	}

	template<typename RetType, typename Type>
	auto make_normalize(const vec4_t<Type>& v) {
		using RetT = vec4_t<RetType>;
		RetType invLen = div<RetType>(static_cast<RetType>(1), length<RetType>(v));
		return RetT{ invLen * v[0], invLen * v[1], invLen * v[2], invLen * v[3]	};
	}

	template<typename RetType, typename Type>
	auto make_normalize(const vec3_t<Type>& v) {
		using RetT = vec3_t<RetType>;
		RetType invLen = div<RetType>(static_cast<RetType>(1), length<RetType>(v));
		return RetT{ invLen * v[0], invLen * v[1], invLen * v[2]};
	}

	template<typename RetType, typename Type>
	auto make_normalize(const vec2_t<Type>& v) {
		using RetT = vec2_t<RetType>;
		RetType invLen = div<RetType>(static_cast<RetType>(1), length<RetType>(v));
		return RetT{ invLen * v[0], invLen * v[1] };
	}

	template<typename RetType, typename Type>
	auto make_normalize(const quat_base<Type>& v) {
		using RetT = quat_base<RetType>;
		RetType invLen = div<RetType>(static_cast<RetType>(1), length<RetType>(v));
		return RetT{ invLen * v[0], invLen * v[1], invLen * v[2], invLen * v[3] };
	}

	template<typename RetType, typename Type>
	auto make_normalize(const complex_base<Type>& v) {
		using RetT = complex_base<RetType>;
		RetType invLen = div<RetType>(static_cast<RetType>(1), length<RetType>(v));
		return RetT{ invLen * v[0], invLen * v[1] };
	}
#pragma endregion
	
}//namespace Math
}//namespace Zee