#pragma once
#include "ZSpecializationTypes.h"

namespace Zee {
namespace Math {
	//vec == vec error handling
	template<typename LeftType, typename RightType, int LDim, int RDim>
	bool operator==(const vec_base<LeftType, LDim>& lhs, const vec_base<RightType, RDim>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	//vec != vec error handling
	template<typename LeftType, typename RightType, int LDim, int RDim>
	bool operator!=(const vec_base<LeftType, LDim>& lhs, const vec_base<RightType, RDim>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	//mat == mat error handling
	template<typename LeftType, typename RightType, int LRow, int LCol, int RRow, int RCol>
	bool operator==(const mat_base<LeftType, LRow, LCol>& lhs, const mat_base<RightType, RRow, RCol>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	//mat != mat error handling
	template<typename LeftType, typename RightType, int LRow, int LCol, int RRow, int RCol>
	bool operator!=(const mat_base<LeftType, LRow, LCol>& lhs, const mat_base<RightType, RRow, RCol>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return false;
	}

	//vec + vec error handling
	template<typename LeftType, typename RightType, int LDim, int RDim>
	vec_base<LeftType, LDim>
		operator+(const vec_base<LeftType, LDim>& lhs, const vec_base<RightType, RDim>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return vec_base<LeftType, LDim>();
	}

	//vec - vec error handling
	template<typename LeftType, typename RightType, int LDim, int RDim>
	vec_base<LeftType, LDim>
		operator-(const vec_base<LeftType, LDim>& lhs, const vec_base<RightType, RDim>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return vec_base<LeftType, LDim>();
	}

	//vec / vec error handling
	template<typename LeftType, typename RightType, int LDim, int RDim>
	vec_base<LeftType, LDim>
		operator/(const vec_base<LeftType, LDim>& lhs, const vec_base<RightType, RDim>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return vec_base<LeftType, LDim>();
	}

	//vec * vec error handling
	template<typename LeftType, typename RightType, int LDim, int RDim>
	vec_base<LeftType, LDim>
		operator*(const vec_base<LeftType, LDim>& lhs, const vec_base<RightType, RDim>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return vec_base<LeftType, LDim>();
	}

	//mat + mat error handling
	template<typename LeftType, typename RightType, int LRow, int LCol, int RRow, int RCol>
	mat_base<LeftType, LRow, LCol>
		operator+(const mat_base<LeftType, LRow, LCol>& lhs, const mat_base<LeftType, RRow, RCol>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return mat_base<LeftType, LRow, LCol>();
	}

	//mat - mat error handling
	template<typename LeftType, typename RightType, int LRow, int LCol, int RRow, int RCol>
	mat_base<LeftType, LRow, LCol>
		operator-(const mat_base<LeftType, LRow, LCol>& lhs, const mat_base<LeftType, RRow, RCol>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return mat_base<LeftType, LRow, LCol>();
	}

	//mat / mat error handling
	template<typename LeftType, typename RightType, int LRow, int LCol, int RRow, int RCol>
	mat_base<LeftType, LRow, LCol>
		operator/(const mat_base<LeftType, LRow, LCol>& lhs, const mat_base<LeftType, RRow, RCol>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return mat_base<LeftType, LRow, LCol>();
	}

	//vec *= mat error handling
	template<typename LeftType, typename RightType, int Dim, int Row, int Col>
	vec_base<LeftType, Dim>&
		operator*=(vec_base<LeftType, Dim>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		Z_MATH_ASSERT(0, __FUNCTION__ " : vector's dimension must same matrix row and col.");
		return lhs;
	}

	//mat *= vec error handling
	template<typename LeftType, typename RightType, int LRow, int LCol, int Dim>
	mat_base<LeftType, LRow, LCol>&
		operator*=(mat_base<LeftType, LRow, LCol>& lhs, const vec_base<RightType, Dim>& rhs) {
		Z_MATH_ASSERT(0, __FUNCTION__ " : invalid type. vec * mat operation must return vec.");
		return lhs;
	}

	//mat *= quat error handling
	template<typename LeftType, typename RightType, int Row, int Col>
	mat_base<LeftType, Row, Col>&
		operator*=(mat_base<LeftType, Row, Col>& lhs, const quat_base<RightType>& rhs) {
		Z_MATH_ASSERT(0, __FUNCTION__ " : invalid type.");
		return lhs;
	}

	//mat *= mat error handling
	template<typename LeftType, typename RightType, int LRow, int LCol, int RRow, int RCol>
	mat_base<LeftType, LRow, LCol>&
		operator*=(mat_base<LeftType, LRow, LCol>& lhs, const mat_base<RightType, RRow, RCol>& rhs) {
		Z_MATH_ASSERT_SQUARE;
		return lhs;
	}

	//mat /= vec error handling
	template<typename LeftType, typename RightType, int LRow, int LCol, int Dim>
	mat_base<LeftType, LRow, LCol>&
		operator/=(mat_base<LeftType, LRow, LCol>& lhs, const vec_base<RightType, Dim>& rhs) {
		Z_MATH_ASSERT(0, __FUNCTION__ " : invalid type. vec * mat operation must return vec.");
		return lhs;
	}

	//mat /= quat error handling
	template<typename LeftType, typename RightType, int Row, int Col>
	mat_base<LeftType, Row, Col>&
		operator/=(mat_base<LeftType, Row, Col>& lhs, const quat_base<RightType>& rhs) {
		Z_MATH_ASSERT(0, __FUNCTION__ " : invalid type.");
		return lhs;
	}

	//mat /= mat error handling
	template<typename LeftType, typename RightType, int LRow, int LCol, int RRow, int RCol>
	mat_base<LeftType, LRow, LCol>&
		operator/=(mat_base<LeftType, LRow, LCol>& lhs, const mat_base<RightType, RRow, RCol>& rhs) {
		Z_MATH_ASSERT_EQUIVALENT;
		return lhs;
	}

	//quat * vec error handling
	template<typename LeftType, typename RightType, int Dim>
	quat_base<LeftType>
		operator*(const quat_base<LeftType>& lhs, const vec_base<RightType, Dim>& rhs) {
		Z_MATH_ASSERT(0, __FUNCTION__ " : invalid type.");
		return quat_base<LeftType>();
	}

	//quat * mat error handling
	template<typename LeftType, typename RightType, int Row, int Col>
	quat_base<LeftType>
		operator*(const quat_base<LeftType>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		Z_MATH_ASSERT(0, __FUNCTION__ " : invalid type.");
		return quat_base<LeftType>();
	}

	//vec * quat error handling
	template<typename LeftType, typename RightType, int Dim>
	vec_base<LeftType, Dim>
		operator*(const vec_base<LeftType, Dim>& lhs, const quat_base<RightType>& rhs) {
		Z_MATH_ASSERT(0, __FUNCTION__ " : invalid type.");
		return vec_base<LeftType, Dim>();
	}

	//mat * quat error handling
	template<typename LeftType, typename RightType, int Row, int Col>
	quat_base<RightType>
		operator*(const mat_base<LeftType, Row, Col>& lhs, const quat_base<RightType>& rhs) {
		Z_MATH_ASSERT(0, __FUNCTION__ " : invalid type.");
		return quat_base<RightType>();
	}

	//vec == vec
	template<typename LeftType, typename RightType, int Dim>
	bool operator==(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		return Math::is_equal(lhs, rhs);
	}

	//vec != vec
	template<typename LeftType, typename RightType, int Dim>
	bool operator!=(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		return Math::is_not_equal(lhs, rhs);
	}

	//mat == mat
	template<typename LeftType, typename RightType, int Row, int Col>
	bool operator==(const mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		return Math::is_equal(lhs, rhs);
	}

	//mat != mat
	template<typename LeftType, typename RightType, int Row, int Col>
	bool operator!=(const mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		return Math::is_not_equal(lhs, rhs);
	}

	//quat == quat
	template<typename LeftType, typename RightType>
	bool operator==(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		return Math::is_equal(lhs, rhs);
	}

	//quat != quat
	template<typename LeftType, typename RightType, int Dim>
	bool operator!=(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		return Math::is_not_equal(lhs, rhs);
	}

	//vec + vec
	template<typename LeftType, typename RightType, int Dim>
	vec_base<promotion_t<LeftType, RightType>, Dim>
		operator+(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		return Math::add<promotion_t<LeftType, RightType>>(lhs, rhs);
	}

	//vec - vec
	template<typename LeftType, typename RightType, int Dim>
	vec_base<promotion_t<LeftType, RightType>, Dim>
		operator-(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		return Math::sub<promotion_t<LeftType, RightType>>(lhs, rhs);
	}

	//vec / vec
	template<typename LeftType, typename RightType, int Dim>
	vec_base<promotion_t<LeftType, RightType>, Dim>
		operator/(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		return Math::div<promotion_t<LeftType, RightType>>(lhs, rhs);
	}

	//vec * vec
	template<typename LeftType, typename RightType, int Dim>
	vec_base<promotion_t<LeftType, RightType>, Dim>
		operator*(const vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		return Math::mul<promotion_t<LeftType, RightType>>(lhs, rhs);
	}

	//mat + mat
	template<typename LeftType, typename RightType, int Row, int Col>
	mat_base<promotion_t<LeftType, RightType>, Row, Col>
		operator+(const mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		return Math::add<promotion_t<LeftType, RightType>>(lhs, rhs);
	}

	//mat - mat
	template<typename LeftType, typename RightType, int Row, int Col>
	mat_base<promotion_t<LeftType, RightType>, Row, Col>
		operator-(const mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		return Math::sub<promotion_t<LeftType, RightType>>(lhs, rhs);
	}

	//mat / mat
	template<typename LeftType, typename RightType, int Row, int Col>
	mat_base<promotion_t<LeftType, RightType>, Row, Col>
		operator/(const mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		return Math::div<promotion_t<LeftType, RightType>>(lhs, rhs);
	}

	//mat * mat
	template<typename LeftType, typename RightType, int LRow, int LCol, int RRow, int RCol>
	mat_base<promotion_t<LeftType, RightType>, LRow, RCol>
		operator*(const mat_base<LeftType, LRow, LCol>& lhs, const mat_base<RightType, RRow, RCol>& rhs) {
		return Math::mul<promotion_t<LeftType, RightType>>(lhs, rhs);
	}

	//quat + quat
	template<typename LeftType, typename RightType>
	quat_base<promotion_t<LeftType, RightType>>
		operator+(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		return Math::add<promotion_t<LeftType, RightType>>(lhs, rhs);
	}

	//quat - quat
	template<typename LeftType, typename RightType>
	quat_base<promotion_t<LeftType, RightType>>
		operator-(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		return Math::sub<promotion_t<LeftType, RightType>>(lhs, rhs);
	}

	//quat / quat
	template<typename LeftType, typename RightType>
	quat_base<promotion_t<LeftType, RightType>>
		operator/(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		return Math::div<promotion_t<LeftType, RightType>>(lhs, rhs);
	}

	//quat * quat
	template<typename LeftType, typename RightType>
	quat_base<promotion_t<LeftType, RightType>>
		operator*(const quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		return Math::mul<promotion_t<LeftType, RightType>>(lhs, rhs);
	}

	//vec * mat
	template<typename LeftType, typename RightType, int Row, int Col>
	vec_base<promotion_t<LeftType, RightType>, Col>
		operator*(const vec_base<LeftType, Col>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		return Math::mul<promotion_t<LeftType, RightType>>(lhs, rhs);
	}

	//mat * vec
	template<typename LeftType, typename RightType, int Row, int Col>
	vec_base<promotion_t<LeftType, RightType>, Row>
		operator*(const mat_base<LeftType, Row, Col>& lhs, const vec_base<RightType, Row>& rhs) {
		return Math::mul<promotion_t<LeftType, RightType>>(lhs, rhs);
	}

	//vec * scalar
	template<typename LeftType, typename RightType, int Dim>
	vec_base<promotion_t<LeftType, RightType>, Dim>
		operator*(const vec_base<LeftType, Dim>& lhs, const RightType& rhs) {
		return Math::mul<promotion_t<LeftType, RightType>>(lhs, rhs);
	}

	//vec / scalar
	template<typename LeftType, typename RightType, int Dim>
	vec_base<promotion_t<LeftType, RightType>, Dim>
		operator/(const vec_base<LeftType, Dim>& lhs, const RightType& rhs) {
		return Math::div<promotion_t<LeftType, RightType>>(lhs, rhs);
	}

	//scalar * vec
	template<typename LeftType, typename RightType, int Dim>
	vec_base<promotion_t<LeftType, RightType>, Dim>
		operator*(const LeftType& lhs, const vec_base<RightType, Dim>& rhs) {
		return Math::mul<promotion_t<LeftType, RightType>>(lhs, rhs);
	}

	//scalar / vec
	template<typename LeftType, typename RightType, int Dim>
	vec_base<promotion_t<LeftType, RightType>, Dim>
		operator/(const LeftType& lhs, const vec_base<RightType, Dim>& rhs) {
		return Math::div<promotion_t<LeftType, RightType>>(lhs, rhs);
	}


	//mat * scalar
	template<typename LeftType, typename RightType, int Row, int Col>
	mat_base<promotion_t<LeftType, RightType>, Row, Col>
		operator*(const mat_base<LeftType, Row, Col>& lhs, const RightType& rhs) {
		return Math::mul<promotion_t<LeftType, RightType>>(lhs, rhs);
	}

	//mat / scalar
	template<typename LeftType, typename RightType, int Row, int Col>
	mat_base<promotion_t<LeftType, RightType>, Row, Col>
		operator/(const mat_base<LeftType, Row, Col>& lhs, const RightType& rhs) {
		return Math::div<promotion_t<LeftType, RightType>>(lhs, rhs);
	}

	//scalar * mat
	template<typename LeftType, typename RightType, int Row, int Col>
	mat_base<promotion_t<LeftType, RightType>, Row, Col>
		operator*(const LeftType& lhs, const mat_base<RightType, Row, Col>& rhs) {
		return Math::mul<promotion_t<LeftType, RightType>>(lhs, rhs);
	}

	//scalar / mat
	template<typename LeftType, typename RightType, int Row, int Col>
	mat_base<promotion_t<LeftType, RightType>, Row, Col>
		operator/(const LeftType& lhs, const mat_base<RightType, Row, Col>& rhs) {
		return Math::div<promotion_t<LeftType, RightType>>(lhs, rhs);
	}

	//quat * scalar
	template<typename LeftType, typename RightType>
	quat_base<promotion_t<LeftType, RightType>>
		operator*(const quat_base<LeftType>& lhs, const RightType& rhs) {
		return Math::mul<promotion_t<LeftType, RightType>>(lhs, rhs);
	}

	//quat / scalar
	template<typename LeftType, typename RightType>
	quat_base<promotion_t<LeftType, RightType>>
		operator/(const quat_base<LeftType>& lhs, const RightType& rhs) {
		return Math::div<promotion_t<LeftType, RightType>>(lhs, rhs);
	}

	//scalar * quat
	template<typename LeftType, typename RightType>
	quat_base<promotion_t<LeftType, RightType>>
		operator*(const LeftType& lhs, const quat_base<RightType>& rhs) {
		return Math::mul<promotion_t<LeftType, RightType>>(lhs, rhs);
	}

	//scalar / quat
	template<typename LeftType, typename RightType>
	quat_base<promotion_t<LeftType, RightType>>
		operator/(const LeftType& lhs, const quat_base<RightType>& rhs) {
		return Math::div<promotion_t<LeftType, RightType>>(lhs, rhs);
	}

	//combination assign += -= *= /=
	template<typename LeftType, typename RightType, int Dim>
	vec_base<LeftType, Dim>&
		operator+=(vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		return lhs = Math::add<LeftType>(lhs, rhs);
	}

	template<typename LeftType, typename RightType, int Dim>
	vec_base<LeftType, Dim>&
		operator-=(vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		return lhs = Math::sub<LeftType>(lhs, rhs);
	}

	template<typename LeftType, typename RightType, int Dim>
	vec_base<LeftType, Dim>&
		operator/=(vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		return lhs = Math::div<LeftType>(lhs, rhs);
	}

	template<typename LeftType, typename RightType, int Dim>
	vec_base<LeftType, Dim>&
		operator*=(vec_base<LeftType, Dim>& lhs, const vec_base<RightType, Dim>& rhs) {
		return lhs = Math::mul<LeftType>(lhs, rhs);
	}


	template<typename LeftType, typename RightType, int Row, int Col>
	mat_base<LeftType, Row, Col>&
		operator+=(mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		return lhs = Math::add<LeftType>(lhs, rhs);
	}

	template<typename LeftType, typename RightType, int Row, int Col>
	mat_base<LeftType, Row, Col>&
		operator-=(mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		return lhs = Math::sub<LeftType>(lhs, rhs);
	}

	template<typename LeftType, typename RightType, int Row, int Col>
	mat_base<LeftType, Row, Col>&
		operator/=(mat_base<LeftType, Row, Col>& lhs, const mat_base<RightType, Row, Col>& rhs) {
		return lhs = Math::div<LeftType>(lhs, rhs);
	}

	template<typename LeftType, typename RightType, int Square>
	mat_base<LeftType, Square, Square>&
		operator*=(mat_base<LeftType, Square, Square>& lhs, const mat_base<RightType, Square, Square>& rhs) {
		return lhs = Math::mul<LeftType>(lhs, rhs);
	}

	template<typename LeftType, typename RightType>
	quat_base<LeftType>&
		operator+=(quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		return lhs = Math::add<LeftType>(lhs, rhs);
	}

	template<typename LeftType, typename RightType>
	quat_base<LeftType>&
		operator-=(quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		return lhs = Math::sub<LeftType>(lhs, rhs);
	}

	template<typename LeftType, typename RightType>
	quat_base<LeftType>&
		operator/=(quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		return lhs = Math::div<LeftType>(lhs, rhs);
	}

	template<typename LeftType, typename RightType>
	quat_base<LeftType>&
		operator*=(quat_base<LeftType>& lhs, const quat_base<RightType>& rhs) {
		return lhs = Math::mul<LeftType>(lhs, rhs);
	}

	template<typename LeftType, typename RightType, int Dim>
	vec_base<LeftType, Dim>&
		operator*=(vec_base<LeftType, Dim>& lhs, const mat_base<RightType, Dim, Dim>& rhs) {
		return lhs = Math::mul<LeftType>(lhs, rhs);
	}

	template<typename LeftType, typename RightType, int Dim>
	vec_base<LeftType, Dim>&
		operator*=(vec_base<LeftType, Dim>& lhs, const RightType& rhs) {
		return lhs = Math::mul<LeftType>(lhs, rhs);
	}

	template<typename LeftType, typename RightType, int Row, int Col>
	mat_base<LeftType, Row, Col>&
		operator*=(mat_base<LeftType, Row, Col>& lhs, const RightType& rhs) {
		return lhs = Math::mul<LeftType>(lhs, rhs);
	}

	template<typename LeftType, typename RightType>
	quat_base<LeftType>&
		operator*=(quat_base<LeftType>& lhs, const RightType& rhs) {
		return lhs = Math::mul<LeftType>(lhs, rhs);
	}

	template<typename LeftType, typename RightType, int Dim>
	vec_base<LeftType, Dim>&
		operator/=(vec_base<LeftType, Dim>& lhs, const RightType& rhs) {
		return lhs = Math::div<LeftType>(lhs, rhs);
	}

	template<typename LeftType, typename RightType, int Row, int Col>
	mat_base<LeftType, Row, Col>&
		operator/=(mat_base<LeftType, Row, Col>& lhs, const RightType& rhs) {
		return lhs = Math::div<LeftType>(lhs, rhs);
	}

	template<typename LeftType, typename RightType>
	quat_base<LeftType>&
		operator/=(quat_base<LeftType>& lhs, const RightType& rhs) {
		return lhs = Math::div<LeftType>(lhs, rhs);
	}

	//unary - , + 
	template<typename Type, int Dim>
	vec_base<Type, Dim>
		operator-(const vec_base<Type, Dim>& v) {
		return Math::minus<Type>(v);
	}

	template<typename Type, int Dim>
	const vec_base<Type, Dim>&
		operator+(const vec_base<Type, Dim>& v) {
		return v;
	}

	template<typename Type>
	quat_base<Type>
		operator-(const quat_base<Type>& v) {
		return Math::minus<Type>(v);
	}

	template<typename Type>
	const quat_base<Type>&
		operator+(const quat_base<Type>& v) {
		return v;
	}

	template<typename Type, int Row, int Col>
	mat_base<Type, Row, Col>
		operator-(const mat_base<Type, Row, Col>& v) {
		return Math::minus<Type>(v);
	}

	template<typename Type, int Row, int Col>
	const mat_base<Type, Row, Col>&
		operator+(const mat_base<Type, Row, Col>& v) {
		return v;
	}
}//namespace Math
}//namespace Zee