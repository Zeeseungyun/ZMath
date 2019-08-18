#pragma once
#include "../Base/ZVec_Base.h"

namespace Zee {
namespace Math {
namespace Vector2 {
	template<typename RetType, typename Type>
	vec2_t<RetType> tripleProduct(const vec2_t<Type>& v0, const vec2_t<Type>& v1, const vec2_t<Type>& v2);

	template<typename RetType, typename LeftType, typename RightType>
	RetType cross(const vec2_t<LeftType>& lhs, const vec2_t<RightType>& rhs);

	template<typename RetType, typename VecType, typename MatType>
	vec2_t<RetType> make_transform(const vec2_t<VecType>& v, const mat3x3_t<MatType>& m);

	template<typename RetType, typename VecType, typename MatType>
	vec2_t<RetType> make_transform_normal(const vec2_t<VecType>& v, const mat3x3_t<MatType>& m);
}//namespace Vector2 

	template<typename ScalarType>
	struct vec_base<ScalarType, 2> {
		using scalar_type = ScalarType;
		static constexpr int dim = 2;

	public: // constructor
		constexpr vec_base(const vec_base&) = default;
		constexpr vec_base(vec_base&&) = default;

		constexpr vec_base() : data{ scalar_type() } {
			static_assert(is_arithmetic<scalar_type>::value, "scalar type must arithemtic type.");
		}

		template<typename ...Args>
		constexpr vec_base(scalar_type a, Args ...args) : data{ static_cast<scalar_type>(a), static_cast<scalar_type>(args)... } {
		}

		template<typename other_scalar>
		vec_base(const vec2_t<other_scalar>& other) { this->operator=(other); }

	public://assign operator
		constexpr vec_base& operator=(const vec_base&) = default;
		constexpr vec_base& operator=(vec_base&&) = default;

		template<typename other_scalar>
		vec_base& operator=(const vec2_t<other_scalar>& other) {
			return Math::assign(*this, other);
		}

	public: // coversions
		operator scalar_type*() { return data; }
		operator const scalar_type*() const { return data; }

		template<typename other_scalar>
		explicit operator vec2_t<other_scalar>() const {
			return vec2_t<other_scalar>{
				static_cast<other_scalar>(data[0]), static_cast<other_scalar>(data[1])
			};
		}

	public: //vector common functions
		template<typename other_scalar>
		bool is_equal(const vec2_t<other_scalar>& rhs) const { return common::is_equal(*this, rhs); }
		bool is_zero() const { return common::is_zero(*this); }
		bool is_unit() const { return common::is_unit(*this); }

		template<typename other_scalar, typename EpsT>
		bool is_near_equal(const vec2_t<other_scalar>& rhs, const scalar_type e) const { return common::is_near_equal(*this, rhs, e); }
		template<typename EpsT>
		bool is_near_zero(const EpsT& e) const { return common::is_near_zero(*this, e); }
		template<typename EpsT>
		bool is_near_unit(const EpsT& e) const { return common::is_near_unit(*this, e); }

		scalar_type length() const { return common::length(*this); }
		scalar_type lengthSq() const { return common::lengthSq(*this); }

		vec_base get_normalize() const { return common::make_normalize(*this); }
		void set_normalize() { *this = common::make_normalize(*this); }

		template<typename other_scalar>
		scalar_type distance(const vec2_t<other_scalar>& rhs) const { return common::distance(*this, rhs); }
		template<typename other_scalar>
		scalar_type distanceSq(const vec2_t<other_scalar>& rhs) const { return common::distanceSq(*this, rhs); }
		template<typename other_scalar>
		scalar_type dot(const vec2_t<other_scalar>& rhs) const { return common::dot(*this, rhs); }
		template<typename other_scalar>
		bool in_bounds(const vec2_t<other_scalar>& bound) const { return common::in_bounds(*this, bound); }

	public: // vec2_t personal functions
		template<typename other_scalar>
		scalar_type cross(const vec2_t<other_scalar>& rhs) const { return personal::cross(*this, rhs); }

		template<typename other_scalar>
		vec_base get_transform(const mat3x3_t<other_scalar>& m) { return personal::make_transform(*this, m); }

		template<typename other_scalar>
		void set_transform(const mat3x3_t<other_scalar>& m) { *this = personal::make_transform(*this, m); }

		template<typename other_scalar>
		vec_base get_transform_normal(const mat3x3_t<other_scalar>& m) { return personal::make_transform_normal(*this, m); }

		template<typename other_scalar>
		void set_transform_normal(const mat3x3_t<other_scalar>& m) { *this = personal::make_transform_normal(*this, m); }

	public:
		union {
			scalar_type data[dim];
			struct {
				scalar_type x, y;
			};
		};

	public:
		struct constant;
		struct common;
		struct personal;
	};

	template<typename T>
	struct vec_base<T, 2>::constant {
		static constexpr vec_base zero   = vec_base();
		static constexpr vec_base one    = vec_base(1, 1);
		static constexpr vec_base unit_x = vec_base(1, 0);
		static constexpr vec_base unit_y = vec_base(0, 1);
		static constexpr vec_base basis[2] = { unit_x, unit_y };
	};

	template<typename T>
	struct vec_base<T, 2>::common {
		template<typename U>
		static bool is_equal(const vec_base<T, 2>& lhs, const vec_base<U, 2>& rhs) { return Math::is_equal(lhs, rhs); }
		static bool is_zero(const vec_base<T, 2>& v) { return Math::is_zero(v); }
		static bool is_unit(const vec_base<T, 2>& v) { return Vector::is_unit(v); }

		template<typename U, typename EpsT>
		static bool is_near_equal(const vec_base<T, 2>& lhs, const vec_base<U, 2>& rhs, const EpsT& e) { return Math::is_near_equal(lhs, rhs, e); }
		template<typename EpsT>
		static bool is_near_zero(const vec_base<T, 2>& v, const EpsT& e) { return Math::is_near_zero(v, e); }
		template<typename EpsT>
		static bool is_near_unit(const vec_base<T, 2>& v, const EpsT& e) { return Vector::is_near_unit(v, e); }

		static T length(const vec_base<T, 2>& v) { return Math::length<T>(v); }
		static T lengthSq(const vec_base<T, 2>& v) { return Math::lengthSq<T>(v); }

		template<typename U>
		static T distanceSq(const vec_base<T, 2>& lhs, const vec_base<U, 2>& rhs) {
			return Math::distanceSq<T>(lhs, rhs);
		}

		template<typename U>
		static T distance(const vec_base<T, 2>& lhs, const vec_base<U, 2>& rhs) {
			return Math::distance<T>(lhs, rhs);
		}

		template<typename U>
		static T dot(const vec_base<T, 2>& lhs, const vec_base<U, 2>& rhs) {
			return Math::dot<T>(lhs, rhs);
		}

		template<typename U>
		static bool in_bounds(const vec_base<T, 2>& lhs, const vec_base<U, 2>& bound) {
			return Math::is_less_equal(lhs, bound) && Math::is_greater_equal(lhs, Math::minus<U>(bound));
		}

		static vec_base<T, 2> make_normalize(const vec_base<T, 2>& v) { 
			return Math::make_normalize<T>(v); 
		}
	};


	template<typename T>
	struct vec_base<T, 2>::personal {
		template<typename U>
		static vec2_t<T> make_transform(const vec2_t<T>& v, const mat3x3_t<U>& m) {
			return Vector2::make_transform<T>(v, m);
		}

		template<typename U>
		static vec2_t<T> make_transform_normal(const vec2_t<T>& v, const mat3x3_t<U>& m) {
			return Vector2::make_transform_normal<T>(v, m);
		}

		template<typename U>
		static T cross(const vec2_t<T>& lhs, const vec2_t<U>& rhs) {
			return Vector2::cross<T>(lhs,rhs);
		}

		template<typename U>
		static vec2_t<T> tripleProduct(const vec2_t<U>& v0, const vec2_t<U>& v1, const vec2_t<U>& v2) {
			return Vector2::tripleProduct<T>(v0, v1, v2);
		}
	};

namespace Vector2 {
	template<typename RetType, typename LeftType, typename RightType>
	RetType cross(const vec2_t<LeftType>& lhs, const vec2_t<RightType>& rhs) {
		return sub<RetType>(mul<RetType>(lhs[0], rhs[1]), mul<RetType>(rhs[0], lhs[1]));
	}

	template<typename RetType, typename Type>
	vec2_t<RetType> tripleProduct(const vec2_t<Type>& v0, const vec2_t<Type>& v1, const vec2_t<Type>& v2) {
		return sub<RetType>(mul<RetType>(v1, dot<RetType>(v2, v0)), mul<RetType>(v0, dot<RetType>(v2, v1)));
	}

	template<typename RetType, typename VecType, typename MatType>
	vec2_t<RetType> make_transform(const vec2_t<VecType>& v, const mat3x3_t<MatType>& m) {
		return vec2_t<RetType>{
			mul<RetType>(v[0], m[0][0]) + mul<RetType>(v[1], m[1][0]) + (RetType)m[2][0],
			mul<RetType>(v[0], m[0][1]) + mul<RetType>(v[1], m[1][1]) + (RetType)m[2][1]
		};
	}

	template<typename RetType, typename VecType, typename MatType>
	vec2_t<RetType> make_transform_normal(const vec2_t<VecType>& v, const mat3x3_t<MatType>& m) {
		return vec2_t<RetType>{
			mul<RetType>(v[0], m[0][0]) + mul<RetType>(v[1], m[1][0]),
			mul<RetType>(v[0], m[0][1]) + mul<RetType>(v[1], m[1][1])
		};
	}
}//namespace Vector2
}//namespace Math 
}//namespace Zee