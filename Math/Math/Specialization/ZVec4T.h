#pragma once
#include "../Base/ZVec_Base.h"

namespace Zee {
namespace Math {
namespace Vector4 {
	template<typename RetType, typename VecType, typename MatType>
	vec4_t<RetType> make_transform(const vec4_t<VecType>& v, const mat4x4_t<MatType>& m);

	template<typename RetType, typename VecType, typename MatType>
	vec4_t<RetType> make_transform_norm(const vec4_t<VecType>& v, const mat4x4_t<MatType>& m);
}//namespace Vector4

	template<typename ScalarType>
	struct vec_base<ScalarType, 4> {
		using scalar_type = ScalarType;
		static constexpr int dim = 4;

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
		vec_base(const vec4_t<other_scalar>& other) { this->operator=(other); }

	public: // assign operator
		vec_base& operator=(const vec_base&) = default;
		vec_base& operator=(vec_base&&) = default;

		template<typename other_scalar>
		vec_base& operator=(const vec4_t<other_scalar>& other) {
			return Math::assign(*this, other);
		}

	public: // coversions
		operator scalar_type*() { return data; }
		operator const scalar_type*() const { return data; }
		
		template<typename other_scalar>
		explicit operator quat_base<other_scalar>() const {
			return quat_base<other_scalar> {
				static_cast<other_scalar>(data[0]), static_cast<other_scalar>(data[1]),
				static_cast<other_scalar>(data[2]), static_cast<other_scalar>(data[3])
			};
		}

		template<typename other_scalar>
		explicit operator vec4_t<other_scalar>() const {
			return vec4_t<other_scalar> {
				static_cast<other_scalar>(data[0]), static_cast<other_scalar>(data[1]),
				static_cast<other_scalar>(data[2]), static_cast<other_scalar>(data[3])
			};
		}

	public: // common member functions
		template<typename other_scalar>
		bool is_equal(const vec_base<other_scalar, dim>& rhs) const { return common::is_equal(*this, rhs); }
		bool is_zero() const { return common::is_zero(*this); }
		bool is_unit() const { return common::is_unit(*this); }

		template<typename other_scalar, typename EpsT>
		bool is_near_equal(const vec_base<other_scalar, dim>& rhs, const EpsT& e) const { return common::is_near_equal(*this, rhs, e); }
		template<typename EpsT>
		bool is_near_zero(const EpsT& e) const { return common::is_near_zero(*this, e); }
		template<typename EpsT>
		bool is_near_unit(const EpsT& e) const { return common::is_near_unit(*this, e); }

		scalar_type length() const { return common::length(*this); }
		scalar_type lengthSq() const { return common::lengthSq(*this); }

		vec_base get_normalize() const { return common::make_normalize(*this); }
		void set_normalize() { *this = common::make_normalize(*this); }

		template<typename other_scalar>
		scalar_type distance(const vec_base< other_scalar, dim>& rhs) const { return common::distance(*this, rhs); }
		template<typename other_scalar>
		scalar_type distanceSq(const vec_base< other_scalar, dim>& rhs) const { return common::distanceSq(*this, rhs); }
		template<typename other_scalar>
		scalar_type dot(const vec_base<other_scalar, dim>& rhs) const { return common::dot(*this, rhs); }

		template<typename other_scalar>
		bool in_bounds(const vec_base<other_scalar, dim>& bound) const { return common::in_bounds(*this, bound); }

	public: // vec4_t personal functions
		template<typename other_scalar>
		vec_base get_transform(const mat4x4_t<other_scalar>& m) { return personal::make_transform(*this, m); }
		template<typename other_scalar>
		void set_transform(const mat4x4_t<other_scalar>& m) { *this = personal::make_transform(*this, m); }

		template<typename other_scalar>
		vec_base get_transform_normal(const mat4x4_t<other_scalar>& m) { return personal::make_transform_normal(*this, m); }
		template<typename other_scalar>
		void set_transform_normal(const mat4x4_t<other_scalar>& m) { *this = personal::make_transform_normal(*this, m); }

	public:
		union {
			scalar_type data[dim];
			struct {
				scalar_type x, y, z, w;
			};
			struct {
				vec3_t<scalar_type> v;
				scalar_type w;
			};
		};

	public:
		struct constant;
		struct common;
		struct personal;
	};

	template<typename T>
	struct vec_base<T, 4>::constant {
		static constexpr vec4_t<T> zero   = vec4_t<T>();
		static constexpr vec4_t<T> one    = vec4_t<T>(1, 1, 1, 1);
		static constexpr vec4_t<T> unit_x = vec4_t<T>(1, 0, 0, 0);
		static constexpr vec4_t<T> unit_y = vec4_t<T>(0, 1, 0, 0);
		static constexpr vec4_t<T> unit_z = vec4_t<T>(0, 0, 1, 0);
		static constexpr vec4_t<T> unit_w = vec4_t<T>(0, 0, 0, 1);
		static constexpr vec4_t<T> basis[4] = { unit_x, unit_y, unit_z, unit_w };
	};

	template<typename T>
	struct vec_base<T, 4>::common {
		template<typename U>
		static bool is_equal(const vec_base<T, 4>& lhs, const vec_base<U, 4>& rhs) { return Math::is_equal(lhs, rhs); }
		static bool is_zero(const vec_base<T, 4>& v) { return Math::is_zero(v); }
		static bool is_unit(const vec_base<T, 4>& v) { return Vector::is_unit(v); }

		template<typename U, typename EpsT>
		static bool is_near_equal(const vec_base<T, 4>& lhs, const vec_base<U, 4>& rhs, const EpsT& e) { return Math::is_near_equal(lhs, rhs, e); }
		template<typename EpsT>
		static bool is_near_zero(const vec_base<T, 4>& v, const EpsT& e) { return Math::is_near_zero(v, e); }
		template<typename EpsT>
		static bool is_near_unit(const vec_base<T, 4>& v, const EpsT& e) { return Vector::is_near_unit(v, e); }

		static T length(const vec_base<T, 4>& v) { return Math::length<T>(v); }
		static T lengthSq(const vec_base<T, 4>& v) { return Math::lengthSq<T>(v); }

		template<typename U>
		static T distanceSq(const vec_base<T, 4>& lhs, const vec_base<U, 4>& rhs) {
			return Math::distanceSq<T>(lhs, rhs);
		}

		template<typename U>
		static T distance(const vec_base<T, 4>& lhs, const vec_base<U, 4>& rhs) {
			return Math::distance<T>(lhs, rhs);
		}

		template<typename U>
		static T dot(const vec_base<T, 4>& lhs, const vec_base<U, 4>& rhs) {
			return Math::dot<T>(lhs, rhs);
		}

		template<typename U>
		static bool in_bounds(const vec_base<T, 4>& lhs, const vec_base<U, 4>& bound) {
			return Math::is_less_equal(lhs, bound) && Math::is_greater_equal(lhs, Math::minus<U>(bound));
		}

		static vec_base<T, 4> make_normalize(const vec_base<T, 4>& v) {
			return Math::make_normalize<T>(v);
		}
	};

	template<typename T>
	struct vec_base<T, 4>::personal {		
		template<typename U>
		static vec4_t<T> make_transform(const vec4_t<T>& v, const mat4x4_t<U>& m) {
			return Vector4::make_transform<T>(v, m);
		}
		template<typename U>
		static vec4_t<T> make_transform_normal(const vec4_t<T>& v, const mat4x4_t<U>& m) {
			return Vector4::make_transform_norm<T>(v, m);
		}
	};

namespace Vector4 {
	template<typename RetType, typename VecType, typename MatType>
	vec4_t<RetType> make_transform(const vec4_t<VecType>& v, const mat4x4_t<MatType>& m) {
		return mul<RetType>(v, m);
	}

	template<typename RetType, typename VecType, typename MatType>
	vec4_t<RetType> make_transform_norm(const vec4_t<VecType>& v, const mat4x4_t<MatType>& m) {
		return vec4_t<RetType>{
			mul<RetType>(v[0], m[0][0]) + mul<RetType>(v[1], m[1][0]) + mul<RetType>(v[2], m[2][0]),
			mul<RetType>(v[0], m[0][1]) + mul<RetType>(v[1], m[1][1]) + mul<RetType>(v[2], m[2][1]),
			mul<RetType>(v[0], m[0][2]) + mul<RetType>(v[1], m[1][2]) + mul<RetType>(v[2], m[2][2]),
			(RetType)v[3],
		};
	}
}//namespace Vector4 
}//namespace Math 
}//namespace Zee