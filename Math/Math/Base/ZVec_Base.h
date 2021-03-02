#pragma once
#include "../Common/ZTypeTraits.h"
#include "../Common/ZCommonFunctions.h"

namespace Zee {
namespace Math {
namespace Vector {
	template<typename Type, int Dim>
	bool is_unit(const vec_base<Type, Dim>& v);

	template<typename Type, int Dim, typename EpsT>
	bool is_near_unit(const vec_base<Type, Dim>& v, const EpsT& e);
}//namespace Vector

	template<typename ScalarType, int Dimension>
	struct vec_base {
		using scalar_type = ScalarType;
		static constexpr int dim = Dimension;

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
		vec_base(const vec_base<other_scalar, dim>& other) { this->operator=(other); }

	public: // assign operator
		vec_base& operator=(const vec_base&) = default;
		vec_base& operator=(vec_base&&) = default;

		template<typename other_scalar>
		vec_base& operator=(const vec_base<other_scalar, dim>& other) {
			return Math::assign(*this, other);
		}

	public: // coversions
		operator scalar_type*() { return data; }
		operator const scalar_type*() const { return data; }

		template<typename other_scalar>
		explicit operator vec_base<other_scalar, dim>() const {
			vec_base<other_scalar, dim> ret;
			for (int i = 0; i != dim; i++) { 
				ret[i] = static_cast<other_scalar>(data[i]); 
			}
			return ret;
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

	public:
		scalar_type data[dim];

	public:
		struct constant;
		struct common;
		struct personal;
	};

	template<typename T, int Dim>
	struct vec_base<T, Dim>::constant {
		static constexpr vec_base zero = vec_base();
	};

	template<typename T, int Dim>
	struct vec_base<T, Dim>::common {
		template<typename U>
		static bool is_equal(const vec_base<T, Dim>& lhs, const vec_base<U, Dim>& rhs) { return Math::is_equal(lhs, rhs); }
		static bool is_zero(const vec_base<T, Dim>& v) { return Math::is_zero(v); }
		static bool is_unit(const vec_base<T, Dim>& v) { return Vector::is_unit(v); }

		template<typename U, typename EpsT>
		static bool is_near_equal(const vec_base<T, Dim>& lhs, const vec_base<U, Dim>& rhs, EpsT e) { return Math::is_near_equal(lhs, rhs, e); }
		template<typename EpsT>
		static bool is_near_zero(const vec_base<T, Dim>& v, const EpsT& e) { return Math::is_near_zero(v, e); }
		template<typename EpsT>
		static bool is_near_unit(const vec_base<T, Dim>& v, const EpsT& e) { return Vector::is_near_unit(v,e); }

		static T length(const vec_base<T, Dim>& v) { return Math::length<T>(v); }
		static T lengthSq(const vec_base<T, Dim>& v) { return Math::lengthSq<T>(v); }

		template<typename U>
		static T distanceSq(const vec_base<T, Dim>& lhs, const vec_base<U, Dim>& rhs) { 
			return Math::distanceSq<T>(lhs, rhs); 
		}

		template<typename U>
		static T distance(const vec_base<T, Dim>& lhs, const vec_base<U, Dim>& rhs) {
			return Math::distance<T>(lhs, rhs);
		}

		template<typename U>
		static T dot(const vec_base<T, Dim>& lhs, const vec_base<U, Dim>& rhs) {
			return Math::dot<T>(lhs, rhs);
		}

		template<typename U>
		static bool in_bounds(const vec_base<T, Dim>& lhs, const vec_base<U, Dim>& bound) {
			return Math::is_less_equal(lhs, bound) && Math::is_greater_equal(lhs, Math::minus<U>(bound));
		}

		static vec_base<T, Dim> make_normalize(const vec_base<T, Dim>& v) {	
			return Math::make_normalize<T>(v);
		}
	};

	template<typename T, int Dimension>
	struct vec_base<T, Dimension>::personal {
		//empty.
	};
namespace Vector {
	template<typename Type, int Dim>
	bool is_unit(const vec_base<Type, Dim>& v) {
		return is_equal(Math::length<Type>(v), 1);
	}

	template<typename Type, int Dim, typename EpsT>
	bool is_near_unit(const vec_base<Type, Dim>& v, const EpsT& e) {
		return is_near_equal(Math::length<Type>(v), 1, e);
	}
}//namespace Vector 
}//namespace Math 
}//namespace Zee