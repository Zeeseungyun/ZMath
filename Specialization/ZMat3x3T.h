#pragma once
#include "../Base/ZMat_Base.h"

namespace Zee {
namespace Math {
namespace Matrix3x3 {
	template<typename RetType, typename QuatType>
	mat3x3_t<RetType> make_matrix(const quat_t<QuatType>& q);
}//namespace Matrix3x3

	template<typename ScalarType>
	struct mat_base< ScalarType, 3, 3> {
		static constexpr int row = 3; static constexpr int col = 3;
		static constexpr bool is_square_matrix_v = row == col;

		using scalar_type = ScalarType;
		using element_type = vec_base<scalar_type, col>;

	public: // constructor
		constexpr mat_base(const mat_base&) = default;
		constexpr mat_base(mat_base&&) = default;

		constexpr mat_base() : data{ element_type() } {
			static_assert(is_arithmetic<scalar_type>::value, "scalar type must arithemtic type.");
		}

		template<typename other_scalar>
		constexpr mat_base(
			const vec_base<other_scalar, col>& r0, 
			const vec_base<other_scalar, col>& r1 = vec_base<other_scalar, col>(), 
			const vec_base<other_scalar, col>& r2 = vec_base<other_scalar, col>())
			: data{ static_cast<element_type>(r0), static_cast<element_type>(r1), static_cast<element_type>(r2) } 
		{
			static_assert(is_arithmetic<other_scalar>::value, "scalar type must arithemtic type.");
		}

		template<typename other_scalar>
		constexpr mat_base(
			const other_scalar& m00    , const other_scalar& m01 = 0, const other_scalar& m02 = 0,
			const other_scalar& m10 = 0, const other_scalar& m11 = 0, const other_scalar& m12 = 0,
			const other_scalar& m20 = 0, const other_scalar& m21 = 0, const other_scalar& m22 = 0) 
			: data{ element_type{static_cast<scalar_type>(m00), static_cast<scalar_type>(m01), static_cast<scalar_type>(m02)},
					element_type{static_cast<scalar_type>(m10), static_cast<scalar_type>(m11), static_cast<scalar_type>(m12)},
					element_type{static_cast<scalar_type>(m20), static_cast<scalar_type>(m21), static_cast<scalar_type>(m22)}, } 
		{
			static_assert(is_arithmetic<other_scalar>::value, "scalar type must arithemtic type.");
		}

		template<typename other_scalar>
		mat_base(const mat_base<other_scalar, row, col>& other) { this->operator=(other); }
		
	public: // assign operator
		mat_base& operator=(const mat_base&) = default;
		mat_base& operator=(mat_base&&) = default;

		template<typename other_scalar>
		mat_base& operator=(const mat_base<other_scalar, row, col>& other) {
			return Math::assign(*this, other);
		}

	public: // coversions
		operator element_type*() { return data; }
		operator const element_type*() const { return data; }
		
		template<typename other_scalar>
		explicit operator mat_base<other_scalar, row, col>() const {
			return mat_base<other_scalar, row, col>{
				static_cast<vec_base<other_scalar, col>>(data[0]),
				static_cast<vec_base<other_scalar, col>>(data[1]),
				static_cast<vec_base<other_scalar, col>>(data[2])
			};
		}

	public: // common member functions
		template<typename other_scalar>
		bool is_equal(const mat_base<other_scalar, row, col>& rhs) const { return common::is_equal(*this, rhs); }
		bool is_zero() const { return common::is_zero(*this); }
		bool is_identity() const { return common::is_identity(*this); }

		template<typename other_scalar, typename EpsT>
		bool is_near_equal(const mat_base<other_scalar, row, col>& rhs, const EpsT& e) const { return common::is_near_equal(*this, rhs, e); }
		template<typename EpsT>
		bool is_near_zero(const EpsT& e) const { return common::is_near_zero(*this, e); }
		template<typename EpsT>
		bool is_near_identity(const EpsT& e) const { return common::is_near_identity(*this, e); }

		mat_base<scalar_type, col, row> get_transpose() const { return common::make_transpose(*this); }
		mat_base<scalar_type, row, col> get_inverse(scalar_type& det) const { return common::make_inverse(*this, det); }
		mat_base<scalar_type, row, col> get_inverse()  const { return common::make_inverse(*this); }
		mat_base<scalar_type, row, col> get_adjoint() const { return common::make_adjoint(*this); }

		void set_transpose() { *this = common::make_transpose(*this); }
		void set_inverse(scalar_type& det) { *this = common::make_inverse(*this, det); }
		void set_inverse() { *this = common::make_inverse(*this); }
		void set_adjoint() { *this = common::make_adjoint(*this); }

		scalar_type determinant() const { return common::determinant(*this); }

	public: // mat3x3_t personal functions
		template<typename other_scalar>
		void set_matrix(const quat_t<other_scalar>& q) { *this = personal::make_matrix(q); }

	public:
		element_type data[row];

	public: //misc
		struct constant;
		struct common;
		struct personal;
	};

	template<typename T>
	struct mat_base<T, 3, 3>::constant {
		static constexpr mat_base zero = mat_base();
		static constexpr mat_base identity = mat_base(
			1, 0, 0, 
			0, 1, 0,
			0, 0, 1
		);
	};

	template<typename T>
	struct mat_base<T, 3, 3>::common {
		template<typename other_scalar>
		static bool is_equal(const mat_base<T, 3, 3>& lhs, const mat_base<other_scalar, 3, 3>& rhs) { return Math::is_equal(lhs, rhs); }
		static bool is_zero(const mat_base<T, 3, 3>& v) { return Math::is_zero(v); }
		static bool is_identity(const mat_base<T, 3, 3>& v) { return Matrix::is_identity(v); }

		template<typename other_scalar, typename EpsT>
		static bool is_near_equal(const mat_base<T, 3, 3>& lhs, const mat_base<other_scalar, 3, 3>& rhs, EpsT e) { return Math::is_near_equal(lhs, rhs, e); }
		template<typename EpsT>
		static bool is_near_zero(const mat_base<T, 3, 3>& v, const EpsT& e) { return Math::is_near_zero(v, e); }
		template<typename EpsT>
		static bool is_near_identity(const mat_base<T, 3, 3>& v, const EpsT& e) { return Matrix::is_near_identity(v, e); }

		static T determinant(const mat_base<T, 3, 3>& m) { return Matrix::determinant<T>(m); }
		static mat_base<T, 3, 3> make_adjoint(const mat_base<T, 3, 3>& m) { return Matrix::make_adjoint<T>(m); }
		static mat_base<T, 3, 3> make_transpose(const mat_base<T, 3, 3>& m) { return Matrix::make_transpose<T>(m); }
		static mat_base<T, 3, 3> make_inverse(const mat_base<T, 3, 3>& m) { return Matrix::make_inverse<T>(m); }
		static mat_base<T, 3, 3> make_inverse(const mat_base<T, 3, 3>& m, T& det) { return Matrix::make_inverse<T>(m, det); }
	};

	template<typename T>
	struct mat_base<T, 3, 3>::personal {
		template<typename U>
		static mat_base<T, 3, 3> make_matrix(const quat_t<U>& q) { return Matrix3x3::make_matrix<T>(q); }
	};
	
namespace Matrix3x3 {
	template<typename RetType, typename QuatType>
	mat3x3_t<RetType> make_matrix(const quat_t<QuatType>& q) {
		return Quaternion::make_mat3x3<RetType>(q);
	}

}//namespace Matrix3x3
}//namespace Math
}//namespace Zee