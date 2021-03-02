#pragma once
#include "../Base/ZMat_Base.h"

namespace Zee {
namespace Math {
namespace Matrix4x4 {
	template<typename RetType, typename QuatType>
	mat4x4_t<RetType> make_rotation(const quat_t<QuatType>& q);

	template<typename RetType, typename RadT>
	mat4x4_t<RetType> make_rotation_x(const RadT& rad);
	template<typename RetType, typename RadT>
	mat4x4_t<RetType> make_rotation_y(const RadT& rad);
	template<typename RetType, typename RadT>
	mat4x4_t<RetType> make_rotation_z(const RadT& rad);

	template<typename RetType, typename RadT>
	mat4x4_t<RetType> make_rotation_yaw_pitch_roll(const RadT& yaw, const RadT& pitch, const RadT& roll);
	template<typename RetType, typename RadT>
	mat4x4_t<RetType> make_rotation_yaw_pitch_roll(const vec3_t<RadT>& xYaw_yPitch_zRoll);

	template<typename RetType, typename Type>
	mat4x4_t<RetType> make_scaling(const Type& x, const Type& y, const Type& z);
	template<typename RetType, typename Type>
	mat4x4_t<RetType> make_scaling(const vec3_t<Type>& xyz);

	template<typename RetType, typename Type>
	mat4x4_t<RetType> make_translation(const Type& x, const Type& y, const Type& z);
	template<typename RetType, typename Type>
	mat4x4_t<RetType> make_translation(const vec3_t<Type>& xyz);

	template<typename RetType, typename Type>
	mat4x4_t<RetType> make_transform(const vec3_t<Type>& scale, const quat_t<Type>& rotation, const vec3_t<Type>& translation);
	template<typename RetType, typename Type>
	mat4x4_t<RetType> make_inverse_transform(const vec3_t<Type>& scale, const quat_t<Type>& rotation, const vec3_t<Type>& translation);

	template<typename RetType, typename Type>
	bool decompose(const mat4x4_t<Type>& m, vec3_t<RetType>& out_scale, quat_t<RetType>& out_rotation, vec3_t<RetType>& out_translation);
	
	template<typename RetType, typename Type>
	vec3_t<RetType> extract_scaling_rapid(const mat4x4_t<Type>& m);
	template<typename RetType, typename Type>
	quat_t<RetType> extract_rotation_rapid(const mat4x4_t<Type>& m);
	template<typename RetType, typename Type>
	quat_t<RetType> extract_rotation_rapid(const mat4x4_t<Type>& m, vec3_t<RetType>& out_scale);
	template<typename RetType, typename Type>
	vec3_t<RetType> extract_translation(const mat4x4_t<Type>& m);

	template<typename RetType, typename Type>
	vec3_t<RetType> extract_axis_x(const mat4x4_t<Type>& m);
	template<typename RetType, typename Type>
	vec3_t<RetType> extract_axis_y(const mat4x4_t<Type>& m);
	template<typename RetType, typename Type>
	vec3_t<RetType> extract_axis_z(const mat4x4_t<Type>& m);

	template<typename RetType, typename Type>
	vec3_t<RetType> extract_yaw_pitch_roll(const mat4x4_t<Type>& m);
}//namespace Matrix4x4

	template<typename ScalarType>
	struct mat_base< ScalarType, 4, 4> {
		static constexpr int row = 4; static constexpr int col = 4;
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
			const vec_base<other_scalar, col>& r2 = vec_base<other_scalar, col>(),
			const vec_base<other_scalar, col>& r3 = vec_base<other_scalar, col>())
			: data{ static_cast<element_type>(r0), static_cast<element_type>(r1), static_cast<element_type>(r2), static_cast<element_type>(r3) }
		{
			static_assert(is_arithmetic<other_scalar>::value, "scalar type must arithemtic type.");
		}

		template<typename other_scalar>
		constexpr mat_base(
			const other_scalar& m00    , const other_scalar& m01 = 0, const other_scalar& m02 = 0, const other_scalar& m03 = 0,
			const other_scalar& m10 = 0, const other_scalar& m11 = 0, const other_scalar& m12 = 0, const other_scalar& m13 = 0,
			const other_scalar& m20 = 0, const other_scalar& m21 = 0, const other_scalar& m22 = 0, const other_scalar& m23 = 0,
			const other_scalar& m30 = 0, const other_scalar& m31 = 0, const other_scalar& m32 = 0, const other_scalar& m33 = 0 )
			: data{ element_type{static_cast<scalar_type>(m00), static_cast<scalar_type>(m01), static_cast<scalar_type>(m02), static_cast<scalar_type>(m03)},
					element_type{static_cast<scalar_type>(m10), static_cast<scalar_type>(m11), static_cast<scalar_type>(m12), static_cast<scalar_type>(m13)},
					element_type{static_cast<scalar_type>(m20), static_cast<scalar_type>(m21), static_cast<scalar_type>(m22), static_cast<scalar_type>(m23)}, 
					element_type{static_cast<scalar_type>(m30), static_cast<scalar_type>(m31), static_cast<scalar_type>(m32), static_cast<scalar_type>(m33)}}
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
				static_cast<vec_base<other_scalar, col>>(data[2]),
				static_cast<vec_base<other_scalar, col>>(data[3])
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

	public: // mat4x4_t personal functions
		template<typename other_scalar>
		void set_rotation(const quat_t<other_scalar>& q) { *this = personal::make_rotation(q); }
		template<typename rad_t>
		void set_rotation_x(const rad_t& rad) { *this = personal::make_rotation_x(rad); }
		template<typename rad_t>
		void set_rotation_y(const rad_t& rad) { *this = personal::make_rotation_y(rad); }
		template<typename rad_t>
		void set_rotation_z(const rad_t& rad) { *this = personal::make_rotation_z(rad); }
		template<typename rad_t>
		void set_rotation_yaw_pitch_roll(const rad_t& yaw, const rad_t& pitch, const rad_t& roll) { 
			*this = personal::make_rotation_yaw_pitch_roll(yaw, pitch, roll); 
		}
		template<typename rad_t>
		void set_rotation_yaw_pitch_roll(const vec3_t<rad_t>& xYaw_yPitch_zRoll) {
			*this = personal::make_rotation_yaw_pitch_roll(xYaw_yPitch_zRoll);
		}
		template<typename other_scalar>
		void set_scaling(const vec3_t<other_scalar>& xyz) {
			*this = personal::make_scaling(xyz);
		}
		template<typename other_scalar>
		void set_scaling(const other_scalar& x, const other_scalar& y, const other_scalar& z) {
			*this = personal::make_scaling(x, y, z);
		}
		template<typename other_scalar>
		void set_translation(const vec3_t<other_scalar>& xyz) {
			*this = personal::make_translation(xyz);
		}
		template<typename other_scalar>
		void set_translation(const other_scalar& x, const other_scalar& y, const other_scalar& z) {
			*this = personal::make_translation(x, y, z);
		}

		bool decompose(vec3_t<scalar_type>& out_scale, quat_t<scalar_type>& out_rotation, vec3_t<scalar_type>& out_translation) {
			return personal::decompose(out_scale, out_rotation, out_translation);
		}

		vec3_t<scalar_type> get_scaling_rapid() const {
			return personal::extract_scaling_rapid(*this);
		}

		quat_t<scalar_type> get_rotation_rapid(vec3_t<scalar_type>& out_scale) const {
			return personal::extract_rotation_rapid(*this, out_scale);
		}

		quat_t<scalar_type> get_rotation_rapid() const {
			return personal::extract_rotation_rapid(*this);
		}

		vec3_t<scalar_type> get_yaw_pitch_roll() const {
			return personal::extract_yaw_pitch_roll(*this);
		}

		const vec3_t<scalar_type>& get_translation() const { return data[3].v; }
		const vec3_t<scalar_type>& get_axis_x() const { return data[0].v; }
		const vec3_t<scalar_type>& get_axis_y() const { return data[1].v; }
		const vec3_t<scalar_type>& get_axis_z() const { return data[2].v; }

	public:
		element_type data[row];
		
	public: //misc
		struct constant;
		struct common;
		struct personal;
	};

	template<typename T>
	struct mat_base<T, 4, 4>::constant {
		static constexpr mat_base zero = mat_base();
		static constexpr mat_base identity = mat_base(
			vec4_t<T>::constant::basis[0],
			vec4_t<T>::constant::basis[1],
			vec4_t<T>::constant::basis[2],
			vec4_t<T>::constant::basis[3]
		);
	};

	template<typename T>
	struct mat_base<T, 4, 4>::common {
		template<typename other_scalar>
		static bool is_equal(const mat_base<T, 4, 4>& lhs, const mat_base<other_scalar, 4, 4>& rhs) { return Math::is_equal(lhs, rhs); }
		static bool is_zero(const mat_base<T, 4, 4>& v) { return Math::is_zero(v); }
		static bool is_identity(const mat_base<T, 4, 4>& v) { return Matrix::is_identity(v); }

		template<typename other_scalar, typename EpsT>
		static bool is_near_equal(const mat_base<T, 4, 4>& lhs, const mat_base<other_scalar, 4, 4>& rhs, EpsT e) { return Math::is_near_equal(lhs, rhs, e); }
		template<typename EpsT>
		static bool is_near_zero(const mat_base<T, 4, 4>& v, const EpsT& e) { return Math::is_near_zero(v, e); }
		template<typename EpsT>
		static bool is_near_identity(const mat_base<T, 4, 4>& v, const EpsT& e) { return Matrix::is_near_identity(v, e); }

		static T determinant(const mat_base<T, 4, 4>& m) { return Matrix::determinant<T>(m); }
		static mat_base<T, 4, 4> make_adjoint(const mat_base<T, 4, 4>& m) { return Matrix::make_adjoint<T>(m); }
		static mat_base<T, 4, 4> make_transpose(const mat_base<T, 4, 4>& m) { return Matrix::make_transpose<T>(m); }
		static mat_base<T, 4, 4> make_inverse(const mat_base<T, 4, 4>& m) { return Matrix::make_inverse<T>(m); }
		static mat_base<T, 4, 4> make_inverse(const mat_base<T, 4, 4>& m, T& det) { return Matrix::make_inverse<T>(m, det); }
	};

	template<typename T>
	struct mat_base<T, 4, 4>::personal {
		template<typename U>
		static mat_base<T, 4, 4> make_rotation(const quat_t<U>& q) { return Matrix4x4::make_rotation<T>(q); }
		template<typename RadT>
		static mat_base<T, 4, 4> make_rotation_x(const RadT& rad) { return Matrix4x4::make_rotation_x<T>(rad); }
		template<typename RadT>
		static mat_base<T, 4, 4> make_rotation_y(const RadT& rad) { return Matrix4x4::make_rotation_y<T>(rad); }
		template<typename RadT>
		static mat_base<T, 4, 4> make_rotation_z(const RadT& rad) { return Matrix4x4::make_rotation_z<T>(rad); }

		template<typename RadT>
		static mat_base<T, 4, 4> make_rotation_yaw_pitch_roll(const RadT& yaw, const RadT& pitch, const RadT& roll) {
			return Matrix4x4::make_rotation_yaw_pitch_roll<T>(yaw, pitch, roll);
		}

		template<typename RadT>
		static mat_base<T, 4, 4> make_rotation_yaw_pitch_roll(const vec3_t<RadT>& xYaw_yPitch_zRoll) {
			return Matrix4x4::make_rotation_yaw_pitch_roll<T>(xYaw_yPitch_zRoll);
		}

		template<typename U>
		static mat_base<T, 4, 4> make_scaling(const U& x, const U& y, const U& z) {
			return Matrix4x4::make_scaling<T>(x,y,z);
		}

		template<typename U>
		static mat_base<T, 4, 4> make_scaling(const vec3_t<U>& xyz) {
			return Matrix4x4::make_scaling<T>(xyz);
		}

		template<typename U>
		static mat_base<T, 4, 4> make_translation(const U& x, const U& y, const U& z) {
			return Matrix4x4::make_translation<T>(x, y, z);
		}

		template<typename U>
		static mat_base<T, 4, 4> make_translation(const vec3_t<U>& xyz) {
			return Matrix4x4::make_translation<T>(xyz);
		}

		template<typename U>
		static mat_base<T, 4, 4> make_transform(const vec3_t<U>& scale, const quat_t<U>& rotation, const vec3_t<U>& translation) {
			return Matrix4x4::make_transform<T>(scale, rotation, translation);
		}

		template<typename U>
		static mat_base<T, 4, 4> make_inverse_transform(const vec3_t<U>& scale, const quat_t<U>& rotation, const vec3_t<U>& translation) {
			return Matrix4x4::make_inverse_transform<T>(scale, rotation, translation);
		}

		template<typename U>
		static bool decompose(const mat_base<U, 4, 4>& m, vec3_t<T>& out_scale, quat_t<T>& out_rotation, vec3_t<T>& out_translation) {
			return Matrix4x4::decompose<T>(m, out_scale, out_rotation, out_translation);
		}

		template<typename U>
		static vec3_t<T> extract_scaling_rapid(const mat_base<U, 4, 4>& m) {
			return Matrix4x4::extract_scaling_rapid<T>(m);
		}

		template<typename U>
		static vec3_t<T> extract_rotation_rapid(const mat_base<U, 4, 4>& m , vec3_t<T>& out_scale) {
			return Matrix4x4::extract_rotation_rapid<T>(m, out_scale);
		}

		template<typename U>
		static vec3_t<T> extract_rotation_rapid(const mat_base<U, 4, 4>& m) {
			return Matrix4x4::extract_rotation_rapid<T>(m);
		}

		template<typename U>
		static vec3_t<T> extract_translation(const mat_base<U, 4, 4>& m) {
			return Matrix4x4::extract_translation<T>(m);
		}

		template<typename U>
		static vec3_t<T> extract_axis_x(const mat_base<U, 4, 4>& m) {
			return Matrix4x4::extract_axis_x<T>(m);
		}

		template<typename U>
		static vec3_t<T> extract_axis_y(const mat_base<U, 4, 4>& m) {
			return Matrix4x4::extract_axis_y<T>(m);
		}

		template<typename U>
		static vec3_t<T> extract_axis_z(const mat_base<U, 4, 4>& m) {
			return Matrix4x4::extract_axis_z<T>(m);
		}

		template<typename U>
		static vec3_t<T> extract_yaw_pitch_roll(const mat_base<U, 4, 4>& m) {
			return extract_yaw_pitch_roll<T>(m);
		}
	};
	
namespace Matrix4x4 {
namespace Decompose {
	template<typename T>
	constexpr T epsilon_decompose_v = static_cast<T>(0.0001);

	template<typename T>
	void rank_decompose(int& a, int& b, int& c, const T& x, const T& y, const T& z) {
		if (x < y) {
			if (y < z) {
				a = 2; b = 1; c = 0;
			} else {
				a = 1;
				if (x < z) {
					b = 2; c = 0;
				} else {
					b = 0; c = 2;
				}
			}
		} else {
			if (x < z) {
				a = 2; b = 0; c = 1;
			} else {
				a = 0;
				if (y < z) {
					b = 2; c = 1;
				} else {
					b = 1; c = 2;
				}
			}
		}
	}

}//namespace Decompose
	template<typename RetType, typename QuatType>
	mat4x4_t<RetType> make_rotation(const quat_t<QuatType>& q) {
		return Quaternion::make_mat4x4<RetType>(q);
	}

	template<typename RetType, typename RadT>
	mat4x4_t<RetType> make_rotation_x(const RadT& rad) {
		RetType c = (RetType)std::cos(rad);
		RetType s = (RetType)std::sin(rad);
		return mat4x4_t<RetType>{
			1, 0, 0, 0,
			0, c,-s, 0,
			0, s, c, 0,
			0, 0, 0, 1
		};
	}

	template<typename RetType, typename RadT>
	mat4x4_t<RetType> make_rotation_y(const RadT& rad) {
		RetType c = (RetType)std::cos(rad);
		RetType s = (RetType)std::sin(rad);
		return mat4x4_t<RetType>{
			c , 0, s, 0,
			0 , 1, 0, 0,
			-s, 0, c, 0,
			0 , 0, 0, 1
		};
	}

	template<typename RetType, typename RadT>
	mat4x4_t<RetType> make_rotation_z(const RadT& rad) {
		RetType c = (RetType)std::cos(rad);
		RetType s = (RetType)std::sin(rad);
		return mat4x4_t<RetType>{
			c,-s, 0, 0,
			s, c, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1
		};
	}

	template<typename RetType, typename RadT>
	mat4x4_t<RetType> make_rotation_yaw_pitch_roll(const RadT& yaw, const RadT& pitch, const RadT& roll) {
		return Quaternion::make_mat4x4<RetType>(Quaternion::make_rotation_yaw_pitch_roll<RetType>(yaw,pitch,roll));
	}

	template<typename RetType, typename RadT>
	mat4x4_t<RetType> make_rotation_yaw_pitch_roll(const vec3_t<RadT>& xYaw_yPitch_zRoll) {		
		return make_rotation_yaw_pitch_roll<RetType>(xYaw_yPitch_zRoll.x, xYaw_yPitch_zRoll.y, xYaw_yPitch_zRoll.z);
	}

	template<typename RetType, typename Type>
	mat4x4_t<RetType> make_scaling(const Type& x, const Type& y, const Type& z) {
		return mat4x4_t<RetType>{
			x, 0, 0, 0,
			0, y, 0, 0,
			0, 0, z, 0,
			0, 0, 0, 1,
		};
	}

	template<typename RetType, typename Type>
	mat4x4_t<RetType> make_scaling(const vec3_t<Type>& xyz) {
		return make_scaling<RetType>(xyz.x, xyz.y, xyz.z);
	}

	template<typename RetType, typename Type>
	mat4x4_t<RetType> make_translation(const Type& x, const Type& y, const Type& z) {
		return mat4x4_t<RetType>{
			vec4_t<RetType>::constant::basis[0],
			vec4_t<RetType>::constant::basis[1],
			vec4_t<RetType>::constant::basis[2],
			vec4_t<RetType>(x,y,z,1)
		};
	}

	template<typename RetType, typename Type>
	mat4x4_t<RetType> make_translation(const vec3_t<Type>& xyz) {
		return make_translation<RetType>(xyz.x, xyz.y, xyz.z);
	}

	template<typename RetType, typename Type>
	mat4x4_t<RetType> make_transform(const vec3_t<Type>& scale, const quat_t<Type>& rotation, const vec3_t<Type>& translation) {
		RetType xx = mul<RetType>(rotation.x, rotation.x); RetType yy = mul<RetType>(rotation.y, rotation.y); RetType zz = mul<RetType>(rotation.z, rotation.z);
		RetType xy = mul<RetType>(rotation.x, rotation.y); RetType xz = mul<RetType>(rotation.x, rotation.z); RetType yz = mul<RetType>(rotation.y, rotation.z);
		RetType wx = mul<RetType>(rotation.w, rotation.x); RetType wy = mul<RetType>(rotation.w, rotation.y); RetType wz = mul<RetType>(rotation.w, rotation.z);

		return mat4x4_t<RetType>{
			(RetType)1 - mul<RetType>(mul<RetType>(2, yy + zz), scale.x), mul<RetType>(mul<RetType>(2, xy + wz), scale.x)             , mul<RetType>(mul<RetType>(2, xz - wy), scale.x)            , 0 ,
			mul<RetType>(mul<RetType>(2, xy - wz), scale.y)             , (RetType)1 - mul<RetType>(mul<RetType>(2, xx + zz), scale.y), mul<RetType>(mul<RetType>(2, yz + wx), scale.y)            , 0 ,
			mul<RetType>(mul<RetType>(2, xz + wy), scale.z)             , mul<RetType>(mul<RetType>(2, yz - wx), scale.z)             ,(RetType)1 - mul<RetType>(mul<RetType>(2, xx + yy), scale.z), 0 ,
			translation.x, translation.y, translation.z, 1
		};

		//(RetType)1 - mul<RetType>(mul<RetType>(2, yy + zz), scale.x);
		//mul<RetType>(mul<RetType>(2, xy - wz), scale.y);
		//mul<RetType>(mul<RetType>(2, xz + wy), scale.z);
		//
		//mul<RetType>(mul<RetType>(2, xy + wz), scale.x);
		//(RetType)1 - mul<RetType>(mul<RetType>(2, xx + zz), scale.y);
		//mul<RetType>(mul<RetType>(2, yz - wx), scale.z);
		//
		//mul<RetType>(mul<RetType>(2, xz - wy), scale.x);
		//mul<RetType>(mul<RetType>(2, yz + wx), scale.y);
		//(RetType)1 - mul<RetType>(mul<RetType>(2, xx + yy), scale.z);
	}

	template<typename RetType, typename Type>
	mat4x4_t<RetType> make_inverse_transform(const vec3_t<Type>& scale, const quat_t<Type>& rotation, const vec3_t<Type>& translation) {

		vec3_t<RetType> invScale = div<RetType>(1, scale);
		quat_t<RetType> cQ = Quaternion::make_conjugate<RetType>(rotation);
		vec3_t<RetType> invPos = minus<RetType>(mul<RetType>(vec3_t<RetType>::personal::make_rotation(translation, cQ), invScale));
		
		RetType xx = cQ.x*cQ.x; RetType yy = cQ.y*cQ.y; RetType zz = cQ.z*cQ.z;
		RetType xy = cQ.x*cQ.y; RetType xz = cQ.x*cQ.z; RetType yz = cQ.y*cQ.z;
		RetType wx = cQ.w*cQ.x; RetType wy = cQ.w*cQ.y; RetType wz = cQ.w*cQ.z;
		
		return mat4x4_t<RetType>{
			((RetType)1 - (RetType)2 * (yy + zz))* invScale.x, ((RetType)2 * (xy + wz))* invScale.y             , ((RetType)2 * (xz - wy))* invScale.z    , 0 ,
			((RetType)2 * (xy - wz)) * invScale.x            , ((RetType)1 - (RetType)2 * (xx + zz))* invScale.y, ((RetType)2 * (yz + wx))* invScale.z    , 0 ,
			((RetType)2 * (xz + wy)) * invScale.x            , ((RetType)2 * (yz - wx))* invScale.y             , ((RetType)1 - 2 * (xx + yy))* invScale.z, 0 ,
			invPos.x, invPos.y, invPos.z, 1
		};
	}

	template<typename RetType, typename Type>
	bool decompose(const mat4x4_t<Type>& m, vec3_t<RetType>& out_scale, quat_t<RetType>& out_rotation, vec3_t<RetType>& out_translation) {
		out_translation = extract_translation<RetType>(m);
		out_scale = extract_scaling_rapid<RetType>(m);
		mat3x3_t<RetType> m33 = {
			m[0].v, m[1].v, m[2].v
		};

		int a, b, c;
		Decompose::rank_decompose(a, b, c, out_scale[0], out_scale[0], out_scale[0]);

		if (out_scale[a] < Decompose::epsilon_decompose_v<RetType>)
			m33[a] = vec3_t<RetType>::constant::basis[a];
		m33[a].set_normalize();

		if (out_scale[b] < Decompose::epsilon_decompose_v<RetType>) {
			RetType abs_x = std::abs(m33[a][0]);
			RetType abs_y = std::abs(m33[a][1]);
			RetType abs_z = std::abs(m33[a][2]);
			int aa, bb, cc;
			Decompose::rank_decompose(aa, bb, cc, abs_x, abs_y, abs_z);
			m33[b] = m33[a].cross(vec3_t<RetType>::constant::basis[cc]);
		}
		m33[b].set_normalize();

		if (out_scale[c] < Decompose::epsilon_decompose_v<RetType>)
			m33[c] = m33[a].cross(m33[b]);
		m33[c].set_normalize();

		RetType det = m33.determinant();
		if (det < 0) {
			out_scale[a] = -out_scale[a];
			m33[a] = Math::minus<RetType>(m33[a]);
			det = -det;
		}

		det -= (RetType)1;
		det *= det;

		if (Decompose::epsilon_decompose_v<RetType> < det)
			return false;

		out_rotation = Quaternion::make_rotation_matrix<RetType>(m33);
		return true;
	}

	template<typename RetType, typename Type>
	vec3_t<RetType> extract_scaling_rapid(const mat4x4_t<Type>& m) {
		return vec3_t<RetType>{
			Math::length<RetType>(m[0].v),
			Math::length<RetType>(m[1].v),
			Math::length<RetType>(m[2].v)
		};
	}

	template<typename RetType, typename Type>
	quat_t<RetType> extract_rotation_rapid(const mat4x4_t<Type>& m, vec3_t<RetType>& out_scale) {
		auto out_scale = extract_scaling_rapid<RetType>(m);
		mat3x3_t<RetType> m33 = {
			div<RetType>(m[0].v, out_scale[0]),
			div<RetType>(m[1].v, out_scale[1]),
			div<RetType>(m[2].v, out_scale[2])
		};
		return Quaternion::make_rotation_matrix<RetType>(m33);
	}

	template<typename RetType, typename Type>
	quat_t<RetType> extract_rotation_rapid(const mat4x4_t<Type>& m) {
		vec3_t<RetType> unused;
		return extract_rotation_rapid<RetType>(m, unused);
	}

	template<typename RetType, typename Type>
	vec3_t<RetType> extract_translation(const mat4x4_t<Type>& m) {
		return static_cast<vec3_t<RetType>>(m[3].v);
	}

	template<typename RetType, typename Type>
	vec3_t<RetType> extract_axis_x(const mat4x4_t<Type>& m) {
		return static_cast<vec3_t<RetType>>(m[0].v);
	}
	template<typename RetType, typename Type>
	vec3_t<RetType> extract_axis_y(const mat4x4_t<Type>& m) {
		return static_cast<vec3_t<RetType>>(m[1].v);
	}
	template<typename RetType, typename Type>
	vec3_t<RetType> extract_axis_z(const mat4x4_t<Type>& m) {
		return static_cast<vec3_t<RetType>>(m[2].v);
	}

	template<typename RetType, typename Type>
	vec3_t<RetType> extract_yaw_pitch_roll(const mat4x4_t<Type>& m) {
		return Quaternion::make_yaw_pitch_roll<RetType>(extract_rotation_rapid<RetType>(m));
	}
}//namespace Matrix4x4
}//namespace Math
}//namespace Zee