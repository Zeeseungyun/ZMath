#pragma once
#include "../Common/ZTypeTraits.h"
#include "../Common/ZCommonFunctions.h"

namespace Zee {
namespace Math {
namespace Quaternion {
	//use only mat3x3_t or mat4x4_t.
	template<typename RetType, typename Type, int Square>
	quat_t<RetType> make_rotation_matrix(const mat_base<Type, Square, Square>& m);

	template<typename Type, typename EpsT>
	bool is_identity(const quat_t<Type>& q);

	template<typename Type, typename EpsT>
	bool is_near_identity(const quat_t<Type>& q, const EpsT& e);

	template<typename RetType, typename Type>
	quat_t<RetType> make_conjugate(const quat_t<Type>& q);

	template<typename RetType, typename Type>
	quat_t<RetType> make_inverse(const quat_t<Type>& q);

	template<typename RetType, typename RadT>
	quat_t<RetType> make_rotation_yaw_pitch_roll(const RadT& yaw, const RadT& pitch, const RadT& roll);

	template<typename RetType, typename RadT>
	quat_t<RetType> make_rotation_yaw_pitch_roll(const vec3_t<RadT>& xYaw_yPitch_zRoll);

	template<typename RetType, typename Type, typename RadT>
	quat_t<RetType> make_rotation_axis(const vec3_t<Type>& axis, const RadT& rad);

	template<typename RetType, typename Type, typename RadT>
	quat_t<RetType> make_rotation_norm(const vec3_t<Type>& norm, const RadT& rad);

	template<typename RetType, typename Type>
	vec3_t<RetType> make_axis_x(const quat_t<Type>& q);
	template<typename RetType, typename Type>
	vec3_t<RetType> make_axis_y(const quat_t<Type>& q);
	template<typename RetType, typename Type>
	vec3_t<RetType> make_axis_z(const quat_t<Type>& q);
	template<typename RetType, typename Type>
	vec3_t<RetType> make_axis_x_transpose(const quat_t<Type>& q);
	template<typename RetType, typename Type>
	vec3_t<RetType> make_axis_y_transpose(const quat_t<Type>& q);
	template<typename RetType, typename Type>
	vec3_t<RetType> make_axis_z_transpose(const quat_t<Type>& q);

	//*return xYaw_yPitch_zRoll.
	template<typename RetType, typename Type>
	vec3_t<RetType> make_yaw_pitch_roll(const quat_t<Type>& q);

	template<typename RetType, typename Type>
	mat3x3_t<RetType> make_mat3x3(const quat_t<Type>& q);
	template<typename RetType, typename Type>
	mat4x4_t<RetType> make_mat4x4(const quat_t<Type>& q);

	template<typename RetType, typename Type, typename FloatT>
	quat_t<RetType> slerp(const quat_t<Type>& q0, const quat_t<Type>& q1, const FloatT& t);

	template<typename RetType, typename Type>
	void extract_axis_angle(const quat_t<Type>& q, vec3_t<RetType>& out_axis, RetType& out_rad);

}//namespace Quaternion

	template<typename ScalarType>
	struct quat_base {
		using scalar_type = ScalarType;

	public:
		struct constant;
		struct common;
		struct personal;

	public: // constructor
		constexpr quat_base(const quat_base&) = default;
		constexpr quat_base(quat_base&&) = default;

		constexpr quat_base() : data{ scalar_type() } {
			static_assert(is_arithmetic<scalar_type>::value, "scalar type must arithemtic type.");
		}

		template<typename ...Args>
		constexpr quat_base(scalar_type a, Args ...args) : data{ static_cast<scalar_type>(a), static_cast<scalar_type>(args)... } {
		}

		template<typename other_scalar>
		quat_base(const quat_base<other_scalar>& other) { this->operator=(other); }

	public: // assign operator
		constexpr quat_base& operator=(const quat_base&) = default;
		constexpr quat_base& operator=(quat_base&&) = default;

		template<typename other_scalar>
		quat_base& operator=(const quat_base<other_scalar>& other) {
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

	public: //quaternion common functions
		template<typename other_scalar>
		bool is_equal(const quat_t<other_scalar>& rhs) const { return common::is_equal(*this, rhs); }
		template<typename other_scalar>
		bool is_near_equal(const quat_t<other_scalar>& rhs, const scalar_type e) const { return common::is_near_equal(*this, rhs, e); }

		bool is_zero() const { return common::is_zero(*this); }
		bool is_near_zero(const scalar_type e) const { return common::is_near_zero(*this, e); }

		bool is_identity() const { return common::is_identity(*this); }
		bool is_near_identity(const scalar_type& e) const { return common::is_near_identity(*this, e); }

		scalar_type length() const { return common::length(*this); }
		scalar_type lengthSq() const { return common::lengthSq(*this); }

		quat_t<scalar_type> get_normalize() const { return common::make_normalize(*this); }
		void set_normalize() { *this = common::make_normalize(*this); }

		template<typename other_scalar>
		scalar_type distance(const quat_t<other_scalar>& rhs) const { return common::distance(*this, rhs); }
		template<typename other_scalar>
		scalar_type distanceSq(const quat_t<other_scalar>& rhs) const { return common::distanceSq(*this, rhs); }
		template<typename other_scalar>
		scalar_type dot(const quat_t<other_scalar>& rhs) const { return common::dot(*this, rhs); }
		template<typename other_scalar>
		bool in_bounds(const quat_t<other_scalar>& bound) const { return common::in_bounds(*this, bound); }

	public: //quaternion personal functions
		quat_t<scalar_type> get_conjugate() const { return personal::make_conjugate(*this); }
		quat_t<scalar_type> get_inverse() const { return personal::make_inverse(*this); }
		void set_conjugate() { *this = personal::make_conjugate(*this); }
		void set_inverse() { *this = personal::make_inverse(*this); }

		template<typename other_scalar>
		void set_rotation_yaw_pitch_roll(other_scalar yaw, other_scalar pitch, other_scalar roll) {
			*this = personal::make_rotation_yaw_pitch_roll(yaw, pitch, roll);
		}

		template<typename other_scalar>
		void set_rotation_yaw_pitch_roll(const vec3_t<other_scalar>& xYaw_yPitch_zRoll) {
			*this = personal::make_rotation_yaw_pitch_roll(xYaw_yPitch_zRoll);
		}

		template<typename other_scalar, typename rad_t>
		void set_rotation_axis(const vec3_t<other_scalar>& axis, const rad_t& rad) {
			*this = personal::make_rotation_axis(axis, rad);
		}

		template<typename other_scalar, typename rad_t>
		void set_rotation_norm(const vec3_t<other_scalar>& norm, const rad_t& rad) {
			*this = personal::make_rotation_norm(norm, rad);
		}

		// matrix has no scaling and must rotation matrix. t is ignored. (m[3][0~3])
		template<typename other_scalar>
		void set_rotation_matrix(const mat4x4_t<other_scalar>& m) {
			*this = personal::make_rotation_matrix(m);
		}

		// matrix has no scaling and must rotation matrix.
		template<typename other_scalar>
		void set_rotation_matrix(const mat3x3_t<other_scalar>& m) {
			*this = personal::make_rotation_matrix(m);
		}

		mat3x3_t<scalar_type> get_mat3x3() const { return personal::make_mat3x3(*this);	}
		mat4x4_t<scalar_type> get_mat4x4() const { return personal::make_mat4x4(*this); }

		vec3_t<scalar_type> get_yaw_pitch_roll() const { return personal::make_yaw_pitch_roll(*this); }

		const vec3_t<scalar_type>& get_axis() const { return v; }
		scalar_type get_angle() const { return (scalar_type)std::acos(mul<scalar_type>(2, s)); }

		void get_axis_angle(vec3_t<scalar_type>& out_angle, scalar_type& out_rad) const {
			personal::extract_axis_angle(out_angle, out_rad);
		}

	public:
		union {
			scalar_type data[4];
			struct { scalar_type x, y, z, w; };
			struct { vec_base<scalar_type, 3> v;  scalar_type s; };
		};
	};

	template<typename T>
	struct quat_base<T>::constant {
		static constexpr quat_t<T> zero     = quat_t<T>();
		static constexpr quat_t<T> one      = quat_t<T>(1, 1, 1, 1);
		static constexpr quat_t<T> identity = quat_t<T>(0, 0, 0, 1);
	};
	
	template<typename T>
	struct quat_base<T>::common {
		template<typename U>
		static bool is_equal(const quat_t<T>& lhs, const quat_t<U>& rhs) { return Math::is_equal(lhs, rhs); }
		static bool is_zero(const quat_t<T>& v) { return Math::is_zero(v); }
		static bool is_identity(const quat_t<T>& v) { return Quaternion::is_identity(v); }

		template<typename U, typename EpsT>
		static bool is_near_equal(const quat_t<T>& lhs, const quat_t<U>& rhs, const EpsT& e) { return Math::is_near_equal(lhs, rhs, e); }
		template<typename EpsT>
		static bool is_near_zero(const quat_t<T>& v, const EpsT& e) { return Math::is_near_zero(v, e); }
		template<typename EpsT>
		static bool is_near_identity(const quat_t<T>& v, const EpsT& e) { return Quaternion::is_near_identity(v, e); }

		static T length(const quat_t<T>& v) { return Math::length<T>(v); }
		static T lengthSq(const quat_t<T>& v) { return Math::lengthSq<T>(v); }

		template<typename U>
		static T distanceSq(const quat_t<T>& lhs, const quat_t<U>& rhs) {
			return Math::distanceSq<T>(lhs, rhs);
		}

		template<typename U>
		static T distance(const quat_t<T>& lhs, const quat_t<U>& rhs) {
			return Math::distance<T>(lhs, rhs);
		}

		template<typename U>
		static T dot(const quat_t<T>& lhs, const quat_t<U>& rhs) {
			return Math::dot<T>(lhs, rhs);
		}

		template<typename U>
		static bool in_bounds(const quat_t<T>& lhs, const quat_t<U>& bound) {
			return Math::is_less_equal(lhs, bound) && Math::is_greater_equal(lhs, Math::minus<U>(bound));
		}

		static quat_base make_normalize(const quat_base& v) { return Math::make_normalize<T>(v); }
	};

	template<typename T>
	struct quat_base<T>::personal {
		static quat_t<T> make_conjugate(const quat_t<T>& q) { return Quaternion::make_conjugate<T>(q); }
		static quat_t<T> make_inverse(const quat_t<T>& q) { return Quaternion::make_inverse<T>(q); }

		template<typename U>
		static quat_t<T> make_rotation_yaw_pitch_roll(const U& yaw, const U& pitch, const U& roll) {
			return Quaternion::make_rotation_yaw_pitch_roll<T>(yaw, pitch, roll);
		}

		template<typename U>
		static quat_t<T> make_rotation_yaw_pitch_roll(const vec3_t<U>& xYaw_yPitch_zRoll) {
			return Quaternion::make_rotation_yaw_pitch_roll<T>(xYaw_yPitch_zRoll);
		}

		template<typename U, typename RadT>
		static quat_t<T> make_rotation_axis(const vec3_t<U>& axis, const RadT& rad) {
			return Quaternion::make_rotation_axis<T>(axis, rad);
		}

		template<typename U, typename RadT>
		static quat_t<T> make_rotation_norm(const vec3_t<U>& norm, const RadT& rad) {
			return Quaternion::make_rotation_norm<T>(norm, rad);
		}

		// matrix has no scaling and must rotation matrix. t is ignored. (m[3][0~3])
		template<typename U>
		static quat_t<T> make_rotation_matrix(const mat4x4_t<U>& m) {
			return Quaternion::make_rotation_matrix<T>(m);
		}

		// matrix has no scaling and must rotation matrix.
		template<typename U>
		static quat_t<T> make_rotation_matrix(const mat3x3_t<U>& m) {
			return Quaternion::make_rotation_matrix<T>(m);
		}

		static vec3_t<T> make_axis_x(const quat_t<T>& q) { return Quaternion::make_axis_x<T>(q); }
		static vec3_t<T> make_axis_y(const quat_t<T>& q) { return Quaternion::make_axis_y<T>(q); }
		static vec3_t<T> make_axis_z(const quat_t<T>& q) { return Quaternion::make_axis_z<T>(q); }
		static vec3_t<T> make_axis_x_transpose(const quat_t<T>& q) { return Quaternion::make_axis_x_transpose<T>(q); }
		static vec3_t<T> make_axis_y_transpose(const quat_t<T>& q) { return Quaternion::make_axis_y_transpose<T>(q); }
		static vec3_t<T> make_axis_z_transpose(const quat_t<T>& q) { return Quaternion::make_axis_z_transpose<T>(q); }

		//*return xYaw_yPitch_zRoll.
		static vec3_t<T> make_yaw_pitch_roll(const quat_t<T>& q) { return Quaternion::make_yaw_pitch_roll<T>(q); }

		static mat3x3_t<T> make_mat3x3(const quat_t<T>& q) { return Quaternion::make_mat3x3<T>(q);}
		static mat4x4_t<T> make_mat4x4(const quat_t<T>& q) { return Quaternion::make_mat4x4<T>(q);}

		template<typename U, typename FloatT>
		static quat_t<T> slerp(const quat_t<U>& q0, const quat_t<U>& q1, const FloatT& t) {
			return Quaternion::slerp<T>(q0, q1, t);
		}

		template<typename U>
		static void extract_axis_angle(const quat_t<U>& q, vec3_t<T>& out_axis, T& out_rad) {
			return Quaternion::extract_axis_angle<T>(q, out_axis, out_rad);
		}
	};

namespace Quaternion {
	template<typename RetType, typename Type>
	quat_t<RetType> make_rotation_matrix(const mat_base<Type, 1, 1>& m) {
		static_assert(0, "invalid matrix type");
		return quat_t<RetType>();
	}

	template<typename RetType, typename Type>
	quat_t<RetType> make_rotation_matrix(const mat_base<Type, 2, 2>& m) {
		static_assert(0, "invalid matrix type");
		return quat_t<RetType>();
	}

	//use only mat3x3_t or mat4x4_t.
	template<typename RetType, typename Type, int Square>
	quat_t<RetType> make_rotation_matrix(const mat_base<Type, Square, Square>& m) {
		quat_t<RetType> q;

		RetType r22 = (RetType)m[2][2];
		if (r22 <= 0)  // x^2 + y^2 >= z^2 + w^2
		{
			RetType dif10 = sub<RetType>(m[1][1], m[0][0]);
			RetType omr22 = sub<RetType>(1, r22);
			if (dif10 <= 0)  // x^2 >= y^2
			{
				RetType fourXSqr = omr22 - dif10;
				RetType inv4x = div<RetType>(0.5, std::sqrt(fourXSqr));
				q[0] = fourXSqr * inv4x;
				q[1] = mul<RetType>(m[0][1] + m[1][0], inv4x);
				q[2] = mul<RetType>(m[0][2] + m[2][0], inv4x);
				q[3] = mul<RetType>(m[1][2] - m[2][1], inv4x);
			}
			else  // y^2 >= x^2
			{
				RetType fourYSqr = omr22 + dif10;
				RetType inv4y = div<RetType>(0.5, std::sqrt(fourYSqr));
				q[0] = mul<RetType>(m[0][1] + m[1][0], inv4y);
				q[1] = fourYSqr * inv4y;
				q[2] = mul<RetType>(m[1][2] + m[2][1], inv4y);
				q[3] = mul<RetType>(m[2][0] - m[0][2], inv4y);
			}
		}
		else  // z^2 + w^2 >= x^2 + y^2
		{
			RetType sum10 = add<RetType>(m[1][1], m[0][0]);
			RetType opr22 = add<RetType>(1, r22);
			if (sum10 <= 0)  // z^2 >= w^2
			{
				RetType fourZSqr = opr22 - sum10;
				RetType inv4z = div<RetType>(0.5, std::sqrt(fourZSqr));
				q[0] = mul<RetType>(m[0][2] + m[2][0], inv4z);
				q[1] = mul<RetType>(m[1][2] + m[2][1], inv4z);
				q[2] = fourZSqr * inv4z;
				q[3] = mul<RetType>(m[0][1] - m[1][0], inv4z);
			}
			else  // w^2 >= z^2
			{
				RetType fourWSqr = opr22 + sum10;
				RetType inv4w = div<RetType>(0.5, std::sqrt(fourWSqr));
				q[0] = mul<RetType>(m[1][2] - m[2][1], inv4w);
				q[1] = mul<RetType>(m[2][0] - m[0][2], inv4w);
				q[2] = mul<RetType>(m[0][1] - m[1][0], inv4w);
				q[3] = fourWSqr * inv4w;
			}
		}

		return make_normalize<RetType>(q);
	}

	template<typename Type, typename EpsT>
	bool is_identity(const quat_t<Type>& q) {
		return is_equal(q, quat_t<Type>::constant::identity);
	}

	template<typename Type, typename EpsT>
	bool is_near_identity(const quat_t<Type>& q, const EpsT& e) {
		return is_near_equal(q, quat_t<Type>::constant::identity, e);
	}

	template<typename RetType, typename Type>
	quat_t<RetType> make_conjugate(const quat_t<Type>& q) {
		return quat_t<RetType>{ -q[0], -q[1], -q[2], q[3] };
	}

	template<typename RetType, typename Type>
	quat_t<RetType> make_inverse(const quat_t<Type>& q) {
		return make_normalize<RetType>(make_conjugate<RetType>(q));
	}

	template<typename RetType, typename RadT>
	quat_t<RetType> make_rotation_yaw_pitch_roll(const RadT& yaw, const RadT& pitch, const RadT& roll) {
		RetType yaw_half   = mul<RetType>(yaw  , 0.5);
		RetType pitch_half = mul<RetType>(pitch, 0.5);
		RetType roll_half  = mul<RetType>(roll , 0.5);

		RetType c1 = std::cos(yaw_half);	RetType c2 = std::cos(roll_half); RetType c3 = std::cos(pitch_half);
		RetType s1 = std::sin(yaw_half);	RetType s2 = std::sin(roll_half); RetType s3 = std::sin(pitch_half);

		return quat_t<RetType>{
			-(s1 * s2 * c3 + c1 * c2 * s3),
			-(s1 * c2 * c3 - c1 * s2 * s3),
			-(c1 * s2 * c3 - s1 * c2 * s3),
			  c1 * c2 * c3 + s1 * s2 * s3
		};
	}

	template<typename RetType, typename RadT>
	quat_t<RetType> make_rotation_yaw_pitch_roll(const vec3_t<RadT>& xYaw_yPitch_zRoll) {
		return make_rotation_yaw_pitch_roll<RetType>(xYaw_yPitch_zRoll.x, xYaw_yPitch_zRoll.y, xYaw_yPitch_zRoll.z);
	}

	template<typename RetType, typename Type, typename RadT>
	quat_t<RetType> make_rotation_axis(const vec3_t<Type>& axis, const RadT& rad) {
		return make_rotation_norm(make_normalize<RetType>(axis), rad);
	}

	template<typename RetType, typename Type, typename RadT>
	quat_t<RetType> make_rotation_norm(const vec3_t<Type>& norm, const RadT& rad) {
		RetType rad_half = mul<RetType>(rad, 0.5);
		RetType c = std::cos(rad_half);
		RetType s = std::sin(rad_half);
		return quat_t<RetType>{
			mul<RetType>(norm[0], s),
			mul<RetType>(norm[1], s),
			mul<RetType>(norm[2], s),
			c
		};
	}

	template<typename RetType, typename Type>
	vec3_t<RetType> make_axis_x(const quat_t<Type>& q) {
		return vec3_t<RetType>{
			(RetType)1 - mul<RetType>(2 ,q.y * q.y + q.z * q.z) ,
			mul<RetType>(2, q.x * q.y + q.w * q.z),
			mul<RetType>(2, q.x * q.z - q.w * q.y)
		};
	}

	template<typename RetType, typename Type>
	vec3_t<RetType> make_axis_y(const quat_t<Type>& q) {		
		return vec3_t<RetType>{
			mul<RetType>(2, q.x * q.y - q.w * q.z),
			(RetType)1 - mul<RetType>(2, q.x * q.x + q.z * q.z),
			mul<RetType>(2, q.y * q.z + q.w * q.x)
		};
	}

	template<typename RetType, typename Type>
	vec3_t<RetType> make_axis_z(const quat_t<Type>& q) {
		return vec3_t<RetType>{
			mul<RetType>(2, q.x * q.z + q.w * q.y),
			mul<RetType>(2, q.y * q.z - q.w * q.x),
			(RetType)1 - mul<RetType>(2, q.x * q.x + q.y * q.y)
		};
	}
	
	template<typename RetType, typename Type>
	vec3_t<RetType> make_axis_x_transpose(const quat_t<Type>& q) {
		return vec3_t<RetType>{
			(RetType)1 - mul<RetType>(2, q.y * q.y + q.z * q.z),
			mul<RetType>(2, q.x * q.y - q.w * q.z),
			mul<RetType>(2, q.x * q.z + q.w * q.y)
		};
	}

	template<typename RetType, typename Type>
	vec3_t<RetType> make_axis_y_transpose(const quat_t<Type>& q) {		
		return vec3_t<RetType>{
			mul<RetType>(2, q.x * q.y + q.w * q.z),
			(RetType)1 - mul<RetType>(2, q.x * q.x + q.z * q.z),
			mul<RetType>(2, q.y * q.z - q.w * q.x)
		};
	}

	template<typename RetType, typename Type>
	vec3_t<RetType> make_axis_z_transpose(const quat_t<Type>& q) {		
		return vec3_t<RetType>{
			mul<RetType>(2, q.x * q.z - q.w * q.y),
			mul<RetType>(2, q.y * q.z + q.w * q.x),
			(RetType)1 - mul<RetType>(2, q.x * q.x + q.y * q.y)
		};
	}

	//*return xYaw_yPitch_zRoll.
	template<typename RetType, typename Type>
	vec3_t<RetType> make_yaw_pitch_roll(const quat_t<Type>& q) {
		quat_t<RetType> c = make_conjugate<RetType>(q);
		RetType xx = c.x*c.x; RetType yy = c.y*c.y; RetType zz = c.z*c.z; RetType ww = c.w*c.w;

		RetType y = (RetType)std::atan2((RetType)2 * (c.x * c.z + c.w * c.y), (-xx - yy + zz + ww));
		RetType p = (RetType)std::asin((RetType)-2 * (c.y * c.z - c.w * c.x) / (xx + yy + zz + ww));
		RetType r = (RetType)std::atan2((RetType)2 * (c.x * c.y + c.w * c.z), (-xx + yy - zz + ww));

		return vec3_t<RetType>{ y, p, r };
	}

	template<typename RetType, typename Type>
	mat3x3_t<RetType> make_mat3x3(const quat_t<Type>& q) {		
		RetType xx = mul<RetType>(q.x, q.x); RetType yy = mul<RetType>(q.y, q.y); RetType zz = mul<RetType>(q.z, q.z);
		RetType xy = mul<RetType>(q.x, q.y); RetType xz = mul<RetType>(q.x, q.z); RetType yz = mul<RetType>(q.y, q.z);
		RetType wx = mul<RetType>(q.w, q.x); RetType wy = mul<RetType>(q.w, q.y); RetType wz = mul<RetType>(q.w, q.z);

		return mat3x3_t<RetType>{
			vec3_t<RetType>{(RetType)1 - (RetType)2 * (yy + zz) , (RetType)2 * (xy + wz)		      , (RetType)2 * (xz - wy)		        },
			vec3_t<RetType>{(RetType)2 * (xy - wz)				, (RetType)1 - (RetType)2 * (xx + zz) , (RetType)2 * (yz + wx)		        },
			vec3_t<RetType>{(RetType)2 * (xz + wy)				, (RetType)2 * (yz - wx)			  , (RetType)1 - (RetType)2 * (xx + yy) }
		};
	}

	template<typename RetType, typename Type>
	mat4x4_t<RetType> make_mat4x4(const quat_t<Type>& q) {	
		RetType xx = mul<RetType>(q.x, q.x); RetType yy = mul<RetType>(q.y, q.y); RetType zz = mul<RetType>(q.z, q.z);
		RetType xy = mul<RetType>(q.x, q.y); RetType xz = mul<RetType>(q.x, q.z); RetType yz = mul<RetType>(q.y, q.z);
		RetType wx = mul<RetType>(q.w, q.x); RetType wy = mul<RetType>(q.w, q.y); RetType wz = mul<RetType>(q.w, q.z);

		return mat4x4_t<RetType>{
			vec4_t<RetType>{(RetType)1 - (RetType)2 * (yy + zz) , (RetType)2 * (xy + wz)		      , (RetType)2 * (xz - wy)		       , 0 },
			vec4_t<RetType>{(RetType)2 * (xy - wz)				, (RetType)1 - (RetType)2 * (xx + zz) , (RetType)2 * (yz + wx)		       , 0 },
			vec4_t<RetType>{(RetType)2 * (xz + wy)				, (RetType)2 * (yz - wx)			  , (RetType)1 - (RetType)2 * (xx + yy), 0 },
			vec4_t<RetType>{0,0,0,1}
		};
	}

	template<typename RetType, typename Type, typename FloatT>
	quat_t<RetType> slerp(const quat_t<Type>& q0, const quat_t<Type>& q1, const FloatT& t) {
		quat_t<RetType> qt;
		RetType cosom, sinom, scale0, scale1;

		// calc cosine
		cosom = dot<RetType>(q0, q1);

		// adjust signs (if necessary)
		if (is_less(cosom, 0)) {
			cosom = -cosom;
			qt = minus<RetType>(q1);
		}
		else {
			qt = static_cast<quat_t<RetType>>(q1);
		}

		// calculate coefficients
		
		if (is_greater(sub<RetType>(1, cosom), 0.00001)) {
			// standard case (slerp)
			RetType omega = std::acos(cosom);
			sinom = std::sin(omega);
			scale0 = div<RetType>(std::sin(sub<RetType>(1, t) * omega), sinom);
			scale1 = div<RetType>(std::sin(mul<RetType>(t, omega)), sinom);
		}
		else {
			// "from" and "to" quaternions are very close 
			// ... so we can do a linear interpolation
			scale0 = sub<RetType>(1, t);
			scale1 = (RetType)t;
		}

		return make_normalize(add<RetType>(mul<RetType>(q0, scale0), mul<RetType>(qt, scale1)));
	}

	template<typename RetType, typename Type>
	void extract_axis_angle(const quat_t<Type>& q, vec3_t<RetType>& out_axis, RetType& out_rad) {
		assign(out_axis, q.v);
		out_rad = mul<RetType>(2, std::acos(q.s));
	}

}//namespace Quaternion
}//namespace Math 
}//namespace Zee