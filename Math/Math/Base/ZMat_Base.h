#pragma once
#include "../Common/ZTypeTraits.h"
#include "../Common/ZCommonFunctions.h"

namespace Zee {
namespace Math {
namespace Matrix {
	template<typename Type, int Row, int Col>
	bool is_identity(const mat_base<Type, Row, Col>& m);

	template<typename Type, int Row, int Col, typename EpsT>
	bool is_near_identity(const mat_base<Type, Row, Col>& m, const EpsT& e);

	template<typename RetType, typename Type, int Row, int Col>
	RetType determinant(const mat_base<Type, Row, Col>& m);

	template<typename RetType, typename Type, int Row, int Col>
	auto make_transpose(const mat_base<Type, Row, Col>& m);

	template<typename RetType, typename Type, int Row, int Col>
	auto make_inverse(const mat_base<Type, Row, Col>& m);

	template<typename RetType, typename Type, int Row, int Col>
	auto make_inverse(const mat_base<Type, Row, Col>& m, RetType& det);

	template<typename RetType, typename Type, int Row, int Col>
	auto make_adjoint(const mat_base<Type, Row, Col>& m);

}//namespace Matrix

	template<typename ScalarType, int Row, int Col>
	struct mat_base {
		static constexpr int row = Row; static constexpr int col = Col;
		static constexpr bool is_square_matrix_v = row == col;

		using scalar_type = ScalarType;
		using element_type = vec_base<scalar_type, col>;

	public: // constructor
		constexpr mat_base(const mat_base&) = default;
		constexpr mat_base(mat_base&&) = default;

		constexpr mat_base() : data{ element_type() } {
			static_assert(is_arithmetic<scalar_type>::value, "scalar type must arithemtic type.");
		}

		template<typename ...Args>
		constexpr mat_base(const element_type& a, Args ...args) : data{ static_cast<element_type>(a), static_cast<element_type>(args)... } {
		}

		template<typename other_scalar, typename ...Args>
		constexpr mat_base(const vec_base<other_scalar, col>& a, Args ...args) : data{ static_cast<element_type>(a), static_cast<element_type>(args)... } {
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
			mat_base<other_scalar, row, col> ret;
			Math::assign(ret, *this);
			return ret;
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

		//정방행렬만 가능. 정방행렬이 아닌 경우 0행렬로 반환됨.
		mat_base<scalar_type, row, col> get_inverse(scalar_type& det) const { return common::make_inverse(*this, det); }
		//정방행렬만 가능. 정방행렬이 아닌 경우 0행렬로 반환됨.
		mat_base<scalar_type, row, col> get_inverse()  const { return common::make_inverse(*this); }
		//정방행렬만 가능. 정방행렬이 아닌 경우 0행렬로 반환됨.
		mat_base<scalar_type, row, col> get_adjoint() const { return common::make_adjoint(*this); }

		//정방행렬만 가능.
		//void set_transpose() { *this = common::make_transpose(*this); } 
		//정방행렬만 가능. 정방행렬이 아닌 경우 0행렬로 설정됨.
		void set_inverse(scalar_type& det) { *this = common::make_inverse(*this, det); }
		//정방행렬만 가능. 정방행렬이 아닌 경우 0행렬로 설정됨.
		void set_inverse() { *this = common::make_inverse(*this); }
		//정방행렬만 가능. 정방행렬이 아닌 경우 0행렬로 설정됨.
		void set_adjoint() { *this = common::make_adjoint(*this); }

	public:
		element_type data[row];

	public: //misc
		struct constant;
		struct common;
		struct personal;
	};

	template<typename T, int Row, int Col>
	struct mat_base<T, Row, Col>::constant {
		static constexpr mat_base zero = mat_base();
	};

	template<typename T, int Row, int Col>
	struct mat_base<T, Row, Col>::common {
		template<typename other_scalar>
		static bool is_equal(const mat_base<T, Row, Col>& lhs, const mat_base<other_scalar, row, col>& rhs) { return Math::is_equal(lhs, rhs); }
		static bool is_zero(const mat_base<T, Row, Col>& v) { return Math::is_zero(v); }
		static bool is_identity(const mat_base<T, Row, Col>& v) { return Matrix::is_identity(v); }

		template<typename other_scalar, typename EpsT>
		static bool is_near_equal(const mat_base<T, Row, Col>& lhs, const mat_base<other_scalar, row, col>& rhs, EpsT e) { return Math::is_near_equal(lhs, rhs, e); }
		template<typename EpsT>
		static bool is_near_zero(const mat_base<T, Row, Col>& v, const EpsT& e) { return Math::is_near_zero(v, e); }
		template<typename EpsT>
		static bool is_near_identity(const mat_base<T, Row, Col>& v, const EpsT& e) { return Matrix::is_near_identity(v, e); }

		static T determinant(const mat_base<T, Row, Col>& m) { return Matrix::determinant<T>(m); }
		static mat_base<T, Row, Col> make_adjoint(const mat_base<T, Row, Col>& m) { return Matrix::make_adjoint<T>(m); }
		static mat_base<T, Col, Row> make_transpose(const mat_base<T, Row, Col>& m) { return Matrix::make_transpose<T>(m); }
		static mat_base<T, Row, Col> make_inverse(const mat_base<T, Row, Col>& m) { return Matrix::make_inverse<T>(m); }
		static mat_base<T, Row, Col> make_inverse(const mat_base<T, Row, Col>& m, T& det) { return Matrix::make_inverse<T>(m, det); }
	};

	template<typename T, int Row, int Col>
	struct mat_base<T, Row, Col>::personal {
		//empty.
	};

namespace Matrix {
	template<typename Type, int Row, int Col>
	bool is_identity(const mat_base<Type, Row, Col>& m) {
		return false;
	}
	
	template<typename Type, int Row, int Col, typename EpsT>
	bool is_near_identity(const mat_base<Type, Row, Col>& m, const EpsT& e) {
		return false;
	}

	template<typename Type, int Square>
	bool is_identity(const mat_base<Type, Square, Square>& m) {
		for (int r = 0; r != Square; ++r) {
			for (int c = 0; c != Square; ++c) {
				if (r == c) {
					if (Math::is_not_equal(m[r][c], 1)) return false;
				} else {
					if (Math::is_not_equal(m[r][c], 0)) return false;
				}
			}
		}
		return true;
	}

	template<typename Type, int Square, typename EpsT>
	bool is_near_identity(const mat_base<Type, Square, Square>& m, const EpsT& e) {
		for (int r = 0; r != Square; ++r) {
			for (int c = 0; c != Square; ++c) {
				if (r == c) {
					if (Math::is_not_near_equal(m[r][c], 1, e)) return false;
				} else {
					if (Math::is_not_near_equal(m[r][c], 0, e)) return false;
				}
			}
		}
		return true;
	}

	template<typename RetType, typename Type, int Row, int Col>
	RetType determinant(const mat_base<Type, Row, Col>& m) {
		//TODO:: error handling?
		return (RetType)0;
	}

	template<typename RetType, typename Type>
	RetType determinant(const mat_base<Type, 1, 1>& m) {
		return static_cast<RetType>(m[0][0]);
	}

	template<typename RetType, typename Type>
	RetType determinant(const mat2x2_t<Type>& m) {
		return static_cast<RetType>(
			Math::mul<RetType>(m[0][0], m[1][1]) - Math::mul<RetType>(m[0][1], m[1][0])
		);
	}

	template<typename RetType, typename Type>
	RetType determinant(const mat3x3_t<Type>& m) {
		return static_cast<RetType>(
			(Math::mul<RetType>(Math::mul<RetType>(m[0][0], m[1][1]), m[2][2]) + Math::mul<RetType>(Math::mul<RetType>(m[0][1], m[1][2]), m[2][0]) + Math::mul<RetType>(Math::mul<RetType>(m[0][2], m[1][0]), m[2][1])) -
			(Math::mul<RetType>(Math::mul<RetType>(m[0][0], m[1][2]), m[2][1]) + Math::mul<RetType>(Math::mul<RetType>(m[0][1], m[1][0]), m[2][2]) + Math::mul<RetType>(Math::mul<RetType>(m[0][2], m[1][1]), m[2][0]))
		);
	}

	template<typename RetType, typename Type>
	RetType determinant(const mat4x4_t<Type>& m) {
		return static_cast<RetType>(
			(Math::mul<RetType>(m[0][0], m[1][1]) - Math::mul<RetType>(m[1][0], m[0][1])) * (Math::mul<RetType>(m[2][2], m[3][3]) - Math::mul<RetType>(m[3][2], m[2][3])) -
			(Math::mul<RetType>(m[0][0], m[1][2]) - Math::mul<RetType>(m[1][0], m[0][2])) * (Math::mul<RetType>(m[2][1], m[3][3]) - Math::mul<RetType>(m[3][1], m[2][3])) +
			(Math::mul<RetType>(m[0][0], m[1][3]) - Math::mul<RetType>(m[1][0], m[0][3])) * (Math::mul<RetType>(m[2][1], m[3][2]) - Math::mul<RetType>(m[3][1], m[2][2])) +
			(Math::mul<RetType>(m[0][1], m[1][2]) - Math::mul<RetType>(m[1][1], m[0][2])) * (Math::mul<RetType>(m[2][0], m[3][3]) - Math::mul<RetType>(m[3][0], m[2][3])) -
			(Math::mul<RetType>(m[0][1], m[1][3]) - Math::mul<RetType>(m[1][1], m[0][3])) * (Math::mul<RetType>(m[2][0], m[3][2]) - Math::mul<RetType>(m[3][0], m[2][2])) +
			(Math::mul<RetType>(m[0][2], m[1][3]) - Math::mul<RetType>(m[1][2], m[0][3])) * (Math::mul<RetType>(m[2][0], m[3][1]) - Math::mul<RetType>(m[3][0], m[2][1]))
		);
	}

	template<typename RetType, typename Type, int Square>
	RetType determinant(const mat_base<Type, Square, Square>& m) {
		mat_base<Type, Square - 1, Square - 1> minorMat;
		Type* p = nullptr;

		RetType det = RetType();

		for (int outer = 0; outer != Square; ++outer) {
			p = &minorMat[0][0];
			for (int i = 1; i != Square; ++i) {
				for (int j = 0; j != Square; ++j) {
					if (outer == j) {
						continue;
					}
					*p++ = m[i][j];
				}
			}

			RetType d = determinant<RetType>(minorMat);
			d = (outer & 1) ? -d : d;

			det += Math::mul<RetType>(m[0][outer], d);
		}

		return det;
	}

	template<typename RetType, typename Type, int Row, int Col>
	RetType cofactor(const mat_base<Type, Row, Col>& m, int row, int col) {
		//TODO:: error handling?
		return (RetType)0;
	}

	template<typename RetType, typename Type>
	RetType cofactor(const mat_base<Type, 1, 1>& m, int row, int col) {
		return static_cast<RetType>(m[row][col]);
	}

	template<typename RetType, typename Type, int Square>
	RetType cofactor(const mat_base<Type, Square, Square>& m, int row, int col) {
		static_assert(Square > 1, __FUNCTION__ " : invalid Component Count.");
		mat_base<Type, Square - 1, Square - 1> minorMat;
		Type* p = &minorMat[0][0];

		for (int r = 0; r != Square; ++r) {
			if (r == row) continue;

			for (int c = 0; c != Square; ++c) {
				if (c == col) continue;
				*p++ = m[r][c];
			}
		}

		RetType cofactor = determinant<RetType>(minorMat);
		cofactor = ((row + col) & 1) ? -cofactor : cofactor;
		return cofactor;
	}

	template<typename RetType, typename Type, int Row, int Col>
	auto make_adjoint(const mat_base<Type, Row, Col>& m) {
		//TODO:: error handling?
		return mat_base<RetType, Row, Col>();
	}

	template<typename RetType, typename Type, int Square>
	auto make_adjoint(const mat_base<Type, Square, Square>& m) {
		mat_base<RetType, Square, Square> ret;

		for (int r = 0; r != Square; ++r) {
			for (int c = 0; c != Square; ++c) {
				//transpose.
				ret[c][r] = cofactor<RetType>(m, r, c);
			}
		}

		return ret;
	}

	template<typename RetType, typename Type, int Row, int Col>
	auto make_inverse(const mat_base<Type, Row, Col>& m, RetType& det) {
		//TODO:: error handling?
		return mat_base<RetType, Row, Col>();
	}

	template<typename RetType, typename Type, int Row, int Col>
	auto make_inverse(const mat_base<Type, Row, Col>& m) {
		//TODO:: error handling?
		return mat_base<RetType, Row, Col>();
	}

	template<typename RetType, typename Type>
	auto make_inverse(const mat2x2_t<Type>& m, RetType& det) {
		using VecT = vec2_t<RetType>;
		using RetT = mat2x2_t<RetType>;

		det = determinant<RetType>(m);

		if (Math::is_zero(det)) {
			return RetT(); //all elements zero.
		}

		RetType invDet = Math::div<RetType>(1, det);

		return RetT{
			VecT{Math::mul<RetType>( m[1][1], invDet) , Math::mul<RetType>(-m[0][1], invDet)},
			VecT{Math::mul<RetType>(-m[1][0], invDet) , Math::mul<RetType>( m[0][0], invDet)}
		};
	}

	template<typename RetType, typename Type>
	auto make_inverse(const mat3x3_t<Type>& m, RetType& det) {
		using VecT = vec3_t<Type>;
		using RetT = mat3x3_t<RetType>;
	
		det = determinant<RetType>(m);
	
		if (Math::is_zero(det)) {
			return RetT(); //all elements zero.
		}
	
		RetType invDet = Math::div<RetType>(1, det);
	
		return RetT{
			VecT{  (Math::mul<RetType>(m[1][1], m[2][2]) - Math::mul<RetType>(m[2][1], m[1][2])) * invDet ,
			      -(Math::mul<RetType>(m[0][1], m[2][2]) - Math::mul<RetType>(m[2][1], m[0][2])) * invDet ,
			       (Math::mul<RetType>(m[0][1], m[1][2]) - Math::mul<RetType>(m[1][1], m[0][2])) * invDet },	
			VecT{ -(Math::mul<RetType>(m[1][0], m[2][2]) - Math::mul<RetType>(m[2][0], m[1][2])) * invDet ,
				   (Math::mul<RetType>(m[0][0], m[2][2]) - Math::mul<RetType>(m[2][0], m[0][2])) * invDet ,
				  -(Math::mul<RetType>(m[0][0], m[1][2]) - Math::mul<RetType>(m[1][0], m[0][2])) * invDet },	
			VecT{  (Math::mul<RetType>(m[1][0], m[2][1]) - Math::mul<RetType>(m[2][0], m[1][1])) * invDet ,
			      -(Math::mul<RetType>(m[0][0], m[2][1]) - Math::mul<RetType>(m[2][0], m[0][1])) * invDet ,
			       (Math::mul<RetType>(m[0][0], m[1][1]) - Math::mul<RetType>(m[1][0], m[0][1])) * invDet }
		};
	}

	template<typename RetType, typename Type>
	auto make_inverse(const mat4x4_t<Type>& m, RetType& det) {
		using RetT = mat4x4_t<RetType>;
		using VecT = vec4_t<RetType>;

		RetType sA = Math::mul<RetType>(m[0][0], m[1][1]) - Math::mul<RetType>(m[1][0], m[0][1]);
		RetType sB = Math::mul<RetType>(m[0][0], m[1][2]) - Math::mul<RetType>(m[1][0], m[0][2]);
		RetType sC = Math::mul<RetType>(m[0][0], m[1][3]) - Math::mul<RetType>(m[1][0], m[0][3]);
		RetType sD = Math::mul<RetType>(m[0][1], m[1][2]) - Math::mul<RetType>(m[1][1], m[0][2]);
		RetType sE = Math::mul<RetType>(m[0][1], m[1][3]) - Math::mul<RetType>(m[1][1], m[0][3]);
		RetType sF = Math::mul<RetType>(m[0][2], m[1][3]) - Math::mul<RetType>(m[1][2], m[0][3]);

		RetType cA = Math::mul<RetType>(m[2][0], m[3][1]) - Math::mul<RetType>(m[3][0], m[2][1]);
		RetType cB = Math::mul<RetType>(m[2][0], m[3][2]) - Math::mul<RetType>(m[3][0], m[2][2]);
		RetType cC = Math::mul<RetType>(m[2][0], m[3][3]) - Math::mul<RetType>(m[3][0], m[2][3]);
		RetType cD = Math::mul<RetType>(m[2][1], m[3][2]) - Math::mul<RetType>(m[3][1], m[2][2]);
		RetType cE = Math::mul<RetType>(m[2][1], m[3][3]) - Math::mul<RetType>(m[3][1], m[2][3]);
		RetType cF = Math::mul<RetType>(m[2][2], m[3][3]) - Math::mul<RetType>(m[3][2], m[2][3]);

		det = sA * cF - sB * cE + sC * cD + sD * cC - sE * cB + sF * cA;

		if (Math::is_zero(det)) {
			return RetT(); //all elements zero.
		}

		RetType invDet = Math::div<RetType>(1, det);

		return RetT{
			VecT{
			(Math::mul<RetType>( m[1][1], cF) - Math::mul<RetType>(m[1][2], cE) + Math::mul<RetType>(m[1][3], cD))* invDet ,
			(Math::mul<RetType>(-m[0][1], cF) + Math::mul<RetType>(m[0][2], cE) - Math::mul<RetType>(m[0][3], cD))* invDet ,
			(Math::mul<RetType>( m[3][1], sF) - Math::mul<RetType>(m[3][2], sE) + Math::mul<RetType>(m[3][3], sD))* invDet ,
			(Math::mul<RetType>(-m[2][1], sF) + Math::mul<RetType>(m[2][2], sE) - Math::mul<RetType>(m[2][3], sD))* invDet },

			VecT{
			(Math::mul<RetType>(-m[1][0], cF) + Math::mul<RetType>(m[1][2], cC) - Math::mul<RetType>(m[1][3], cB))* invDet ,
			(Math::mul<RetType>( m[0][0], cF) - Math::mul<RetType>(m[0][2], cC) + Math::mul<RetType>(m[0][3], cB))* invDet ,
			(Math::mul<RetType>(-m[3][0], sF) + Math::mul<RetType>(m[3][2], sC) - Math::mul<RetType>(m[3][3], sB))* invDet ,
			(Math::mul<RetType>( m[2][0], sF) - Math::mul<RetType>(m[2][2], sC) + Math::mul<RetType>(m[2][3], sB))* invDet },

			VecT{
			(Math::mul<RetType>( m[1][0], cE) - Math::mul<RetType>(m[1][1], cC) + Math::mul<RetType>(m[1][3], cA))* invDet ,
			(Math::mul<RetType>(-m[0][0], cE) + Math::mul<RetType>(m[0][1], cC) - Math::mul<RetType>(m[0][3], cA))* invDet ,
			(Math::mul<RetType>( m[3][0], sE) - Math::mul<RetType>(m[3][1], sC) + Math::mul<RetType>(m[3][3], sA))* invDet ,
			(Math::mul<RetType>(-m[2][0], sE) + Math::mul<RetType>(m[2][1], sC) - Math::mul<RetType>(m[2][3], sA))* invDet },

			VecT{
			(Math::mul<RetType>(-m[1][0], cD) + Math::mul<RetType>(m[1][1], cB) - Math::mul<RetType>(m[1][2], cA))* invDet ,
			(Math::mul<RetType>( m[0][0], cD) - Math::mul<RetType>(m[0][1], cB) + Math::mul<RetType>(m[0][2], cA))* invDet ,
			(Math::mul<RetType>(-m[3][0], sD) + Math::mul<RetType>(m[3][1], sB) - Math::mul<RetType>(m[3][2], sA))* invDet ,
			(Math::mul<RetType>( m[2][0], sD) - Math::mul<RetType>(m[2][1], sB) + Math::mul<RetType>(m[2][2], sA))* invDet}
		};
	}

	template<typename RetType, typename Type, int Square>
	auto make_inverse(const mat_base<Type, Square, Square>& m, RetType& det) {
		using RetT = mat_base<RetType, Square, Square>;
		RetT adjMat = make_adjoint<RetType>(m);
		det = RetType(0);

		for (int i = 0; i != Square; ++i) {
			det += mul<RetType>(adjMat[i][0], m[0][i]);
		}

		if (Math::is_zero(det)) {
			return RetT();// all elements zero.
		}

		return Math::mul<RetType>(adjMat, Math::div<RetType>(1, det));
	}

	template<typename RetType, typename Type>
	auto make_inverse(const mat2x2_t<Type>& m) {
		RetType unused;
		return make_inverse<RetType>(m, unused);
	}

	template<typename RetType, typename Type>
	auto make_inverse(const mat3x3_t<Type>& m) {
		RetType unused;
		return make_inverse<RetType>(m, unused);
	}

	template<typename RetType, typename Type>
	auto make_inverse(const mat4x4_t<Type>& m) {
		RetType unused;
		return make_inverse<RetType>(m, unused);
	}

	template<typename RetType, typename Type, int Square>
	auto make_inverse(const mat_base<Type, Square, Square>& m) {
		RetType unused;
		return make_inverse<RetType>(m, unused);
	}

	template<typename RetType, typename Type>
	mat2x2_t<RetType> make_transpose(const mat2x2_t<Type>& m) {
		using VecT = vec2_t<RetType>;
		using RetT = mat2x2_t<RetType>;
		return mat2x2_t<RetType>{
			m[0][0], m[1][0],
			m[0][1], m[1][1]
		};
	}

	template<typename RetType, typename Type>
	auto make_transpose(const mat3x3_t<Type>& m) {
		using VecT = vec3_t <RetType>;
		using RetT = mat3x3_t<RetType>;
		return RetT{
			VecT {m[0][0], m[1][0], m[2][0] },
			VecT {m[0][1], m[1][1], m[2][1] },
			VecT {m[0][2], m[1][2], m[2][2] }
		};
	}

	template<typename RetType, typename Type>
	auto make_transpose(const mat4x4_t<Type>& m) {
		using VecT = vec4_t <RetType>;
		using RetT = mat4x4_t<RetType>;
		return RetT{
			VecT {m[0][0], m[1][0], m[2][0], m[3][0] },
			VecT {m[0][1], m[1][1], m[2][1], m[3][1] },
			VecT {m[0][2], m[1][2], m[2][2], m[3][2] },
			VecT {m[0][3], m[1][3], m[2][3], m[3][3] }
		};
	}

	template<typename RetType, typename Type, int Row, int Col>
	auto make_transpose(const mat_base<Type, Row, Col>& m) {
		mat_base<RetType, Col, Row> ret;
		for (int r = 0; r != Row; ++r) {
			for (int c = 0; c != Col; ++c) {
				ret[c][r] = static_cast<RetType>(m[r][c]);
			}
		}
		return ret;
	}

}//namespace Matrix
}//namespace Math 
}//namespace Zee