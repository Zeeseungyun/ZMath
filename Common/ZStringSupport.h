#pragma once
#include <string>
#include "ZTypePredefines.h"
namespace Zee {
namespace Math {
namespace String {
	template<typename CharType> //default char.
	struct ToString {
		using type = std::string;
		template<typename T>
		static std::string to_string(const T& v) {
			return std::to_string(v);
		}
	};

	template<> 
	struct ToString<wchar_t> {
		using type = std::wstring;
		template<typename T>
		static std::wstring to_string(const T& v) {
			return std::to_wstring(v);
		}
	};
}//namespace String 
	template<typename CharType, typename QuatType>
	std::basic_string<CharType> to_string(const quat_base<QuatType>& q) {
		std::basic_string<CharType> str;
		str += "quat[" 
			+ String::ToString<CharType>::to_string(q[0]) + ", "
			+ String::ToString<CharType>::to_string(q[1]) + ", "
			+ String::ToString<CharType>::to_string(q[2]) + ", "
			+ String::ToString<CharType>::to_string(q[3]) 
			+ "]";
		return str;
	}

	template<typename CharType, typename VecType, int Dim>
	std::basic_string<CharType> to_string(const vec_base<VecType, Dim>& v) {
		std::basic_string<CharType> str;
		str += "vec[";
		for (int d = 0; d != Dim - 1; ++d) {
			str += String::ToString<CharType>::to_string(v[d]) + ", ";
		}
		str += String::ToString<CharType>::to_string(v[Dim - 1]);
		str += "]";
		return str;
	}

	template<typename CharType, typename MatType, int Row, int Col>
	std::basic_string<CharType> to_string(const mat_base<MatType, Row, Col>& m) {
		std::basic_string<CharType> str;
		str += "mat{\n";
		for (int r = 0; r != Row - 1; ++r) {
			str += to_string<CharType>(m[r]) + ", \n";
		}
		str += to_string<CharType>(m[Row - 1]);
		str += "\n}";
		return str;
	}

}//namespace Math 
}//namespace Zee 

