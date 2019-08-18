#pragma once
#include "../Common/ZTypeTraits.h"
#include "../Common/ZCommonFunctions.h"

namespace Zee {
namespace Math {
	template<typename ScalarType>
	struct complex_base {
		using scalar_type = ScalarType;

		struct constant;

		constexpr complex_base() : data{ scalar_type() } {
			static_assert(is_arithmetic<scalar_type>::value, "scalar type must arithemtic type.");
		}

		template<typename ...Args>
		constexpr complex_base(scalar_type a, Args ...args) : data{ static_cast<scalar_type>(a), static_cast<scalar_type>(args)... } {
		}

		union {
			scalar_type data[2];
			struct { scalar_type x, y; };
		};
	};

	template<typename ScalarType>
	struct complex_base<ScalarType>::constant {
		static constexpr complex_base zero = complex_base();
		static constexpr complex_base identity = complex_base(0, 1);
	};

}//namespace Math 
}//namespace Zee