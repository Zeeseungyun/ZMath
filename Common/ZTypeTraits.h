#pragma once
#include <type_traits>
#include "ZTypePredefines.h"
namespace Zee {
namespace Math {
	using sint8  = signed char;
	using sint16 = signed short;
	using sint32 = signed int;
	using sint64 = signed long long int;

	using uint8  = unsigned char;
	using uint16 = unsigned short;
	using uint32 = unsigned int;
	using uint64 = unsigned long long int;

	using std::is_arithmetic;
	using std::is_floating_point;
	using std::is_same;
	using std::remove_cv;
	using std::remove_reference;
	using std::remove_cv_t;
	using std::remove_reference_t;

namespace Details {
	template<typename T, int D>			constexpr bool is_vec_support(const vec_base<T, D>& unused) { return true; }
	template<typename T>				constexpr bool is_vec_support(const T& unused) { return false; }

	template<typename T, int R , int C> constexpr bool is_mat_support(const mat_base<T, R, C>& unused) { return true; }
	template<typename T>				constexpr bool is_mat_support(const T& unused) { return false; }

	template<typename T>				constexpr bool is_quat_support(const quat_base<T>& unused) { return true; }
	template<typename T>				constexpr bool is_quat_support(const T& unused) { return false; }

	template<typename T>				constexpr bool is_complex_support(const complex_base<T>& unused) { return true; }
	template<typename T>				constexpr bool is_complex_support(const T& unused) { return false; }

	template<typename T, bool is_arithmetic = is_arithmetic<T>::value>
	struct inner_type {
	};

	template<typename T> struct inner_type <T, true> { using type = T; };
	template<typename T> struct inner_type <T, false> { using type = typename T::scalar_type; };

}//namespace Details 

	template<typename T> struct remove_cv_reference { using type = T; };
	template<typename T> struct remove_cv_reference<T&> { using type = T; };
	template<typename T> struct remove_cv_reference<T&&> { using type = T; };
	template<typename T> struct remove_cv_reference<volatile T> { using type = T; };
	template<typename T> struct remove_cv_reference<const T>    { using type = T; };
	template<typename T> struct remove_cv_reference<const T&>   { using type = T; };
	template<typename T> struct remove_cv_reference<const T&&>   { using type = T; };
	template<typename T> struct remove_cv_reference<volatile T&> { using type = T; };
	template<typename T> struct remove_cv_reference<volatile T&&> { using type = T; };
	template<typename T> struct remove_cv_reference<const volatile T> { using type = T; };
	template<typename T> struct remove_cv_reference<const volatile T&> { using type = T; };
	template<typename T> struct remove_cv_reference<const volatile T&&> { using type = T; };

	template<typename T>
	using remove_cv_reference_t = typename remove_cv_reference<T>::type;

	template<typename T> struct is_vector		: std::bool_constant<Details::is_vec_support(T())> {};
	template<typename T> struct is_matrix		: std::bool_constant<Details::is_mat_support(T())> {};
	template<typename T> struct is_quaternion	: std::bool_constant<Details::is_quat_support(T())> {};
	template<typename T> struct is_complex		: std::bool_constant<Details::is_complex_support(T())> {};	
	
	template<typename T, typename ...Args> struct promotion_type {
		using type = decltype(typename Details::inner_type<remove_cv_reference_t<T>>::type() + typename promotion_type<Args...>::type());
	};
	
	template<typename T>
	struct promotion_type<T> { using type = typename Details::inner_type<remove_cv_reference_t<T>>::type; };

	template<typename ...Args>
	using promotion_t = typename promotion_type<Args...>::type;

}//namespace Math
}//namespace Zee