#include <iostream>
#include "ZMath.h"
#include "Common/ZStringSupport.h"

int main() {
	using namespace Zee;
	using std::cout;
	using std::endl;
	Zee::Math::vec_base<float, 2> vf2;
	Zee::Math::vec_base<double, 2> vd2;
	Zee::Math::promotion_t< Zee::Math::vec_base<double, 2>, Zee::Math::vec_base<float, 2>> dd;
	Zee::Math::promotion_t<float, double> ffasdf;

	double() + int();


	vec2 v = { 5, -5 };
	cout << " v : " << Math::to_string<char>(v) << endl;
	v = Math::clamp<float>(v, 0, 3);
	cout << " after clamp : " << Math::to_string<char>(v) << endl;

	v = Math::vec2_t<int>{ 10, 10 };
	cout << " other scalar ctor : " <<Math::to_string<char>(v) << endl;
	cout << " conversion : " << Math::to_string<char>(static_cast<Math::vec2_t<int>>(v)) <<endl;

	constexpr mat2 identity = { vec2::constant::basis[0], vec2::constant::basis[1] };
	cout << "constexpr identity : " << Math::to_string<char>(identity) << endl;

	mat2 m = { identity[0] * 2.0f, identity[1] * 0.5f };
	cout << "m : " <<Math::to_string<char>(m) << endl;
	m.determinant();
	m.get_inverse();
	
	mat3 m3;
	m3[0] = vec3::constant::unit_x;
	m3[1] = vec3::constant::unit_x;
	m3[2] = vec3::constant::unit_x;
	auto trans = m3.get_transpose();
	cout << "trans : " << Math::to_string<char>(trans) << endl;

	vec2 v2 = v * m;
	cout << "v2 (v * m) : " <<Math::to_string<char>(v2) << endl;
	vec3 v3 = vec3::constant::unit_y;
	vec3 v4 = vec3::constant::unit_z;
	cout << " v3 cross v4 : " << 
		Math::to_string<char>(vec3::personal::cross(v3, v4)) << endl;

	system("pause");
	return 0;
}