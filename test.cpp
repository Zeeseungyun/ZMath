#include <iostream>
#include "ZMath.h"
#include "Common/ZStringSupport.h"

int main() {
	using namespace Zee;
	using std::cout;
	using std::endl;

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

	vec2 v2 = v * m;
	cout << "v2 (v * m) : " <<Math::to_string<char>(v2) << endl;
	vec3 v3 = vec3::constant::unit_y;
	vec3 v4 = vec3::constant::unit_z;
	cout << " v3 cross v4 : " << 
		Math::to_string<char>(vec3::personal::cross(v3, v4)) << endl;

	system("pause");
	return 0;
}