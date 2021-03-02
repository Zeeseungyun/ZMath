#pragma once
#include <bitset>

namespace Zee {
namespace Math {

	template<int Row, int Col >
	struct bools {
		static constexpr int row = Row; static constexpr int col = Col;

		const bools<col, 0>& operator[](int idx) const { return data[idx]; }
		bools<col, 0>& operator[](int idx) { return data[idx]; }

	public:
		bools<col, 0> data[row];
	};

	template<int Dimension>
	struct bools<Dimension, 0> : public std::bitset<Dimension> {
		static constexpr int dim = Dimension;
		static const bools all_true; 
		static const bools all_false; 

		bools(const std::bitset<dim>& r) : std::bitset<dim>(r) { }
		bools(std::bitset<dim>&& r) : std::bitset<dim>(std::move(r)) { }
		bools& operator=(const std::bitset<dim>& r) { std::bitset<dim>::operator=(r); return *this; }

		bools() : bools(all_false) {}
		bools(const bools&) = default;
		bools(bools&&) = default;
		bools& operator=(const bools&) = default;
		bools& operator=(bools&&) = default;

		operator std::bitset<dim>() { return static_cast<std::bitset<dim>>(*this); }
		operator const std::bitset<dim>() const { return static_cast<std::bitset<dim>>(*this); }
	};

	template<int Dim>
	const bools<Dim, 0> bools<Dim, 0>::all_true = std::bitset<dim>().flip();

	template<int Dim>
	const bools<Dim, 0> bools<Dim, 0>::all_false = std::bitset<dim>();

	template<int Row, int Col>
	bool operator==(const bools<Row, Col>& lhs, const bools<Row, Col>& rhs) {
		for (int r = 0; r != Row; ++r) {
			if (lhs[r] != rhs[r]) return false;
		}
		return true;
	}

	template<int Row, int Col>
	bool operator!=(const bools<Row, Col>& lhs, const bools<Row, Col>& rhs) {
		for (int r = 0; r != Row; ++r) {
			if (lhs[r] == rhs[r]) return false;
		}
		return true;
	}

	template<int Dim>
	bool operator==(const bools<Dim, 0>& lhs, const bools<Dim, 0>& rhs) {
		return static_cast<std::bitset<Dim>>(lhs) == static_cast<std::bitset<Dim>>(rhs);
	}

	template<int Dim>
	bool operator!=(const bools<Dim, 0>& lhs, const bools<Dim, 0>& rhs) {
		return static_cast<std::bitset<Dim>>(lhs) != static_cast<std::bitset<Dim>>(rhs);
	}

}//namespace Math 
}//namespace Zee