#include <cmath>
#include <vector>
#include <memory>
#include <cassert>
#include <fstream>
#include <cstring>
#include <iostream>


template <typename T>
class CubicHermiteSpline {
public:
	CubicHermiteSpline(const T * x_ptr,
                       const T * y_ptr,
                       const T * m_ptr,
                       const size_t size):
    size_(size) {
        x_.resize(size_);
        y_.resize(size_);
        m_.resize(size_);
		std::memcpy(x_.data(), x_ptr, sizeof(T) * size_);
		std::memcpy(y_.data(), y_ptr, sizeof(T) * size_);
		std::memcpy(m_.data(), m_ptr, sizeof(T) * size_);
    }

    T get_interpolated_value(const T x) const {
		const size_t idx = binary_search_(x);
		if(idx == size_ - 1)
		    return y_[size_-1];
		const T t = (x - x_[idx]) / (x_[idx+1] - x_[idx]);
		return interp_func_(t, idx);
	}

private:
    bool is_x_in_boundary(const size_t idx, const T x) const {
	    return (x_[idx] <= x) && (x < x_[idx+1]);
	}

	size_t binary_search_(const T x) const {
	    assert((x_[0] <= x) && (x <= x_[size_-1]));
		size_t idx_l = 0, idx_r = size_ - 1, idx = size_ / 2;
		while (1) {
			if(idx_r - idx_l == 1)  {
				if(is_x_in_boundary(idx, x))
					return idx;
				else
					return (idx + 1);
			}
			if(is_x_in_boundary(idx, x))
				return idx;
			else if(x_[idx+1] <= x) {
				idx_l = idx;
				idx = (idx_r - idx_l) / 2 + idx_l;
			}
			else {
				idx_r = idx;
				idx = (idx_r - idx_l) / 2 + idx_l;
			}
		}
	}

	T interp_func_(const T t, const size_t idx) const {
	    return (2 * std::pow(t, 3) - 3 * std::pow(t, 2) + 1) * y_[idx] +
		       (std::pow(t, 3) - 2 * std::pow(t, 2) + t) * (x_[idx+1] - x_[idx]) * m_[idx] +
		       (-2 * std::pow(t, 3) + 3 * std::pow(t, 2)) * y_[idx+1] +
		       (std::pow(t, 3) - std::pow(t, 2))*(x_[idx+1] - x_[idx]) * m_[idx+1];
	}

    std::vector<T> x_, y_, m_;
	size_t size_;
};


template <typename T>
class MonotoneCubicInterpolation {
public:
	MonotoneCubicInterpolation(const T * x_ptr, const T * y_ptr, const size_t size) {
		std::vector<T> delta(size, 0);
		std::vector<T> m(size, 0);
		for(int i=0; i<size-1; ++i) delta[i] = (y_ptr[i+1] - y_ptr[i]) / (x_ptr[i+1] - x_ptr[i]);
		for(int i=1; i<size-1; ++i) m[i] = (delta[i-1] + delta[i]) / 2;
		m[0] = delta[0];
        m[size-1] = delta[size-2];
		for(int i=0; i<(size-1); ++i) {
			if(std::abs(delta[i]) < keps)
				m[i] = m[i+1] = 0;
		}
        spliner_ptr_ = std::make_unique<CubicHermiteSpline<T>>(x_ptr, y_ptr, m.data(), size);
	}

    T operator()(const T x) const {
        return spliner_ptr_ -> get_interpolated_value(x);
    }

private:
    static constexpr const T keps = 1e-10;
    std::unique_ptr<CubicHermiteSpline<T>> spliner_ptr_;
};


int main(int argc, char* argv[])
{
    using T = double;
	std::vector<T> x = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
	std::vector<T> y(x.size(), 0);
	for(int i=0; i<x.size(); ++i) y[i] = std::sin(x[i]);
	MonotoneCubicInterpolation<T> tester(x.data(), y.data(), x.size());
	std::ofstream file("tester.out");
    const int nPoints = 100;
	for(int i=0; i<nPoints; ++i) {
		auto x = i / static_cast<T>(nPoints);
		file << x << " " << tester(x) << "\n";
	}
	return 0;
}
