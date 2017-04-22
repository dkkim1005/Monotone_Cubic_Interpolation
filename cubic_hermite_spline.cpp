#include <iostream>
#include <cstring>
#include <cmath>
#include <functional>
#include <cassert>
#include <vector>
#include <fstream>


class CubicHermiteSpline
{
public:
	CubicHermiteSpline(const double* xVec, const double* yVec, const double* mVec, const size_t size)
	: xVec_(nullptr), yVec_(nullptr), mVec_(nullptr), size_(0)
	{
		allocate_heap_memory_(xVec, yVec, mVec, size);
	}


	virtual ~CubicHermiteSpline()
	{
		delete_heap_memory_();
	}


	double operator()(const double x) const
	{
		const size_t idx = binary_search_(x);

		if(idx == size_ - 1) {
			return yVec_[size_-1];
		}

		const double t = (x - xVec_[idx])/(xVec_[idx+1] - xVec_[idx]);

		return interp_f_(t, idx);
	}


protected:
	double* xVec_;
	double* yVec_;
	double* mVec_;
	size_t size_;

	explicit CubicHermiteSpline()
	: xVec_(nullptr), yVec_(nullptr), mVec_(nullptr), size_(0)
	{}


	void delete_heap_memory_()
	{
		if(xVec_ != nullptr) {
			delete [] xVec_;
		}
		if(yVec_ != nullptr) {
			delete [] yVec_;
		}
		if(mVec_ != nullptr) {
			delete [] mVec_;
		}
	}


	void allocate_heap_memory_(const double* xVec, const double* yVec, const double* mVec, const size_t size)
	{
		delete_heap_memory_();

		xVec_ = new double [size];
		yVec_ = new double [size];
		mVec_ = new double [size];
		size_ = size;

		std::memcpy(xVec_, xVec,sizeof(double)*size);
		std::memcpy(yVec_, yVec,sizeof(double)*size);
		std::memcpy(mVec_, mVec,sizeof(double)*size);
	}


	size_t binary_search_(const double& x) const
	{
		assert(xVec_[0] <= x and x <= xVec_[size_-1]);
		size_t idx_l = 0, idx_r = size_ - 1, idx = size_/2;

		auto is_x_in_Boundary = [this, &x](const size_t& idx) -> bool
					{
						return this -> xVec_[idx] <= x and x < this -> xVec_[idx + 1];
					};

		while(1)
		{
			if(idx_r - idx_l == 1) 
			{
				if(is_x_in_Boundary(idx)) {
					return idx;
				}
				else {
					return idx + 1;
				}
			}

			if(is_x_in_Boundary(idx)) {
				return idx;
			}
			else if(xVec_[idx+1] <= x)
			{
				idx_l = idx;
				idx = (idx_r - idx_l)/2 + idx_l;
			}
			else
			{
				idx_r = idx;
				idx = (idx_r - idx_l)/2 + idx_l;
			}
		}
	}


	double interp_f_(const double& t, const size_t& idx) const
	{
		return (2*std::pow(t,3) - 3*std::pow(t,2) + 1)*yVec_[idx] +
		(std::pow(t,3) - 2*std::pow(t,2) + t)*(xVec_[idx+1] - xVec_[idx])*mVec_[idx] +
		(-2*std::pow(t,3) + 3*std::pow(t,2))*yVec_[idx+1] +
		(std::pow(t,3) - std::pow(t,2))*(xVec_[idx+1] - xVec_[idx])*mVec_[idx+1];
	}
};


class MonotoneCubicInterpolation : public CubicHermiteSpline
{
public:
	MonotoneCubicInterpolation(const double* xVec, const double* yVec, const size_t size)
	: CubicHermiteSpline()
	{
		std::vector<double> delta(size, 0);
		std::vector<double> mVec(size, 0);

		for(int i=0; i<size-1; ++i) {
			delta[i] = (yVec[i+1] - yVec[i])/(xVec[i+1] - xVec[i]);
		}

		for(int i=1; i<size-1; ++i) {
			mVec[i] = (delta[i-1] + delta[i])/2.;
		}
	
		mVec[0] = delta[0]; mVec[size-1] = delta[size-2];

		for(int i=0; i<size-1; ++i)
		{
			if(std::abs(delta[i]) < 1e-30) {
				mVec[i] = mVec[i+1] = 0.;
			}
		}

		allocate_heap_memory_(xVec, yVec, &mVec[0], size);
	}

	virtual ~MonotoneCubicInterpolation() {}
};



int main(int argc, char* argv[])
{
	std::vector<double> x = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
	std::vector<double> y(11, 0);


	for(int i=0; i<x.size(); ++i) {
		y[i] = std::sin(x[i]);
	}

	MonotoneCubicInterpolation tester(&x[0], &y[0], 11);

	std::ofstream file("tester.out");

	for(int i=0; i<100; ++i) 
	{
		double x = i/100.;
		file<<x<<" "<<tester(x)<<"\n";
	}

	return 0;
}
