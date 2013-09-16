/* 
*************************************************************
Copyright © 2013 Igor Kohanovsky e-mail: Igor.Kohanovsky@gmail.com
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*************************************************************
*/
#pragma once
#ifndef __HELPER1OMP_H__
#define __HELPER1OMP_H__

#include <cmath>
#include <functional>
#include <numeric>
#include <omp.h>
#include "ihelper.h"

namespace mns 
{
	template <typename T>
	class HelperOmp final : public IHelper<T>
	{
	// Implements common vector/matrix operations
	public:
		HelperOmp() {};

		int  GetNumThreads() const;
		void SetNumThreads(int n) const;
		int  GetNumProcs() const;

	private:
		virtual VectorT GetResidualImpl(int n, const SpdMatrixT& a, const VectorT& x, const VectorT& b) const override;
		virtual T GetVectorNorm2Impl(int n, const VectorT& v) const override;

		HelperOmp(const HelperOmp&);
		HelperOmp& operator =(const HelperOmp&);
		HelperOmp& operator =(HelperOmp&&);
	};

	template<typename T> 
	void HelperOmp<T>::SetNumThreads(int n) const
	{
		omp_set_num_threads(n);
	}

	template<typename T> 
	int HelperOmp<T>::GetNumProcs() const
	{
		return omp_get_num_procs();
	}

	template<typename T> 
	int HelperOmp<T>::GetNumThreads() const
	{
		return omp_get_num_threads();
	}

	template<typename T> 
	T HelperOmp<T>::GetVectorNorm2Impl(int n, const VectorT& v) const
	{
		T s = T(0.0);
		int i;

		#pragma omp parallel for reduction(+: s)
		for ( i = 0; i < n; ++i ) 
		{
			s += v[i] * v[i];
		}

		return std::sqrt(s);
	}

/*
#pragma omp parallel  
{ 
	#pragma omp for  
	for (i=0; i<x; i++) 
		fn1(); 

	#pragma omp for  
	for (i=0; i<y; i++)  
		fn2();
} 
*/
	template<typename T> 
	typename HelperOmp<T>::VectorT HelperOmp<T>::GetResidualImpl(int n, const SpdMatrixT& a, const VectorT& x, const VectorT& b) const
	{
		VectorT r(n);
		int j;

		#pragma omp parallel for
		for( int i = 0; i < n; ++i )
		{
			r[i] = T(0.0);
			for( j = 0; j < i; ++j )
			{
				r[i] += a[j + ((Size_T)i) * (i + 1) / 2] * x[j];
			}
			for( j = i; j < n; ++j )
			{
				r[i] += a[i + ((Size_T)j) * (j + 1) / 2] * x[j];
			}
			r[i] = b[i] - r[i];
		}
		return r;
	}


} // end of MNS namespace

#endif // __HELPER1OMP_H__
