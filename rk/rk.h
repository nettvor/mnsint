/* 
*************************************************************
Copyright © 2013 Igor Kohanovsky e-mail: Igor.Kohanovsky@gmail.com
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*************************************************************
*/

#pragma once
#ifndef __RK_H__
#define __RK_H__

#include <cmath>
#include <numeric>
#include "irk.h"

namespace mns 
{
	template <typename T>
	class RK  : public IRK<T> 
	{
	// Computes a Reproducing Kernel
	public:
		RK(int r, T eps) : r_(r), eps_(eps) {};

		//T BFun(T r) const;
		//private
		T ExpBySquaring(T x, int n) const;
		long Fact(int n) const;
		virtual ~RK() {};
	protected:

		RK() {};
	private:
		//T GetPolyCoefficients(int r, T eps, VectorT& a) const;

		

    	RK(const RK&);
		RK& operator =(const RK&);
		RK& operator =(RK&&);

		int r_;
		T eps_;
	};

	template<typename T> 
	long RK<T>::Fact(int n) const
	{
		return (n == 1)? 1: n * fact(n - 1);
	}

	template<typename T> 
	T RK<T>::ExpBySquaring(T x, int n) const
	// Calculates integer power
	{
		if ( n < 0 )
		{
			return ExpBySquaring(T(1.0)/x, -n);
		}
		if ( n == 0 )
		{
			return T(1.0);
		}
		if ( n == 1 )
		{
			return x;
		}
		if ( n%2 == 0 )
		{
			return ExpBySquaring(x*x, n/2);
		}

		return x * ExpBySquaring(x*x, (n-1)/2);
	}

	//template<typename T> 
	//T RK<T>::GetPolyCoefficients(int r, T eps, VectorT& a) const
	//// Calculates coefficients of the Reproducing Kernel polynomial part
	//{
	//	T f = T(1.0);

	//	//if( k2 > 1 )
	//	//{
	//	//	for( i = k1; i <= k2; ++i )
	//	//		f *= i;
	//	//}
	//	return f;
	//}


//
//void nsfa(int r, vec& a)
///* forming a */
//{
// int k, i;
// real s1, s2;
//
// a[r] = 1.;
// if (r == 0) return;
// a[r-1] = 1.;
// if (r == 1) return;
//
// for(k=0; k<=r-2; k++) {
//   s1 = 1.;
//   for(i=1; i<=r-k; i++) s1 *= ((real)2./i);
//   s2 = 1.;
//   for(i=1; i<=r; i++) s2 *= ((real)(k+i)/(r+i));
//   a[k] = s1 * s2;
// }
//  return;
//}
//
//real nsp(int r, real t, vec& a)
///* rk gorner*/
//{
// int i;
// real s;
//
// s = a[0];
// for(i=1; i<=r; i++) s = s*t + a[i];
// return s;
//}

} // end of MNS namespace

#endif /* RK */


