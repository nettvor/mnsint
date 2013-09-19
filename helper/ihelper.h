/* 
*************************************************************
Copyright © 2013 Igor Kohanovsky e-mail: Igor.Kohanovsky@gmail.com
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*************************************************************
*/
#pragma once
#ifndef __IHelper_H__
#define __IHelper_H__

#include <cmath>
#include "../common/defs.h"

namespace mns 
{
	inline Size_T GetIndex(int i, int j) { return (i >= j) ? j + ((Size_T)i) * (i + 1) / 2 : i + ((Size_T)j) * (j + 1) / 2; };

	template <typename T>
	class IHelper
	{
	// Defines interface for common vector/matrix operations
	public:
		typedef typename Defs<T>::VectorT VectorT;
		typedef typename Defs<T>::SpdMatrixT SpdMatrixT;

		inline T PI() const { return 3.141592653589793238L; }
		inline T SQRTPI() const { return std::sqrt(PI()); }

		VectorT GetResidual(int n, const SpdMatrixT& a, const VectorT& x, const VectorT& b) const { return GetResidualImpl(n, a, x, b); };
		T		GetVectorNorm2(int n, const VectorT& v) const { return GetVectorNorm2Impl(n, v); }

		T		GetGamma2(int n) const { return GetGamma2Impl(n); }

		virtual ~IHelper() {};
	protected:
		virtual VectorT GetResidualImpl(int n, const SpdMatrixT& a, const VectorT& x, const VectorT& b) const abstract; 
		virtual T GetVectorNorm2Impl(int n, const VectorT& v) const abstract;

		virtual T GetGamma2Impl(int n) const abstract;

		IHelper() {};
	private:
    	IHelper(const IHelper&);
		IHelper& operator =(const IHelper&);
		IHelper& operator =(IHelper&&);
	};

} // end of mns namespace

#endif // __IHelper_H__
