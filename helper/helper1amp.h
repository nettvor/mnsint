/* 
*************************************************************
Copyright © 2013 Igor Kohanovsky e-mail: Igor.Kohanovsky@gmail.com
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*************************************************************
*/
#pragma once
#ifndef __HELPER1AMP_H__
#define __HELPER1AMP_H__

#include <cmath>
#include <functional>
#include <numeric>
#include <amp.h>
#include "ihelper.h"

using namespace concurrency;

namespace mns 
{
	template <typename T>
	class Helper1 final : public IHelper<T>
	{
	// Implements common vector/matrix operations
	public:
		Helper1() {};
		bool HasAccelerator() const { return accelerator::get_all().size() > 0; };
		bool HasHWAccelerator() const;

	private:
		virtual VectorT GetResidualImpl(int n, const SpdMatrixT& a, const VectorT& x, const VectorT& b) const override;
		virtual T GetVectorNorm2Impl(int n, const VectorT& v) const override;

		Helper1(const Helper1&);
		Helper1& operator =(const Helper1&);
		Helper1& operator =(Helper1&&);
	};


	template<typename T> 
	bool Helper1<T>::HasHWAccelerator() const
	{
		bool res = false;

		std::vector<accelerator> accls = accelerator::get_all();
		accls.erase(std::remove_if(accls.begin(), accls.end(), [](accelerator& a) { 
						return (a.device_path == accelerator::cpu_accelerator) || (a.device_path == accelerator::direct3d_ref); 
					}), accls.end());

		if ( !accls.empty() )
		{
			res = true;
		}
		return res;
	}

// bool can_use_doubles = accelerator().supports_double_precision;
// bool can_use_doubles_with_limits = accelerator().supports_limited_double_precision;

	template<typename T> 
	T Helper1<T>::GetVectorNorm2Impl(int n, const VectorT& v) const
	{
		array_view<const T, 1> a(n, v);
		T s = 0.0;

		return std::sqrt(s);
	}

	template<typename T> 
	typename Helper1<T>::VectorT Helper1<T>::GetResidualImpl(int n, const SpdMatrixT& a, const VectorT& x, const VectorT& b) const
	{
		VectorT r(n);
		
		return r;
	}


} // end of MNS namespace

#endif // __HELPER1AMP_H__
