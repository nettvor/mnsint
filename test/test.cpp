/* 
*************************************************************
Copyright © 2013 Igor Kohanovsky e-mail: Igor.Kohanovsky@gmail.com
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*************************************************************
*/
#include <iostream>
#include <iomanip>

#include <memory>
#include <random>

#include "../service/stopwatch.h"
#include "../spline/rk.h"
#include "../spd/spdchol.h"
#include "../helper/helper1.h"

#define OPENMP
//#define AMP
#ifdef PPL
#include "../helper/helper1ppl.h"
#endif
#ifdef OPENMP
#include "../helper/helper1omp.h"
#endif
#ifdef AMP
#include "../helper/helper1amp.h"
#else
#endif

using std::cin;
using std::cout;
using std::endl;

using namespace mns;

void SetPrintParams(int width, int precision, std::ios::fmtflags fmt=std::ios::fixed);

std::ostream& operator << (std::ostream& os, const mns::Status& obj)
{
   os << static_cast<std::underlying_type<mns::Status>::type>(obj);
   return os;
}

int main(int argc, char* argv[])
{
	cout << endl << "Hit <Return> key to exit..." << endl;
	cin.clear();
	cin.ignore(cin.rdbuf()->in_avail());
	cin.get();
	return 0;
}

void SetPrintParams(int width, int precision, std::ios::fmtflags fmt)
{
	cout.setf(fmt);
	cout.width(width);
	cout.precision(precision);
}


