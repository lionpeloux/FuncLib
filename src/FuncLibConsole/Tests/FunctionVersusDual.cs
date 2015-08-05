// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

using FuncLib;
using FuncLib.DualNumbers;
using FuncLib.Functions;
using FuncLib.Mathematics;

namespace FuncLibConsole.Tests
{
	public static class FunctionVersusDual
	{
		// FIXME Formalize to real unit tests

		public static void TestAll()
		{
			Test1();
		}

		public static void Test1()
		{
			// The forward-mode (DualNumber) and reverse-mode (Function) are completely different implementations, so comparing them is probably a good test.

			Variable x_ = new Variable();
			Variable y_ = new Variable();
			DualNumber x = new DualNumber(1.2345, Vector.Basis(2, 0));
			DualNumber y = new DualNumber(5.4321, Vector.Basis(2, 1));

			DualNumber f = DualNumber.Sin(DualNumber.Cos(x * x * y)) / (DualNumber.Sin(DualNumber.Exp(x * x)) * DualNumber.Cos(y * y));
			Function f_ = Function.Sin(Function.Cos(x_ * x_ * y_)) / (Function.Sin(Function.Exp(x_ * x_)) * Function.Cos(y_ * y_));
			f *= DualNumber.Cos(f);
			f_ *= Function.Cos(f_);
			f /= DualNumber.Cos(x * x * x) * DualNumber.Sin(x * DualNumber.Exp(x / DualNumber.Sqr(DualNumber.Exp(y))));
			f_ /= Function.Cos(x_ * x_ * x_) * Function.Sin(x_ * Function.Exp(x_ / Function.Sqr(Function.Exp(y_))));

			f = DualNumber.Pow(f, DualNumber.Cos(x * x * y) / DualNumber.Exp(x * y));
			f_ = Function.Pow(f_, Function.Cos(x_ * x_ * y_) / Function.Exp(x_ * y_));

			Console.WriteLine(f_.Value(x_ | x.Value, y_ | y.Value));
			Console.WriteLine(f.Value);
			Console.WriteLine();

			Console.WriteLine(f_.Derivative(x_).Value(x_ | x.Value, y_ | y.Value));
			Console.WriteLine(f.Gradient[0]);
			Console.WriteLine();

			Console.WriteLine(f_.Derivative(y_).Value(x_ | x.Value, y_ | y.Value));
			Console.WriteLine(f.Gradient[1]);
			Console.WriteLine();

			Console.WriteLine(f_.Derivative(x_, x_).Value(x_ | x.Value, y_ | y.Value));
			Console.WriteLine(f.Hessian[0, 0]);
			Console.WriteLine();

			Console.WriteLine(f_.Derivative(x_, y_).Value(x_ | x.Value, y_ | y.Value));
			Console.WriteLine(f.Hessian[0, 1]);
			Console.WriteLine();

			Console.WriteLine(f_.Derivative(y_, y_).Value(x_ | x.Value, y_ | y.Value));
			Console.WriteLine(f.Hessian[1, 1]);
			Console.WriteLine();

			//Console.WriteLine(f_.Derivative(x_, x_, x_).Value(x_ | x.Value, y_ | y.Value));
			//Console.WriteLine(f.Third[0, 0, 0]);
			//Console.WriteLine();

			//Console.WriteLine(f_.Derivative(x_, x_, y_).Value(x_ | x.Value, y_ | y.Value));
			//Console.WriteLine(f.Third[0, 0, 1]);
			//Console.WriteLine();

			//Console.WriteLine(f_.Derivative(x_, y_, y_).Value(x_ | x.Value, y_ | y.Value));
			//Console.WriteLine(f.Third[0, 1, 1]);
			//Console.WriteLine();

			//Console.WriteLine(f_.Derivative(y_, y_, y_).Value(x_ | x.Value, y_ | y.Value));
			//Console.WriteLine(f.Third[1, 1, 1]);
			//Console.WriteLine();

			//throw new Exception();
		}
	}
}
