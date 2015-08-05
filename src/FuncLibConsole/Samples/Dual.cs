// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;

using FuncLib;
using FuncLib.DualNumbers;
using FuncLib.Functions;
using FuncLib.Mathematics;
using FuncLib.Optimization.Ipopt;

namespace FuncLibConsole.Samples
{
	public class Dual
	{
		public static void Sample1()
		{
			// Define DualNumber as variables. We specify our knowledge about the partial derivatives of our variables.
			DualNumber x = new DualNumber(2.0, new Vector(new double[] { 1.0, 0.0 }));
			DualNumber y = new DualNumber(3.0, new Vector(new double[] { 0.0, 1.0 }));

			// Computer a new value through some aritmetics. This is our function.
			DualNumber f = x * DualNumber.Exp(x * y);

			// The DualNumber object keeps track of first and second derivatives at the particular point. Show the second order partial derivative.
			Console.WriteLine(f.Hessian[0, 1]);
		}

		public static void Sample2()
		{
			// It's possible to define a traditional FuncLib Function through a DualNumber.

			// Traditional variables (i.e. not DualNumber).
			Variable x = new Variable();
			Variable y = new Variable();

			// Create the Function using DualNumberFunction and specify that it depends on those variables.
			Function f = DualNumberFunction.Create(
				delegate(IDualNumberTransform transform)
				{
					// When this delegate is called a DualNumber is assigned to represent each Variable.
					DualNumber x0 = transform[x];
					DualNumber y0 = transform[y];

					// Compute a new DualNumber just as in the previous sample.
					return x0 * DualNumber.Exp(x0 * y0);
				}, x, y);

			// Now the second order partial derivative may be computer like any other Function. Behind the scenes it just keeps track of the right entry of the Hessian.
			Function g = f.Derivative(x).Derivative(y);

			// Evaluate function for (x,y) = (2,3).
			Console.WriteLine(g.Value(x | 2.0, y | 3.0));
		}

		public static void Sample3()
		{
			// This sample shows how to use DualNumber with one of the optimizers.

			// Traditional variables (i.e. not DualNumber).
			Variable x = new Variable();
			Variable y = new Variable();

			// Create the Function using DualNumberFunction.
			Function f = DualNumberFunction.Create(
				delegate(IDualNumberTransform transform)
				{
					// When this delegate is called a DualNumber is assigned to represent each Variable.
					DualNumber x0 = transform[x];
					DualNumber y0 = transform[y];

					// Compute the actual function using DualNumber. In this case the Rosenbrock function.
					return DualNumber.Sqr(1.0 - x0) + 100.0 * DualNumber.Sqr(y0 - DualNumber.Sqr(x0));
				}, x, y);

			// Now just the usual procedure to start the optimizer.
			IpoptOptimizer o = new IpoptOptimizer();
			o.Variables.Add(x, y);
			o.ObjectiveFunction = f;
			o.PrintLevel = 5;
			Console.WriteLine(o.Run(x | 0.5, y | 0.3).OptimalPoint);
		}
	}
}
