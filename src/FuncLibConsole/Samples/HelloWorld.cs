// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Diagnostics;

using FuncLib.Functions;
using FuncLib.Functions.Compilation;

namespace FuncLibConsole.Samples
{
	public class HelloWorld
	{
		public static void Sample1()
		{
			// Variables.
			Variable x = new Variable();
			Variable y = new Variable();

			// Define function.
			Function f = x * Function.Exp(x * y);

			// Compute second order partial derivative.
			Function g = f.Derivative(x).Derivative(y);

			// Evaluate function for (x, y) = (2, 3).
			Console.WriteLine(g.Value(x | 2.0, y | 3.0));

			// Or without operator overloading.
			Console.WriteLine(g.Value(new Point(new VariableAssignment(x, 2.0), new VariableAssignment(y, 3.0))));

			// Compile to IL code using the variables given.
			CompiledFunction h = Compiler.Compile(g, x, y);

			// Evaluate again. Now the order of the variables is fixed like an ordinary C# method.
			Console.WriteLine(h.Value(2.0, 3.0));

			// A matrix with mixed functions and constants.
			FunctionMatrix a = new FunctionMatrix(new Function[,] { { f, g }, { 1.0, f * g } });

			// Compute inverse and evaluate derivative of first entry.
			Console.WriteLine(FunctionMatrix.Inverse(a)[0, 0].Derivative(x).Value(x | 2.0, y | 3.0));

			// Use of PartialValue.
			Console.WriteLine(g.Value(x | 2.0, y | 3.0));
			Console.WriteLine(g.PartialValue(x | 2.0).Value(y | 3.0));
			Console.WriteLine(g.PartialValue(x | 2.0, y | 3.0).Value());
			Console.WriteLine(g.PartialValue(x | 2.0).Derivative(x).Value(y | 3.0));
		}

		public static void Sample2()
		{
			// Variables.
			Variable x = new Variable();
			Variable y = new Variable();

			// Define function a complicated function.
			Function f = Function.Exp(x * y) + x + 0.5;
			f *= Function.Cos(f);
			f *= Function.Cos(f);
			f *= Function.Cos(f);
			f *= Function.Cos(f);
			f = f.Derivative(x);

            // Evaluate the usual way.
            Console.WriteLine(f.Value(x | 2.0, y | 3.0));

            // Compile the function to IL code.
			CompiledFunction g = Compiler.Compile(f, x, y);

			// Compile the function (including first and second derivatives) to IL code.
			CompiledFunction g2 = Compiler.Compile(f, new Variable[] { x, y }, 1); /* 2 */

			//g = g2;

            // Evaluate.
            Console.WriteLine(g.Value(2.0, 3.0));

            // Or we may still use the Function syntax.
            Console.WriteLine(g.Value(x | 2.0, y | 3.0));

            //return;

            // Measure performance difference.
            int n = 500000;

			Stopwatch sw1 = Stopwatch.StartNew();
			for (int i = 0; i < n; i++)
			{
				f.Value(x | 2.0, y | 3.0);
			}
			Console.WriteLine(sw1.Elapsed);

            Stopwatch sw2 = Stopwatch.StartNew();
			for (int i = 0; i < n; i++)
			{
				g.Value(2.0, 3.0);
				//Evaluation e = new Evaluation(x | 2.0, y | 3.0);
				//g.Value(e);
				//g2.Derivative(x).Value(e);
				//g2.Derivative(y).Value(e);
			}
            Console.WriteLine(sw2.Elapsed);
		}

		public static void Sample3()
		{
			// Variables.
			Variable x = new Variable();

			Function f = Function.Positive(Function.Cos(x));
			//f = f.Derivative(x);
			//CompiledFirstDerivativeFunction g = new CompiledFirstDerivativeFunction(f, x);
			//Function g = new CompiledFunction2(f, x);

			//Console.WriteLine(f.Value(x | 2.0));
			//Console.WriteLine(g.Value(x | 2.0));

			//Console.WriteLine(f.Derivative(x).Value(x | 2.0));
			//Console.WriteLine(g.Derivative(x).Value(x | 2.0));
			for (double t = -5.0; t <= 5.0; t += 0.1)
			{
				//Console.WriteLine("{0:f8} {1:f8} {2}", f.Value(x | t), g.Value(x | t), f.Value(x | t) - g.Value(x | t));
			}
		}

		public static void Sample4()
		{
			// An example showing that FuncLib is very useful even if differentiation or optimization aren't needed.

			// Set up parameters for a two-factor diagonal Vasicek model.
			double kappa_11 = 0.1;
			double kappa_22 = 0.4;
			double y0_1 = 0.4;
			double y0_2 = 0.5;

			// An example from mathematical finance. 

			int n = 100;
			Variable[] x = new Variable[2 * n];
			double dt = 1.0 / n;

			// A diagonal stochastic model.
			FunctionMatrix kappa = new FunctionMatrix(new Function[,] { { 0.1, 0.0 }, { 0.0, 0.4 } });


			// Starting point. Fixed numbers for simplicity.
			//FunctionVector y0 = new FunctionVector(new Function[] { 0.4, 0.5 });


			// Notice that we a keeping the matrix notation to the programmer.
			
			// The initial state vector is just the constant matrix function above.
			//FunctionVector y = y0;

			for (int i = 0; i < n; i++)
			{
				// The Gaussian increments are a two-dimensional vector represented by two consecutive variables.
				//FunctionVector dw = new FunctionVector(new Function[] { x[2 * i], x[2 * i + 1] });

				// Update the state vector, i.e. change reference so that y points to a new updated function vector.
				//y = y + kappa*(theta-y)*dt+sigma*dw;
			}

			// Now y is the reference to the final vector after all the algebra aritmetics.


			// Replace by the compiled vector function. The compilation removes all unnecessary multiplications by zero without the programmer would have to write special code for it.
			//y = Compiler.Compile(y, x);

			// Simulated Gaussian increments.

			//double[] dw = new double[2*n];
			//y.Value(dw);

			// The computation is exactly as fast as if the special diagonal case was implemented manually.
		}
	}
}
