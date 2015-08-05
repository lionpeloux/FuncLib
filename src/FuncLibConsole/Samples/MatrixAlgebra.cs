using System;

using FuncLib.Functions;
using FuncLib.Functions.Compilation;

namespace FuncLibConsole.Samples
{
	public class MatrixAlgebra
	{
		public static void Sample1()
		{
			Variable theta = new Variable();

			// Matrix rotation.
			FunctionMatrix r = new FunctionMatrix(new Function[,] { { Function.Cos(theta), -Function.Sin(theta) }, { Function.Sin(theta), Function.Cos(theta) } });

			Console.WriteLine("Rotations");
			Console.WriteLine(r.Value(theta | 0.0).ToString("f4"));
			Console.WriteLine(r.Value(theta | 0.5 * Math.PI).ToString("f4"));
			Console.WriteLine(r.Value(theta | Math.PI).ToString("f4"));
			Console.WriteLine(r.Value(theta | 1.5 * Math.PI).ToString("f4"));
			Console.WriteLine(r.Value(theta | 2.0 * Math.PI).ToString("f4"));
			Console.WriteLine();

			Console.WriteLine("Inverse rotations");
			FunctionMatrix rinv = FunctionMatrix.Inverse(r);
			Console.WriteLine(rinv.Value(theta | 0.0).ToString("f4"));
			Console.WriteLine(rinv.Value(theta | 0.5 * Math.PI).ToString("f4"));
			Console.WriteLine(rinv.Value(theta | Math.PI).ToString("f4"));
			Console.WriteLine(rinv.Value(theta | 1.5 * Math.PI).ToString("f4"));
			Console.WriteLine(rinv.Value(theta | 2.0 * Math.PI).ToString("f4"));
			Console.WriteLine();

			FunctionMatrix rrinv = r * rinv;
			Console.WriteLine("Rotation * Inverse rotations (including first and second derivatives)");
			Console.WriteLine(rrinv.Value(theta | 0.5 * Math.PI));
			Console.WriteLine(rrinv.Derivative(theta).Value(theta | 0.5 * Math.PI));
			Console.WriteLine(rrinv.Derivative(theta).Derivative(theta).Value(theta | 0.5 * Math.PI));
		}

		public static void Sample2()
		{
			// A really crazy sample to show the versatile capabilities of FuncLib.

			// Variables.
			Variable x = new Variable();
			Variable y = new Variable();

			// Define function.
			Function f = x * Function.Exp(x * y);
			Function g = f.Derivative(x);

			// A matrix with mixed functions and constants.
			FunctionMatrix a = new FunctionMatrix(new Function[,] { { f, g }, { 1.0, f * g } });

			// Compute determinant of the inverse's entry-wise derivative.
			Function d = FunctionMatrix.Determinant(FunctionMatrix.Inverse(a).Derivative(y));

			// Use in a new matrix.
			FunctionMatrix b = new FunctionMatrix(new Function[,] { { d, f }, { 1.0, -d } });

			// Compute inverse again and take the first entry.
			Function h = FunctionMatrix.Inverse(b)[0, 0];

			Console.WriteLine(h.Value(x | 0.4, y | 1.6));

			// Compile this complicated function to remove the potentially slow C# object-structure.
			h = Compiler.Compile(h, x, y);

			Console.WriteLine(h.Value(x | 0.4, y | 1.6));

			// Show generated C# code. Notice that the slow exponential function is only evaluated once.
			Console.WriteLine(((CompiledFunction)h).GeneratedCode);
		}
	}
}
