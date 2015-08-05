// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

using FuncLib.Mathematics;

namespace FuncLib.Functions
{
	[Serializable]
	public sealed class ComplexFunction
	{
		private Function re, im;

		public ComplexFunction(Function re, Function im)
		{
			this.re = re;
			this.im = im;
		}

		public Complex Value(IPoint point)
		{
			return new Complex(re.Value(point), im.Value(point));
		}

		public Complex Value(params VariableAssignment[] assignments)
		{
			return Value(new Point(assignments));
		}

		public static ComplexFunction operator +(ComplexFunction z1, ComplexFunction z2)
		{
			return new ComplexFunction(z1.re + z2.re, z1.im + z2.im);
		}

		public static ComplexFunction operator -(ComplexFunction z)
		{
			return new ComplexFunction(-z.re, -z.im);
		}

		public static ComplexFunction operator -(ComplexFunction z1, ComplexFunction z2)
		{
			return new ComplexFunction(z1.re - z2.re, z1.im - z2.im);
		}

		public static ComplexFunction operator *(ComplexFunction z1, ComplexFunction z2)
		{
			return new ComplexFunction(z1.re * z2.re - z1.im * z2.im, z1.im * z2.re + z1.re * z2.im);
		}

		public static ComplexFunction operator /(ComplexFunction z1, ComplexFunction z2)
		{
			Function c = z2.re * z2.re + z2.im * z2.im;
			return new ComplexFunction((z1.re * z2.re + z1.im * z2.im) / c, (z1.im * z2.re - z1.re * z2.im) / c);
		}

		public static implicit operator ComplexFunction(Function a)
		{
			return new ComplexFunction(a, 0.0);
		}

		public static implicit operator ComplexFunction(Complex z)
		{
			return new ComplexFunction(Complex.Re(z), Complex.Im(z));
		}

		public static implicit operator ComplexFunction(double a)
		{
			return new ComplexFunction(a, 0.0);
		}

		public static Function Re(ComplexFunction z)
		{
			return z.re;
		}

		public static Function Im(ComplexFunction z)
		{
			return z.im;
		}

		public static ComplexFunction Conjugate(ComplexFunction z)
		{
			return new ComplexFunction(z.re, -z.im);
		}

		public static Function Abs(ComplexFunction z)
		{
			return Function.Sqrt(z.re * z.re + z.im * z.im);
		}

		public static ComplexFunction Exp(ComplexFunction z)
		{
			Function c = Function.Exp(z.re);
			return new ComplexFunction(c * Function.Cos(z.im), c * Function.Sin(z.im));
		}

		public static ComplexFunction Sqr(ComplexFunction z)
		{
			return new ComplexFunction(z.re * z.re - z.im * z.im, 2.0 * z.re * z.im);
		}

		public static ComplexFunction Sqrt(ComplexFunction z)
		{
			ReSqrtFunction re = new ReSqrtFunction(z);
			ImSqrtFunction im = new ImSqrtFunction(z);
			ComplexFunction z05 = new ComplexFunction(re, im);
			re.Update(z05);
			im.Update(z05);
			return z05;
		}

		public static Function Arg(ComplexFunction z)
		{
			// Use Complex.Arg through the inner class defined below.
			return new ArgFunction(z);
		}

		public static ComplexFunction Log(ComplexFunction z)
		{
			return new ComplexFunction(0.5 * Function.Log(z.re * z.re + z.im * z.im), Arg(z));
		}

		public static ComplexFunction I
		{
			get
			{
				return new ComplexFunction(0.0, 1.0);
			}
		}

		[Serializable]
		private class ReSqrtFunction : Function
		{
			private ComplexFunction z, z05;

			public ReSqrtFunction(ComplexFunction z)
			{
				this.z = z;
			}

			public void Update(ComplexFunction z05)
			{
				this.z05 = z05;
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				return Complex.Re(Complex.Sqrt(z.Value(evaluation)));
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				throw new NotImplementedException();
			}
		}

		[Serializable]
		private class ImSqrtFunction : Function
		{
			private ComplexFunction z, z05;

			public ImSqrtFunction(ComplexFunction z)
			{
				this.z = z;
			}

			public void Update(ComplexFunction z05)
			{
				this.z05 = z05;
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				return Complex.Im(Complex.Sqrt(z.Value(evaluation)));
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				throw new NotImplementedException();
			}
		}

		[Serializable]
		private class ArgFunction : Function
		{
			private ComplexFunction z;

			public ArgFunction(ComplexFunction z)
			{
				this.z = z;
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				return Complex.Arg(z.Value(evaluation));
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				return (Function.Sqr(z.re) * z.im.Derivative(variable) - z.im * z.re.Derivative(variable)) / (Function.Sqr(z.re) + Function.Sqr(z.im));
			}

			public override Compilation.Expression Compile(Compilation.CodeGenerator compiler)
			{
				// Expression from Complex.Arg? Probably rather messy.
				throw new NotSupportedException();
			}
		}
	}
}
