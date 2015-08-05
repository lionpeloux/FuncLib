// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

using FuncLib.Functions.Compilation;

namespace FuncLib.Functions
{
	public partial class Function
	{
		public static Function Exp(Function f)
		{
			ConstantFunction f0 = f as ConstantFunction;

			if (f0 != null)
			{
				return Math.Exp(f0.Constant);
			}

			return new ExpFunction(f);
		}

		public static Function Log(Function f)
		{
			ConstantFunction f0 = f as ConstantFunction;

			if (f0 != null)
			{
				return Math.Log(f0.Constant);
			}

			return new LogFunction(f);
		}

		public static Function Sqr(Function f)
		{
			ConstantFunction f0 = f as ConstantFunction;

			if (f0 != null)
			{
				return f0.Constant * f0.Constant;
			}

			return new SqrFunction(f);
		}

		public static Function Sqrt(Function f)
		{
			ConstantFunction f0 = f as ConstantFunction;

			if (f0 != null)
			{
				return Math.Sqrt(f0.Constant);
			}

			return new SqrtFunction(f);
		}

		public static Function Pow(Function f, Function g)
		{
			ConstantFunction f0 = f as ConstantFunction;
			ConstantFunction g0 = g as ConstantFunction;

			if (f0 != null && g0 != null)
			{
				return Math.Pow(f0.Constant, g0.Constant);
			}

			if (f0 != null)
			{
				// ...
			}

			if (g0 != null)
			{
				return new PowConstantFunction(f, g0.Constant);
			}

			return new PowFunction(f, g);
		}

		public static Function Cos(Function f)
		{
			ConstantFunction f0 = f as ConstantFunction;

			if (f0 != null)
			{
				return Math.Cos(f0.Constant);
			}

			return new CosFunction(f);
		}

		public static Function Sin(Function f)
		{
			ConstantFunction f0 = f as ConstantFunction;

			if (f0 != null)
			{
				return Math.Sin(f0.Constant);
			}

			return new SinFunction(f);
		}

		[Serializable]
		private class ExpFunction : Function
		{
			private Function f;

			public ExpFunction(Function f)
			{
				this.f = f;
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				return Math.Exp(f.Value(evaluation));
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				// Reuse the same object here.
				return this * f.Derivative(variable);
			}

			protected override Function ComputePartialValue(IPartialValueEvaluation evaluation)
			{
				return Function.Exp(f.PartialValue(evaluation));
			}

			public override Expression Compile(CodeGenerator generator)
			{
				return "Math.Exp(" + generator.GetReference(f) + ")";
			}
		}

		[Serializable]
		private class LogFunction : Function
		{
			private Function f;

			public LogFunction(Function f)
			{
				this.f = f;
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				return Math.Log(f.Value(evaluation));
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				return 1.0 / f * f.Derivative(variable);
			}

			protected override Function ComputePartialValue(IPartialValueEvaluation evaluation)
			{
				return Function.Log(f.PartialValue(evaluation));
			}

			public override Expression Compile(CodeGenerator generator)
			{
				return "Math.Log(" + generator.GetReference(f) + ")";
			}
		}

		[Serializable]
		private class SqrFunction : Function
		{
			private Function f;

			public SqrFunction(Function f)
			{
				this.f = f;
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				double x = f.Value(evaluation);
				return x * x;
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				return 2.0 * f * f.Derivative(variable);
			}

			protected override Function ComputePartialValue(IPartialValueEvaluation evaluation)
			{
				return Function.Sqr(f.PartialValue(evaluation));
			}

			public override Expression Compile(CodeGenerator generator)
			{
				return generator.GetReference(f) + "*" + generator.GetReference(f);
			}
		}

		[Serializable]
		private class SqrtFunction : Function
		{
			private Function f;

			public SqrtFunction(Function f)
			{
				this.f = f;
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				return Math.Sqrt(f.Value(evaluation));
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				// Reuse the same object here.
				return 0.5 / this * f.Derivative(variable);
			}

			protected override Function ComputePartialValue(IPartialValueEvaluation evaluation)
			{
				return Function.Sqrt(f.PartialValue(evaluation));
			}

			public override Expression Compile(CodeGenerator generator)
			{
				return "Math.Sqrt(" + generator.GetReference(f) + ")";
			}
		}

		[Serializable]
		private class PowFunction : Function
		{
			private Function f, g;

			public PowFunction(Function f, Function g)
			{
				this.f = f;
				this.g = g;
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				return Math.Pow(f.Value(evaluation), g.Value(evaluation));
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				return g * Function.Pow(f, g - 1.0) * f.Derivative(variable) + Function.Log(f) * this * g.Derivative(variable);
			}

			protected override Function ComputePartialValue(IPartialValueEvaluation evaluation)
			{
				return Function.Pow(f.PartialValue(evaluation), g.PartialValue(evaluation));
			}

			public override Expression Compile(CodeGenerator generator)
			{
				return "Math.Pow(" + generator.GetReference(f) + "," + generator.GetReference(g) + ")";
			}
		}

		[Serializable]
		private class PowConstantFunction : Function
		{
			private Function f;
			private double a;

			public PowConstantFunction(Function f, double a)
			{
				this.f = f;
				this.a = a;
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				return Math.Pow(f.Value(evaluation), a);
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				return a * Function.Pow(f, a - 1.0) * f.Derivative(variable);
			}

			protected override Function ComputePartialValue(IPartialValueEvaluation evaluation)
			{
				return Function.Pow(f.PartialValue(evaluation), a);
			}

			public override Expression Compile(CodeGenerator generator)
			{
				return "Math.Pow(" + generator.GetReference(f) + "," + a + ")";
			}
		}

		[Serializable]
		private class CosFunction : Function
		{
			private Function f;
			private SinFunction sin;

			public CosFunction(Function f)
				: this(f, null)
			{
			}

			public CosFunction(Function f, SinFunction sin)
			{
				this.f = f;
				this.sin = sin;
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				return Math.Cos(f.Value(evaluation));
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				if (sin == null)
				{
					sin = new SinFunction(f, this);
				}

				return -sin * f.Derivative(variable);
			}

			protected override Function ComputePartialValue(IPartialValueEvaluation evaluation)
			{
				return Function.Cos(f.PartialValue(evaluation));
			}

			public override Expression Compile(CodeGenerator generator)
			{
				return "Math.Cos(" + generator.GetReference(f) + ")";
			}
		}

		[Serializable]
		private class SinFunction : Function
		{
			private Function f;
			private CosFunction cos;

			public SinFunction(Function f)
				: this(f, null)
			{
			}

			public SinFunction(Function f, CosFunction cos)
			{
				this.f = f;
				this.cos = cos;
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				return Math.Sin(f.Value(evaluation));
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				if (cos == null)
				{
					cos = new CosFunction(f, this);
				}

				return cos * f.Derivative(variable);
			}

			protected override Function ComputePartialValue(IPartialValueEvaluation evaluation)
			{
				return Function.Sin(f.PartialValue(evaluation));
			}

			public override Expression Compile(CodeGenerator generator)
			{
				return "Math.Sin(" + generator.GetReference(f) + ")";
			}
		}
	}
}
