// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;
using System.Diagnostics;

using FuncLib.Functions.Compilation;
using FuncLib.Optimization;
using FuncLib.DualFunctions;

namespace FuncLib.Functions
{
	public partial class Function
	{
		private static Function zero;

		static Function()
		{
			// Used very often (e.g. in matrix initialization). It's also important to keep this special value referring to the same object.
			zero = new ConstantFunction(0.0);
		}

		[DebuggerStepThrough]
		public static implicit operator Function(double a)
		{
			if (a == 0.0)
			{
				return zero;
			}

			return new ConstantFunction(a);
		}

		[DebuggerStepThrough]
		public static Function operator +(Function f, Function g)
		{
			ConstantFunction fc = f as ConstantFunction;
			ConstantFunction gc = g as ConstantFunction;

			if (fc != null && gc != null)
			{
				return fc.Constant + gc.Constant;
			}

			if (fc != null && fc.Constant == 0.0)
			{
				return g;
			}

			if (gc != null && gc.Constant == 0.0)
			{
				return f;
			}

			DualFunction fd = f as DualFunction;
			DualFunction gd = g as DualFunction;

			if (fd != null && gd != null)
			{
				return fd + gd;
			}

			if (fd != null && gc != null)
			{
				return fd + gc.Constant;
			}

			if (fc != null && gd != null)
			{
				return fc.Constant + gd;
			}

			return new AddFunction(f, g);
		}

		[DebuggerStepThrough]
		public static Function operator -(Function f, Function g)
		{
			ConstantFunction f0 = f as ConstantFunction;
			ConstantFunction g0 = g as ConstantFunction;

			if (f0 != null && g0 != null)
			{
				return f0.Constant - g0.Constant;
			}

			if (f0 != null && f0.Constant == 0.0)
			{
				return -g;
			}

			if (g0 != null && g0.Constant == 0.0)
			{
				return f;
			}

			return new SubtractFunction(f, g);
		}

		[DebuggerStepThrough]
		public static Function operator -(Function f)
		{
			ConstantFunction f0 = f as ConstantFunction;

			if (f0 != null)
			{
				return -f0.Constant;
			}

			return new UnaryMinus(f);
		}

		[DebuggerStepThrough]
		public static Function operator *(Function f, Function g)
		{
			ConstantFunction f0 = f as ConstantFunction;
			ConstantFunction g0 = g as ConstantFunction;

			if (f0 != null && g0 != null)
			{
				return f0.Constant * g0.Constant;
			}

			if (f0 != null)
			{
				if (f0.Constant == 0.0)
				{
					return 0.0;
				}
				if (f0.Constant == 1.0)
				{
					return g;
				}

				return new MultiplyConstantFunction(g, f0.Constant);
			}

			if (g0 != null)
			{
				if (g0.Constant == 0.0)
				{
					return 0.0;
				}
				if (g0.Constant == 1.0)
				{
					return f;
				}

				return new MultiplyConstantFunction(f, g0.Constant);
			}

			return new MultiplyFunction(f, g);
		}

		[DebuggerStepThrough]
		public static Function operator /(Function f, Function g)
		{
			ConstantFunction f0 = f as ConstantFunction;
			ConstantFunction g0 = g as ConstantFunction;

			if (f0 != null && g0 != null)
			{
				return f0.Constant / g0.Constant;
			}

			if (f0 != null)
			{
				return new ReciprocalFunction(f0.Constant, g);
			}

			if (g0 != null)
			{
				if (g0.Constant == 1.0)
				{
					return f;
				}

				//return new DivideConstantFunction(f, g0.Constant); ??
			}

			return new DivideFunction(f, g);
		}

		public static FunctionConstraint operator <=(Function f, Function g)
		{
			return new FunctionConstraint(f - g, double.NegativeInfinity, 0.0);
		}

		public static FunctionConstraint operator >=(Function f, Function g)
		{
			return g <= f;
		}

		public static FunctionConstraint operator <=(Function f, double a)
		{
			return new FunctionConstraint(f, double.NegativeInfinity, a);
		}

		public static FunctionConstraint operator >=(Function f, double a)
		{
			return new FunctionConstraint(f, a, double.PositiveInfinity);
		}

		public static FunctionConstraint operator <=(double a, Function f)
		{
			return f >= a;
		}

		public static FunctionConstraint operator >=(double a, Function f)
		{
			return f <= a;
		}

		public static FunctionEqualityConstraint operator ==(Function function, double value)
		{
			return new FunctionEqualityConstraint(function, value);
		}

		public static FunctionEqualityConstraint operator !=(Function function, double value)
		{
			// This operator makes no sence in this context but is required by C#.
			throw new InvalidOperationException();
		}

		public static FunctionEqualityConstraint operator ==(double value, Function variable)
		{
			return variable == value;
		}

		public static FunctionEqualityConstraint operator !=(double value, Function variable)
		{
			// This operator makes no sence in this context but is required by C#.
			throw new InvalidOperationException();
		}

		/// <summary>
		/// Computes the sum of one or more functions. Faster than adding many functions using the + operator overload.
		/// </summary>
		public static Function Sum(params Function[] functions)
		{
			Function result;
			if (!TrySumConstantFunction(functions, out result) && !TrySumDualFunction(functions, out result))
			{
				result = new SumFunction(functions);
			}

			return result;
		}

		private static bool TrySumDualFunction(Function[] functions, out Function result)
		{
			// It's faster to do pairwise binary operations for DualFunction instances.
			result = 0.0;
			//result += ...
			return false;
		}

		private static bool TrySumConstantFunction(Function[] functions, out Function result)
		{
			double a = 0.0;
			foreach (Function function in functions)
			{
				if (function == null || !(function is ConstantFunction))
				{
					result = null;
					return false;
				}

				a += ((ConstantFunction)function).Constant;
			}

			result = a;
			return true;
		}

		/// <summary>
		/// Computes the product of one or more functions. Faster than multiplying many functions using the * operator overload.
		/// </summary>
		//public static Function Product(params Function[] functions)
		//{
		//    throw new NotImplementedException();
		//}

		[Serializable]
		private class AddFunction : Function
		{
			private Function f, g;

			public AddFunction(Function f, Function g)
			{
				this.f = f;
				this.g = g;
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				return f.Value(evaluation) + g.Value(evaluation);
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				return f.Derivative(variable) + g.Derivative(variable);
			}

			protected override Function ComputePartialValue(IPartialValueEvaluation evaluation)
			{
				return f.PartialValue(evaluation) + g.PartialValue(evaluation);
			}

			public override Expression Compile(CodeGenerator generator)
			{
				return generator.GetReference(f) + "+" + generator.GetReference(g);
			}
		}

		[Serializable]
		private class SubtractFunction : Function
		{
			private Function f, g;

			public SubtractFunction(Function f, Function g)
			{
				this.f = f;
				this.g = g;
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				return f.Value(evaluation) - g.Value(evaluation);
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				return f.Derivative(variable) - g.Derivative(variable);
			}

			protected override Function ComputePartialValue(IPartialValueEvaluation evaluation)
			{
				return f.PartialValue(evaluation) - g.PartialValue(evaluation);
			}

			public override Expression Compile(CodeGenerator generator)
			{
				return generator.GetReference(f) + "-" + generator.GetReference(g);
			}
		}

		[Serializable]
		private class UnaryMinus : Function
		{
			private Function f;

			public UnaryMinus(Function f)
			{
				this.f = f;
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				return -f.Value(evaluation);
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				return -f.Derivative(variable);
			}

			protected override Function ComputePartialValue(IPartialValueEvaluation evaluation)
			{
				return -f.PartialValue(evaluation);
			}

			public override Expression Compile(CodeGenerator generator)
			{
				return "-" + generator.GetReference(f);
			}
		}

		[Serializable]
		private class MultiplyFunction : Function
		{
			private Function f, g;

			public MultiplyFunction(Function f, Function g)
			{
				this.f = f;
				this.g = g;
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				return f.Value(evaluation) * g.Value(evaluation);
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				return f.Derivative(variable) * g + f * g.Derivative(variable);
			}

			protected override Function ComputePartialValue(IPartialValueEvaluation evaluation)
			{
				return f.PartialValue(evaluation) * g.PartialValue(evaluation);
			}

			public override Expression Compile(CodeGenerator generator)
			{
				return generator.GetReference(f) + "*" + generator.GetReference(g);
			}
		}

		[Serializable]
		private class MultiplyConstantFunction : Function
		{
			private Function f;
			private double a;

			public MultiplyConstantFunction(Function f, double a)
			{
				this.f = f;
				this.a = a;
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				return a * f.Value(evaluation);
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				return a * f.Derivative(variable);
			}

			protected override Function ComputePartialValue(IPartialValueEvaluation evaluation)
			{
				return a * f.PartialValue(evaluation);
			}

			public override Expression Compile(CodeGenerator generator)
			{
				// First + operator creates a string without the explicit conversion.
				return (Expression)a + "*" + generator.GetReference(f);
			}
		}

		[Serializable]
		private class DivideFunction : Function
		{
			private Function f, g;

			public DivideFunction(Function f, Function g)
			{
				this.f = f;
				this.g = g;
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				return f.Value(evaluation) / g.Value(evaluation);
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				//return (f.Derivative(variable) * g - f * g.Derivative(variable)) / Function.Sqr(g);
				return f.Derivative(variable) / g - f * g.Derivative(variable) / Function.Sqr(g);
			}

			protected override Function ComputePartialValue(IPartialValueEvaluation evaluation)
			{
				return f.PartialValue(evaluation) / g.PartialValue(evaluation);
			}

			public override Expression Compile(CodeGenerator generator)
			{
				return generator.GetReference(f) + "/" + generator.GetReference(g);
			}
		}

		[Serializable]
		private class ReciprocalFunction : Function
		{
			private double a;
			private Function f;

			public ReciprocalFunction(double a, Function f)
			{
				this.a = a;
				this.f = f;
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				// Only really the reciprocal function when a=1.
				return a / f.Value(evaluation);
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				return a * -f.Derivative(variable) / Function.Sqr(f);
			}

			protected override Function ComputePartialValue(IPartialValueEvaluation evaluation)
			{
				return a / f.PartialValue(evaluation);
			}

			public override Expression Compile(CodeGenerator generator)
			{
				// First + operator creates a string without the explicit conversion.
				return (Expression)a + "/" + generator.GetReference(f);
			}
		}

		[Serializable]
		private class SumFunction : Function
		{
			private Function[] functions;

			public SumFunction(params Function[] functions)
			{
				this.functions = (Function[])functions.Clone();
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				double value = 0.0;
				foreach (Function function in functions)
				{
					value += function.Value(evaluation);
				}

				return value;
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				List<Function> derivatives = new List<Function>();
				foreach (Function function in functions)
				{
					derivatives.Add(function.Derivative(variable));
				}

				return new SumFunction(derivatives.ToArray());
			}

			public override Expression Compile(CodeGenerator generator)
			{
				return base.Compile(generator);
			}
		}
	}
}
