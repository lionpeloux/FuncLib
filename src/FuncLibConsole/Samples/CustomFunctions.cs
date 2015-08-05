// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using FuncLib.Functions;
using FuncLib.Functions.Compilation;

namespace FuncLibConsole.Samples
{
	public class CustomFunctions
	{
	}

	public class Atan2 : Function
	{
		private Function f, g;

		public Atan2(Function f, Function g)
		{
			this.f = f;
			this.g = g;
		}

		protected override double ComputeValue(IEvaluation evaluation)
		{
			return Math.Atan2(f.Value(evaluation), g.Value(evaluation));
		}

		protected override Function ComputeDerivative(Variable variable)
		{
			// http://en.wikipedia.org/wiki/Atan2#Derivative
			return (g * f.Derivative(variable) - f * g.Derivative(variable)) / (Function.Sqr(g) + Function.Sqr(f));
		}

		public override Expression Compile(CodeGenerator generator)
		{
			// Of course this only needs to be implemented if this function should support compilation.
			return "Math.Atan2(" + generator.GetReference(f) + "," + generator.GetReference(g) + ")";
		}
	}
}
