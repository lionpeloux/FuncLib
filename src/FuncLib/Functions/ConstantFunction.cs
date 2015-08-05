// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Diagnostics;

using FuncLib.Functions.Compilation;

namespace FuncLib.Functions
{
	[Serializable]
	[DebuggerStepThrough]
	[DebuggerDisplay("{Constant}")]
	public class ConstantFunction : Function
	{
		private double a;

		public ConstantFunction(double a)
		{
			this.a = a;
		}

		protected override double ComputeValue(IEvaluation evaluation)
		{
			return a;
		}

		protected override Function ComputeDerivative(Variable variable)
		{
			return 0.0;
		}

		protected override Function ComputePartialValue(IPartialValueEvaluation evaluation)
		{
			return this;
		}

		public override Expression Compile(CodeGenerator generator)
		{
			return a;
		}

		public double Constant
		{
			get
			{
				return a;
			}
		}
	}
}
