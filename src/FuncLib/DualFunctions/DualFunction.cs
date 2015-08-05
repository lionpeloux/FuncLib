// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

using FuncLib.Functions;

namespace FuncLib.DualFunctions
{
	public class DualFunction : Function
	{
		protected override double ComputeValue(IEvaluation evaluation)
		{
			throw new NotImplementedException();
		}

		protected override Function ComputeDerivative(Variable variable)
		{
			throw new NotImplementedException();
		}

		public static DualFunction operator +(DualFunction f, DualFunction g)
		{
			throw new NotImplementedException(); // binary
		}

		public static DualFunction operator +(DualFunction f, double g)
		{
			throw new NotImplementedException(); // unary
		}

		public static DualFunction operator +(double f, DualFunction g)
		{
			throw new NotImplementedException(); // unary
		}
	}
}
