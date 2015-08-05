// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

using FuncLib.Functions;

namespace FuncLib.Optimization
{
	[Serializable]
	public class FunctionEqualityConstraint : FunctionConstraint
	{
		private double value;

		public FunctionEqualityConstraint(Function function, double value)
			: base(function, value, value)
		{
			this.value = value;
		}

		public double Value
		{
            get { return value; }
		}
	}
}
