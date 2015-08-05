// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

using FuncLib.Functions;

namespace FuncLib.Optimization
{
	[Serializable]
	public class FunctionConstraint
	{
		private Function function;
		private double minValue, maxValue;

		public FunctionConstraint(Function function, double minValue, double maxValue)
		{
			this.function = function;
			this.minValue = minValue;
			this.maxValue = maxValue;
		}

		public Function Function
		{
            get { return function; }
		}
		public double MinValue
		{
            get { return minValue; }
		}
		public double MaxValue
		{
            get { return maxValue; }
		}
	}
}
