// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

using FuncLib.Functions;

namespace FuncLib.Optimization
{
	public class VariableEqualityConstraint : VariableConstraint
	{
		private double value;

		public VariableEqualityConstraint(Variable variable, double value)
			: base(variable, value, value)
		{
			this.value = value;
		}

		public double Value
		{
			get
			{
				return value;
			}
		}
	}
}
