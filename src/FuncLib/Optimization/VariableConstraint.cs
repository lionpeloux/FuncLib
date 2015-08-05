// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

using FuncLib.Functions;

namespace FuncLib.Optimization
{
	[Serializable]
	public class VariableConstraint : FunctionConstraint
	{
		private Variable variable;

		public VariableConstraint(Variable variable, double minValue, double maxValue)
			: base(variable, minValue, maxValue)
		{
			this.variable = variable;
		}

		public Variable Variable
		{
            get { return variable; }
		}
	}
}
