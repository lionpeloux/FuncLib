// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Diagnostics;

namespace FuncLib.Functions
{
	[Serializable]
	[DebuggerStepThrough]
	public class VariableAssignment
	{
		private Variable variable;
		private double value;

		public VariableAssignment(Variable variable, double value)
		{
			this.variable = variable;
			this.value = value;
		}

		public Variable Variable
		{
			get
			{
				return variable;
			}
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
