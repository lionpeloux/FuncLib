// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

namespace FuncLib.Functions
{
	[Serializable]
	public class VariableNotAssignedException : Exception
	{
		private Variable variable;

		public VariableNotAssignedException(Variable variable)
			: this(variable, null)
		{
		}

		public VariableNotAssignedException(Variable variable, string message)
			: base(message)
		{
			this.variable = variable;
		}

		public Variable Variable
		{
			get
			{
				return variable;
			}
		}
	}
}
