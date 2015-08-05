// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

namespace FuncLib.Optimization
{
	[Serializable]
	public class ConstraintViolationException : OptimizerException
	{
		private FunctionConstraint constraint;

		public ConstraintViolationException(string message)
			: this(message, null)
		{
		}

		public ConstraintViolationException(string message, FunctionConstraint constraint)
			: base(message)
		{
			this.constraint = constraint;
		}

		public FunctionConstraint Constraint
		{
			get
			{
				return constraint;
			}
		}
	}
}
