// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

namespace FuncLib.Optimization
{
	[Serializable]
	public class InconsistentConstraintException : ConstraintViolationException
	{
		public InconsistentConstraintException(string message, FunctionConstraint constraint)
			: base(message, constraint)
		{
		}
	}
}
