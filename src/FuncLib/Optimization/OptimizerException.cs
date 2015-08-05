// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

namespace FuncLib.Optimization
{
	[Serializable]
	public class OptimizerException : Exception
	{
		public OptimizerException()
		{
		}

		public OptimizerException(string message)
			: base(message)
		{
		}
	}
}
