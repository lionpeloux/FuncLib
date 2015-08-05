// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;

using FuncLib.Functions;

namespace FuncLib.Optimization.Ipopt
{
	[Serializable]
	public class SparseFunctionGradient : List<Variable>
	{
		public SparseFunctionGradient()
		{
		}

		public SparseFunctionGradient(params Variable[] entries)
			: base(entries)
		{
		}

		private SparseFunctionGradient(SparseFunctionGradient gradient)
			: base(gradient)
		{
		}

		public SparseFunctionGradient Clone()
		{
			return new SparseFunctionGradient(this);
		}
	}
}
