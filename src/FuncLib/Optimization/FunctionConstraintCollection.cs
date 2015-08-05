// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;

using FuncLib.Functions;

namespace FuncLib.Optimization
{
	public class FunctionConstraintCollection : List<FunctionConstraint>
	{
		public void Add(Function function, double minValue, double maxValue)
		{
			Add(new FunctionConstraint(function, minValue, maxValue));
		}

		public void Add(params FunctionConstraint[] constraints)
		{
			AddRange(constraints);
		}
	}
}
