// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;

using FuncLib.Functions;

namespace FuncLib.Optimization
{
	public class VariableConstraintCollection : List<VariableConstraint>
	{
		public void Add(Variable variable, double minValue, double maxValue)
		{
			Add(new VariableConstraint(variable, minValue, maxValue));
		}

		public void Add(params VariableConstraint[] constraints)
		{
			AddRange(constraints);
		}
	}
}
