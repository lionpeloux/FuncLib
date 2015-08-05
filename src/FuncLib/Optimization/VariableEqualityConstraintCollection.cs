// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;

using FuncLib.Functions;

namespace FuncLib.Optimization
{
	public class VariableEqualityConstraintCollection : List<VariableEqualityConstraint>
	{
		public void Add(Variable variable, double value)
		{
			Add(new VariableEqualityConstraint(variable, value));
		}

		public void Add(params VariableEqualityConstraint[] constraints)
		{
			AddRange(constraints);
		}
	}
}
