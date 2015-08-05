// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;

using FuncLib.Functions;

namespace FuncLib.Optimization
{
	public class VariableCollection : List<Variable>
	{
		public void Add(params Variable[] variables)
		{
			AddRange(variables);
		}
	}
}
