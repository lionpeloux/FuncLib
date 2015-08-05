// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;

namespace FuncLib.Functions
{
	public class PartialValueEvaluation : Point, IPartialValueEvaluation
	{
		private Dictionary<Function, Function> partialValues;

		public PartialValueEvaluation(IPoint point)
			: this(point.ToDictionary())
		{
		}

		public PartialValueEvaluation(IDictionary<Variable, double> assignments)
			: base(assignments)
		{
			partialValues = new Dictionary<Function, Function>();
		}

		public bool TryGetPartialValue(Function function, out Function partialValue)
		{
			return partialValues.TryGetValue(function, out partialValue);
		}

		public void AddPartialValue(Function function, Function partialValue)
		{
			partialValues.Add(function, partialValue);
		}
	}
}
