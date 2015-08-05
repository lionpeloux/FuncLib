// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace FuncLib.Functions
{
	/// <summary>
	/// Like <see cref="Evaluation" />, but with support for caching of intermediate function values. This usually speeds up
	/// evaluation of complicated function significantly, but uses more memory during the computation.
	/// </summary>
	[DebuggerStepThrough]
	[DebuggerDisplay("{ToString(),nq}")]
	public class CachedEvaluation : Evaluation
	{
		private Dictionary<Function, double> values;

		public CachedEvaluation(IPoint point)
			: this(point.ToDictionary())
		{
		}

		public CachedEvaluation(IDictionary<Variable, double> assignments)
			: base(assignments)
		{
			values = new Dictionary<Function, double>();
		}

		public override bool TryGetValue(Function function, out double value)
		{
			return values.TryGetValue(function, out value);
		}

		public override void AddValue(Function function, double value)
		{
			values.Add(function, value);
		}
	}
}
