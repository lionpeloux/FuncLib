// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace FuncLib.Functions
{
	/// <summary>
	/// A simple implementation of <see chref="IEvaluation" /> without support for caching of intermediate function values.
	/// This is basically just an immutable wrapper around <see cref="Dictionary<Variable, double>" />.
	/// </summary>
	[DebuggerStepThrough]
	[DebuggerDisplay("{ToString(),nq}")]
	public class Evaluation : Point, IEvaluation
	{
		public Evaluation(IPoint point)
			: this(point.ToDictionary())
		{
		}

		public Evaluation(IDictionary<Variable, double> assignments)
			: base(assignments)
		{
		}

		public virtual bool TryGetValue(Function function, out double value)
		{
			// No caching.
			value = double.NaN;
			return false;
		}

		public virtual void AddValue(Function function, double value)
		{
		}
	}
}
