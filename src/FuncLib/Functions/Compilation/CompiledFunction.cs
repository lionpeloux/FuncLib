// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;

namespace FuncLib.Functions.Compilation
{
	[Serializable]
	public abstract class CompiledFunction : Function
	{
		private IDictionary<Variable, Function> derivatives;

		public CompiledFunction(IDictionary<Variable, Function> derivatives)
		{
			this.derivatives = derivatives;
		}

		/// <summary>
		/// Fixed variable order evaluation. Fastest evaluation without the IEvaluation abstraction but without any caching.
		/// </summary>
		public abstract double Value(params double[] values);

		protected override Function ComputeDerivative(Variable variable)
		{
			if (derivatives == null)
			{
				// Indicates that derivatives weren't compiled (up to this order).
				throw new NotSupportedException("Higher order derivatives weren't compiled.");
			}

			if (!derivatives.ContainsKey(variable))
			{
				// Doesn't depend on this variable.
				return 0.0;
			}

			return derivatives[variable];
		}

		/// <summary>
		/// The generated C# code used in the compilation.
		/// </summary>
		public abstract string GeneratedCode
		{
			get;
		}
	}
}
