// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

using FuncLib.Functions;

namespace FuncLib.Optimization
{
	/// <summary>
	/// This interface represents the most basic optimizer, including prepared instances of optimizers. The only
	/// requirement is that it's able to transform an initial point into <see cref="IOptimizerResult" />.
	/// </summary>
	public interface IOptimizer
	{
		/// <summary>
		/// Run the optimizer with the specified initial values.
		/// </summary>
		IOptimizerResult Run(IPoint initialPoint);

		/// <summary>
		/// Run the optimizer with the specified initial values.
		/// </summary>
		IOptimizerResult Run(params VariableAssignment[] initialAssignments);
	}
}
