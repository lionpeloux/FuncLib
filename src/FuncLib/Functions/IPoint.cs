// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;

namespace FuncLib.Functions
{
	/// <summary>
	/// Represents a set of variables with some values assigned.
	/// </summary>
	public interface IPoint
	{
		/// <summary>
		/// Test if a value is assigned to a variable.
		/// </summary>
		bool IsAssigned(Variable variable);

		IDictionary<Variable, double> ToDictionary();

		/// <summary>
		/// Get the value assigned to a variable. Throws VariableNotAssignedException if not present.
		/// </summary>
		double this[Variable variable]
		{
			get;
		}
	}
}
