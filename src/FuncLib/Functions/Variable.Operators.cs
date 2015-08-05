// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

using FuncLib.Optimization;

namespace FuncLib.Functions
{
	public partial class Variable
	{
		public static VariableAssignment operator |(Variable variable, double value)
		{
			return new VariableAssignment(variable, value);
		}

		public static VariableConstraint operator <=(Variable variable, double value)
		{
			return new VariableConstraint(variable, double.NegativeInfinity, value);
		}

		public static VariableConstraint operator >=(Variable variable, double value)
		{
			return new VariableConstraint(variable, value, double.PositiveInfinity);
		}

		public static VariableConstraint operator <=(double value, Variable variable)
		{
			return variable >= value;
		}

		public static VariableConstraint operator >=(double value, Variable variable)
		{
			return variable <= value;
		}

		public static VariableEqualityConstraint operator ==(Variable variable, double value)
		{
			return new VariableEqualityConstraint(variable, value);
		}

		public static VariableEqualityConstraint operator !=(Variable variable, double value)
		{
			// This operator makes no sence in this context but is required by C#.
			throw new InvalidOperationException();
		}

		public static VariableEqualityConstraint operator ==(double value, Variable variable)
		{
			return variable == value;
		}

		public static VariableEqualityConstraint operator !=(double value, Variable variable)
		{
			// This operator makes no sence in this context but is required by C#.
			throw new InvalidOperationException();
		}

		public override int GetHashCode()
		{
			// The hash code from Object is fine, but want to get rid of a warning.
			return base.GetHashCode();
		}

		public override bool Equals(object other)
		{
			// The equality from Object is fine, but want to get rid of a warning.
			return base.Equals(other);
		}
	}
}
