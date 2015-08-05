// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;

using FuncLib.Functions;

namespace FuncLib.Optimization
{
	/// <summary>
	/// Base class for non-linear optimizers with lower and upper bounds on variables (i.e. intervals).
	/// </summary>
	[Serializable]
	public abstract class VariableConstrainedOptimizer : Optimizer
	{
		private VariableConstraintCollection variableConstraints;

		public VariableConstrainedOptimizer()
		{
			variableConstraints = new VariableConstraintCollection();
		}

		/// <summary>
		/// Variables added here are bounded to intervals using inequality operator overloads. The == operator overload may also be used;
		/// in that case the variable is removed from the optimization and replaced by the value specified.
		/// </summary>
		public VariableConstraintCollection VariableConstraints
		{
			get
			{
				return variableConstraints;
			}
		}

		[Serializable]
		protected abstract class InnerVariableConstrainedOptimizer : InnerOptimizer
		{
			private List<VariableConstraint> variableConstraints;

			public InnerVariableConstrainedOptimizer(VariableConstrainedOptimizer optimizer)
				: base(optimizer)
			{
				variableConstraints = new List<VariableConstraint>();

				foreach (VariableConstraint constraint in optimizer.VariableConstraints)
				{
					AddVariableConstraint(constraint);
				}
			}

			public void AddVariableConstraint(VariableConstraint constraint)
			{
				Variable variable = constraint.Variable;

				if (constraint is VariableEqualityConstraint)
				{
					// Not really a constraint. Add using base class.
					AddVariableEqualityConstraint((VariableEqualityConstraint)constraint);
					return;
				}

				if (constraint.MinValue == constraint.MaxValue)
				{
					// Not really a constraint. Add using base class.
					AddVariableEqualityConstraint(new VariableEqualityConstraint(variable, constraint.MinValue));
					return;
				}

				if (VariableEqualityConstraints.ContainsKey(variable))
				{
					// Already defined as an equality constraint. Check for violations.
					if (VariableEqualityConstraints[variable] < constraint.MinValue || VariableEqualityConstraints[variable] > constraint.MaxValue)
					{
						throw new InconsistentConstraintException("Inconsistency found in variable constraints.", constraint);
					}

					// No inconsistencies. Just ignore this constraint since the more restrictive equality constraint is already in force.
					return;
				}

				VariableConstraint constraint0 = variableConstraints.Find(delegate(VariableConstraint c) { return c.Variable == variable; });
				if (constraint0 != null)
				{
					// Modify an existing constraint (i.e. regarding the same variable).
					double minValue = Math.Max(constraint0.MinValue, constraint.MinValue);
					double maxValue = Math.Min(constraint0.MaxValue, constraint.MaxValue);

					if (minValue < maxValue)
					{
						// Remove the existing constraint and add a new with a smaller interval.
						variableConstraints.Remove(constraint0);
						variableConstraints.Add(new VariableConstraint(variable, minValue, maxValue));
					}
					else if (minValue == maxValue)
					{
						// Remove the existing constraint and replace it by a variable equality constraint.
						variableConstraints.Remove(constraint0);
						AddVariableEqualityConstraint(new VariableEqualityConstraint(variable, minValue));
					}
					else
					{
						// The two intervals are disjoint.
						throw new InconsistentConstraintException("Inconsistency found in variable constraints.", constraint);
					}
				}
				else
				{
					variableConstraints.Add(constraint);
				}
			}

			public override void CheckConstraints(IPoint point)
			{
				base.CheckConstraints(point);

				foreach (VariableConstraint constraint in variableConstraints)
				{
					Variable variable = constraint.Variable;
					if (point.IsAssigned(variable) && (point[variable] < constraint.MinValue || point[variable] > constraint.MaxValue))
					{
						throw new ConstraintViolationException("A variable constraint isn't satisfied.", constraint);
					}
				}
			}

			public List<VariableConstraint> VariableConstraints
			{
				get
				{
					return variableConstraints;
				}
			}
		}
	}
}
