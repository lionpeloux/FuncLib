// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;

using FuncLib.Functions;

namespace FuncLib.Optimization
{
	/// <summary>
	/// Base class for non-linear optimizers with non-linear constraints.
	/// </summary>
	[Serializable]
	public abstract class ConstrainedOptimizer : VariableConstrainedOptimizer
	{
		private FunctionConstraintCollection constraints;

		public ConstrainedOptimizer()
		{
			constraints = new FunctionConstraintCollection();
		}

		/// <summary>
		/// Any function or variable may be added here to constrain the optimization. Functions are added using inequality operator overloads.
		/// </summary>
		public FunctionConstraintCollection Constraints
		{
            get { return constraints; }
		}

		[Serializable]
		protected abstract class InnerConstrainedOptimizer : InnerVariableConstrainedOptimizer
		{
			private List<FunctionConstraint> functionConstraints;

			public InnerConstrainedOptimizer(ConstrainedOptimizer optimizer)
				: base(optimizer)
			{
				functionConstraints = new List<FunctionConstraint>();

				foreach (FunctionConstraint constraint in optimizer.Constraints)
				{
					AddFunctionConstraint(constraint);
				}
			}

			public void AddFunctionConstraint(FunctionConstraint constraint)
			{
				Function function = constraint.Function;

				if (constraint is VariableConstraint)
				{
					// Not really a function constraint. Let the base class take care of it.
					AddVariableConstraint((VariableConstraint)constraint);
					return;
				}

				if (function is Variable)
				{
					// Same as above.
					AddVariableConstraint(new VariableConstraint((Variable)function, constraint.MinValue, constraint.MaxValue));
					return;
				}

				FunctionConstraint constraint0 = functionConstraints.Find(delegate(FunctionConstraint c) { return c.Function == function; });
				if (constraint0 != null)
				{
					// Modify an existing constraint (i.e. regarding the same variable).
					double minValue = Math.Max(constraint0.MinValue, constraint.MinValue);
					double maxValue = Math.Min(constraint0.MaxValue, constraint.MaxValue);

					if (minValue <= maxValue)
					{
						// Remove the existing constraint and add a new with a smaller interval.
						functionConstraints.Remove(constraint0);
						functionConstraints.Add(new FunctionConstraint(function, minValue, maxValue));
					}
					else
					{
						// The two intervals are disjoint.
						throw new InconsistentConstraintException("Inconsistency found in constraints.", constraint);
					}
				}
				else
				{
					functionConstraints.Add(constraint);
				}
			}

			public override void CheckConstraints(IPoint point)
			{
				base.CheckConstraints(point);

				foreach (FunctionConstraint constraint in functionConstraints)
				{
					double value = constraint.Function.Value(point);
					if (value < constraint.MinValue || value > constraint.MaxValue)
					{
						throw new ConstraintViolationException("A constraint isn't satisfied.", constraint);
					}
				}
			}

			public List<FunctionConstraint> FunctionConstraints
			{
                get { return functionConstraints; }
			}
		}
	}
}
