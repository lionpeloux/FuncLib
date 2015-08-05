// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;

using FuncLib.Functions;

namespace FuncLib.Optimization
{
	/// <summary>
	/// Base class for non-linear optimizers.
	/// </summary>
	[Serializable]
	public abstract class Optimizer : IOptimizer
	{
		private Function objectiveFunction;
		private double objectiveFunctionScaling;
		private VariableCollection variables;
		private VariableEqualityConstraintCollection variableEqualityConstraints;

		public Optimizer()
		{
			objectiveFunctionScaling = 1.0;
			variables = new VariableCollection();
			variableEqualityConstraints = new VariableEqualityConstraintCollection();
		}

		/// <summary>
		/// Run the optimizer with the specified initial values.
		/// </summary>
		public IOptimizerResult Run(IPoint initialPoint)
		{
			// It's easier to use the prepared optimizer object though we're only using it once in this case.
			return Prepare().Run(initialPoint);
		}

		/// <summary>
		/// Run the optimizer with the specified initial values.
		/// </summary>
		public IOptimizerResult Run(params VariableAssignment[] initialAssignments)
		{
			return Run(new Point(initialAssignments));
		}

		/// <summary>
		/// Use the returned object for multiple runs with different initial values.
		/// </summary>
		public abstract PreparedOptimizer Prepare();

		/// <summary>
		/// Test whether all constraints are satisfied.
		/// </summary>
		public bool CheckConstraints(IEvaluation evaluation)
		{
			// Implement this without exceptions? This method isn't needed initially for some optimizers (e.g. Ipopt), and it may also
			// fail slightly for some after optimization (again Ipopt).

			try
			{
				// Again easier with the prepared object.
				((InnerOptimizer)Prepare()).CheckConstraints(evaluation);
			}
			catch (ConstraintViolationException)
			{
				return false;
			}

			return true;
		}

		/// <summary>
		/// The IEvaluation object used to compute the function value before each iteration. Override this to change the default behavior,
		/// i.e. using a different IEvaluation object than CachedEvaluation.
		/// </summary>
		protected virtual IEvaluation CreateEvaluation(IDictionary<Variable, double> assignments)
		{
			return new CachedEvaluation(assignments);
		}

		/// <summary>
		/// The IPoint object returned as the final optimized point in IOptimizerResult. Override this to change the default behavior, i.e.
		/// using a different IPoint object than Point.
		/// </summary>
		protected virtual IPoint CreateOptimalPoint(IDictionary<Variable, double> assignments)
		{
			return new Point(assignments);
		}

		/// <summary>
		/// The function to optimize.
		/// </summary>
		public Function ObjectiveFunction
		{
			get
			{
				return objectiveFunction;
			}
			set
			{
				objectiveFunction = value;
			}
		}

		/// <summary>
		/// Use a positive number to perform minimization (default) or a negative number to perform maximization.
		/// </summary>
		public double ObjectiveFunctionScaling
		{
			get
			{
				return objectiveFunctionScaling;
			}
			set
			{
				objectiveFunctionScaling = value;
			}
		}

		/// <summary>
		/// The unknown variables to optimize.
		/// </summary>
		public VariableCollection Variables
		{
			get
			{
				return variables;
			}
		}

		/// <summary>
		/// Variables added here using the == operator overload are removed from the optimization and replaced by the values specified.
		/// </summary>
		public VariableEqualityConstraintCollection VariableEqualityConstraints
		{
			get
			{
				return variableEqualityConstraints;
			}
		}

		[Serializable]
		protected abstract class InnerOptimizer : PreparedOptimizer
		{
			private Function objectiveFunction;
			private double objectiveFunctionScaling;
			private List<Variable> variables;
			private Dictionary<Variable, double> variableEqualityConstraints;

			public InnerOptimizer(Optimizer optimizer)
				: base(optimizer)
			{
				objectiveFunction = optimizer.ObjectiveFunction;
				objectiveFunctionScaling = optimizer.ObjectiveFunctionScaling;

				if (objectiveFunction == null)
				{
					throw new OptimizerException("No objective function is specified.");
				}

				variables = new List<Variable>();

				foreach (Variable variable in optimizer.Variables)
				{
					// Avoid duplicates. Not an error, just ignore them.
					if (!variables.Contains(variable))
					{
						variables.Add(variable);
					}
				}

				variableEqualityConstraints = new Dictionary<Variable, double>();

				foreach (VariableEqualityConstraint constraint in optimizer.VariableEqualityConstraints)
				{
					AddVariableEqualityConstraint(constraint);
				}
			}

			public void AddVariableEqualityConstraint(VariableEqualityConstraint constraint)
			{
				Variable variable = constraint.Variable;
				double value = constraint.Value;

				if (variableEqualityConstraints.ContainsKey(variable))
				{
					// Previously defined. Check that it's set to the same value.
					if (variableEqualityConstraints[variable] != value)
					{
						throw new InconsistentConstraintException("Inconsistency found in variable equality constraints.", constraint);
					}
				}
				else if (variables.Contains(variable))
				{
					// Remove the variable from the optimization.
					variables.Remove(variable);
					variableEqualityConstraints.Add(variable, value);
				}
				else
				{
					throw new InconsistentConstraintException("A variable equality constraint is specified on a variable that's not included in the optimization.", constraint);
				}
			}

			public virtual void CheckConstraints(IPoint point)
			{
				foreach (KeyValuePair<Variable, double> constraint in variableEqualityConstraints)
				{
					Variable variable = constraint.Key;
					if (point.IsAssigned(variable) && point[variable] != constraint.Value)
					{
						throw new ConstraintViolationException("A variable equality constraint is violated.");
					}
				}
			}

			public IEvaluation CreateEvaluation(double[] values)
			{
				return Optimizer.CreateEvaluation(CreateDictionary(values));
			}

			public IPoint CreateOptimalPoint(double[] values)
			{
				return Optimizer.CreateOptimalPoint(CreateDictionary(values));
			}

			private IDictionary<Variable, double> CreateDictionary(double[] values)
			{
				int n = variables.Count;

				if (values.Length != n)
				{
					throw new OptimizerException();
				}

				Dictionary<Variable, double> assignments = new Dictionary<Variable, double>(variableEqualityConstraints);
				for (int i = 0; i < n; i++)
				{
					assignments[variables[i]] = values[i];
				}

				return assignments;
			}

			public Function ObjectiveFunction
			{
				get
				{
					return objectiveFunction;
				}
			}

			public double ObjectiveFunctionScaling
			{
				get
				{
					return objectiveFunctionScaling;
				}
			}

			public List<Variable> Variables
			{
				get
				{
					return variables;
				}
			}

			public Dictionary<Variable, double> VariableEqualityConstraints
			{
				get
				{
					return variableEqualityConstraints;
				}
			}
		}
	}
}
