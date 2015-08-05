// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;

using FuncLib.Functions.Compilation;

namespace FuncLib.Functions
{
	/// <summary>
	/// Represents a mathematical function. Any function must derive from this class. Also includes a range
	/// of static method for performing operations on functions.
	/// </summary>
	[Serializable]
	public abstract partial class Function
	{
		private Dictionary<Variable, Function> derivatives;

		/// <summary>
		/// Evaluates the function by assigning values to the variables.
		/// </summary>
		public double Value(IPoint point)
		{
			// Evaluate using CachedEvaluation by default (usually much faster, but uses more memory).
			return Value(new CachedEvaluation(point));
		}

		/// <summary>
		/// Evaluates the function by assigning values to the variables.
		/// </summary>
		public double Value(params VariableAssignment[] assignments)
		{
			return Value(new Point(assignments));
		}

		/// <summary>
		/// Evaluates the function by assigning values to the variables.
		/// </summary>
		public double Value(IEvaluation evaluation)
		{
			double value;
			if (!evaluation.TryGetValue(this, out value))
			{
				value = ComputeValue(evaluation);
				evaluation.AddValue(this, value);
			}

			return value;
		}

		/// <summary>
		/// Evaluates the function by assigning values to the variables. Override this to implement.
		/// </summary>
		protected abstract double ComputeValue(IEvaluation evaluation);

		/// <summary>
		/// Computes the partial derivative with respect to a variable.
		/// </summary>
		public Function Derivative(Variable variable)
		{
			// It's has proved very important for efficient computation that the same object is returned
			// if called repeatably. This is especially true if the function is compiled.

			if (derivatives == null)
			{
				derivatives = new Dictionary<Variable, Function>();
			}

			Function derivative;
			if (!derivatives.TryGetValue(variable, out derivative))
			{
				derivative = ComputeDerivative(variable);
				if (!derivative.IsZero)
				{
					// Only use memory for saving non-zero derivatives.
					derivatives.Add(variable, derivative);
				}
			}

			return derivative;
		}

		/// <summary>
		/// Computes the higher order partial derivative with respect to a variable.
		/// </summary>
		public Function Derivative(Variable variable, int order)
		{
			if (order < 0)
			{
				throw new ArgumentOutOfRangeException("order");
			}

			Function f = this;
			for (int i = 0; i < order; i++)
			{
				f = f.Derivative(variable);
			}
			return f;
		}

		/// <summary>
		/// Computes the higher order partial derivative with respect to a number of variables.
		/// </summary>
		public Function Derivative(params Variable[] variables)
		{
			Function f = this;
			for (int i = 0; i < variables.Length; i++)
			{
				f = f.Derivative(variables[i]);
			}
			return f;
		}

		/// <summary>
		/// Computes the partial derivative with respect to a variable. Override this to implement.
		/// </summary>
		protected abstract Function ComputeDerivative(Variable variable);
		
		/// <summary>
		/// Replaces a number of variables by fixed values. Returns a function of the remaining variables.
		/// </summary>
		public virtual Function PartialValue(IPoint point)
		{
			return PartialValue(new PartialValueEvaluation(point));
		}

		/// <summary>
		/// Replace a number of variables by fixed values. Returns a function of the remaining variables.
		/// </summary>
		public Function PartialValue(params VariableAssignment[] assignments)
		{
			return PartialValue(new Point(assignments));
		}

		/// <summary>
		/// Replace a number of variables by fixed values. Returns a function of the remaining variables.
		/// </summary>
		public Function PartialValue(IPartialValueEvaluation evaluation)
		{
			// Reuse the same object if already computed.
			Function partialValue;
			if (!evaluation.TryGetPartialValue(this, out partialValue))
			{
				partialValue = ComputePartialValue(evaluation);
				evaluation.AddPartialValue(this, partialValue);
			}

			return partialValue;
		}

		/// <summary>
		/// Replace a number of variables by fixed values. Override this to implement.
		/// </summary>
		protected virtual Function ComputePartialValue(IPartialValueEvaluation evaluation)
		{
			// Use slow and inefficient implementation using PartialValueFunction by default. Override to provide a more efficient implementation.
			return new PartialValueFunction(this, evaluation);
		}

		/// <summary>
		/// Returns an expression used by <see cref="CodeGenerator" /> to generate C# code representation of the function. Override this to implement.
		/// </summary>
		public virtual Expression Compile(CodeGenerator generator)
		{
			// Assume that compilation is not supported by default.
			throw new NotSupportedException();
		}

		/// <summary>
		/// Tests if this function is known to be identically zero. Override this to express knowledge about zeroness.
		/// </summary>
		public virtual bool IsZero
		{
			get
			{
				ConstantFunction f0 = this as ConstantFunction;
				return f0 != null && f0.Constant == 0.0;
			}
		}

		[Serializable]
		private class PartialValueFunction : Function
		{
			private Function innerFunction;
			private IDictionary<Variable, double> partialValueAssignments;

			public PartialValueFunction(Function innerFunction, IPartialValueEvaluation partialValueEvaluation)
			{
				this.innerFunction = innerFunction;
				
				partialValueAssignments = partialValueEvaluation.ToDictionary();
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				IDictionary<Variable, double> assignments = evaluation.ToDictionary();

				// Overwrite or extend with values defined for the partial value.
				foreach (KeyValuePair<Variable, double> assignment in partialValueAssignments)
				{
					assignments[assignment.Key] = assignment.Value;
				}

				return innerFunction.Value(new Point(assignments));
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				if (partialValueAssignments.ContainsKey(variable))
				{
					// This variable is replaced by a constant.
					return 0.0;
				}

				return innerFunction.Derivative(variable);
			}
		}
	}
}
