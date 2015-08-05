// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;

using FuncLib.Functions;
using FuncLib.Functions.Compilation;

namespace FuncLib.DualNumbers
{
	/// <summary>
	/// Provides a bridge from <see cref="DualNumber" /> to <see cref="Function" />. Implement this abstract class or create an instance using the static <see cref="Create" /> method.
	/// </summary>
	[Serializable]
	public abstract class DualNumberFunction : Function
	{
		private VariableReadOnlyCollection variables;
		private EvaluationCache<DualNumber> values;

		public DualNumberFunction(Variable[] variables)
		{
			// Check for repeated use of the same variable.
			this.variables = new VariableReadOnlyCollection(variables);

			// Important to have the values cached to avoid recomputing when the value, the gradient, and the Hessian are requested after each other.
			values = new EvaluationCache<DualNumber>();
		}

		/// <summary>
		/// Create an instance of <see cref="DualNumberFunction" /> using a delegate.
		/// </summary>
		public static DualNumberFunction Create(Func<IDualNumberTransform, DualNumber> function, params Variable[] variables)
		{
			return new InnerDualNumberFunction(function, variables);
		}

		//public delegate TResult Func<T, TResult>(T arg); // Part of .NET 3+

		protected override double ComputeValue(IEvaluation evaluation)
		{
			return DualValue(evaluation).Value;
		}

		protected abstract DualNumber ComputeDualValue(IDualNumberTransform transform);

		public DualNumber DualValue(IEvaluation evaluation)
		{
			DualNumber value;
			if (!values.TryGetValue(evaluation, out value))
			{
				Dictionary<Variable, DualNumber> transform = new Dictionary<Variable, DualNumber>();
				for (int i = 0; i < variables.Count; i++)
				{
					transform[variables[i]] = DualNumber.Basis(evaluation[variables[i]], variables.Count, i);
				}

				value = ComputeDualValue(new DualNumberTransform(transform, variables));
				values.Add(evaluation, value);
			}

			return value;
		}

		protected override Function ComputeDerivative(Variable variable)
		{
			if (!variables.Contains(variable))
			{
				// No dependency. The derivative is zero.
				return 0.0;
			}

			return new FirstDerivativeFunction(this, variables.IndexOf(variable));
		}

		private Function SecondDerivative(int variableIndex1, Variable variable2)
		{
			if (!variables.Contains(variable2))
			{
				// No dependency. The derivative is zero.
				return 0.0;
			}

			return new SecondDerivativeFunction(this, variableIndex1, variables.IndexOf(variable2));
		}

		public override Expression Compile(CodeGenerator generator)
		{
			// Compilation support would be really sweet!
			throw new NotSupportedException();
		}

		[Serializable]
		private class InnerDualNumberFunction : DualNumberFunction
		{
			private Func<IDualNumberTransform, DualNumber> function;

			public InnerDualNumberFunction(Func<IDualNumberTransform, DualNumber> function, Variable[] variables)
				: base(variables)
			{
				this.function = function;
			}

			protected override DualNumber ComputeDualValue(IDualNumberTransform transform)
			{
				return function(transform);
			}
		}

		[Serializable]
		private class FirstDerivativeFunction : Function
		{
			private DualNumberFunction parent;
			private int variableIndex;

			public FirstDerivativeFunction(DualNumberFunction parent, int variableIndex)
			{
				this.parent = parent;
				this.variableIndex = variableIndex;
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				return parent.DualValue(evaluation).Gradient[variableIndex];
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				return parent.SecondDerivative(variableIndex, variable);
			}

			public override Expression Compile(CodeGenerator generator)
			{
				throw new NotSupportedException();
			}
		}

		[Serializable]
		private class SecondDerivativeFunction : Function
		{
			private DualNumberFunction parent;
			private int variableIndex1, variableIndex2;

			public SecondDerivativeFunction(DualNumberFunction parent, int variableIndex1, int variableIndex2)
			{
				this.parent = parent;
				this.variableIndex1 = variableIndex1;
				this.variableIndex2 = variableIndex2;
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				return parent.DualValue(evaluation).Hessian[variableIndex1, variableIndex2];
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				// Can't compute higher order derivatives this way since only the gradient and Hessian are known.
				throw new NotSupportedException("Can't compute third order derivatives using DualNumberFunction.");
			}

			public override Expression Compile(CodeGenerator generator)
			{
				throw new NotSupportedException();
			}
		}

		[Serializable]
		private class DualNumberTransform : IDualNumberTransform
		{
			private Dictionary<Variable, DualNumber> transform;
			private VariableReadOnlyCollection variables;

			public DualNumberTransform(Dictionary<Variable, DualNumber> transform, VariableReadOnlyCollection variables)
			{
				this.transform = transform;
				this.variables = variables;
			}

			public DualNumber this[Variable variable]
			{
				get
				{
					DualNumber number;
					if (!transform.TryGetValue(variable, out number))
					{
						throw new VariableNotAssignedException(variable);
					}

					return number;
				}
			}

			public VariableReadOnlyCollection Variables
			{
				get
				{
					return variables;
				}
			}
		}
	}
}
