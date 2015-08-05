// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

using FuncLib.Functions;
using FuncLib.Functions.Compilation;

namespace FuncLib.Optimization.Ipopt
{
	/// <summary>
	/// This wrapper class may be used to express knowledge about the sparsity structure of a function in case the exact structure is not inferred automatically.
	/// </summary>
	[Serializable]
	public class SparseFunction : Function
	{
		private Function innerFunction;
		private SparseFunctionGradient gradient;
		private SparseFunctionHessian hessian;

		public SparseFunction(Function innerFunction, params Variable[] gradientEntries)
			: this(innerFunction, new SparseFunctionGradient(gradientEntries))
		{
		}

		public SparseFunction(Function innerFunction, SparseFunctionGradient gradient)
			: this(innerFunction, gradient, null)
		{
		}

		public SparseFunction(Function innerFunction, SparseFunctionGradient gradient, SparseFunctionHessian hessian)
		{
			this.innerFunction = innerFunction;
			if (gradient != null)
			{
				this.gradient = gradient.Clone();
			}
			if (hessian != null)
			{
				this.hessian = hessian.Clone();
			}
		}

		protected override double ComputeValue(IEvaluation evaluation)
		{
			return innerFunction.Value(evaluation);
		}

		protected override Function ComputeDerivative(Variable variable)
		{
			if (gradient != null && !gradient.Contains(variable))
			{
				// This entry is not explicitly specified as non-zero.
				return 0.0;
			}

			// Wrap in SparseFunctionDerivative to take care of the sparsity of the Hessian.
			return new SparseFunctionDerivative(innerFunction.Derivative(variable), variable, gradient, hessian);
		}

		public override Expression Compile(CodeGenerator generator)
		{
			return innerFunction.Compile(generator);
		}

		[Serializable]
		private class SparseFunctionDerivative : Function
		{
			private Function innerFunction;
			private Variable firstVariable;
			private SparseFunctionGradient gradient;
			private SparseFunctionHessian hessian;

			public SparseFunctionDerivative(Function innerFunction, Variable firstVariable, SparseFunctionGradient gradient, SparseFunctionHessian hessian)
			{
				this.innerFunction = innerFunction;
				this.firstVariable = firstVariable;
				this.gradient = gradient;
				this.hessian = hessian;
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				return innerFunction.Value(evaluation);
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				if (gradient != null && !gradient.Contains(variable))
				{
					// This entry is not explicitly specified as non-zero in the gradient.
					return 0.0;
				}

				if (hessian != null && !hessian.Contains(firstVariable, variable))
				{
				    // This entry is not explicitly specified as non-zero in the Hessian.
				    return 0.0;
				}

				// No more wrapping below this level.
				return innerFunction.Derivative(variable);
			}

			public override Expression Compile(CodeGenerator generator)
			{
				return innerFunction.Compile(generator);
			}
		}
	}
}
