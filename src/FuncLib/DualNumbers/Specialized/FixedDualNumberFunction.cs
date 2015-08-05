// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

using FuncLib.Functions;
using FuncLib.Functions.Compilation;

namespace FuncLib.DualNumbers.Specialized
{
	/// <summary>
	/// A simple wrapper around an instance of <see cref="DualNumber" />. It's assumed that this function is only evaluated at the same
	/// point as that used for computing the <see cref="DualNumber" />. For efficiency, no checks are performed to enforce this.
	/// </summary>
	[Serializable]
	public class FixedDualNumberFunction : Function // DualNumberWrapper, UnsafeDualNumberFunction, SimpleDualNumberFunction?
	{
		private DualNumber number;
		private VariableReadOnlyCollection variables;

		public FixedDualNumberFunction(DualNumber number, params Variable[] variables)
			: this(number, new VariableReadOnlyCollection(variables))
		{
		}

		public FixedDualNumberFunction(DualNumber number, VariableReadOnlyCollection variables)
		{
			if (number == null)
			{
				throw new ArgumentNullException("number");
			}

			this.number = number;
			this.variables = variables;
		}

		protected override double ComputeValue(IEvaluation evaluation)
		{
			// Don't check values assigned to variables. They're assumed match those that were used for computing the DualNumber.
			return number.Value;
		}

		protected override Function ComputeDerivative(Variable variable)
		{
			if (!variables.Contains(variable))
			{
				// No dependency.
				return 0.0;
			}

			int index = variables.IndexOf(variable);

			if (DerivativeIsConstant(index))
			{
				// Generates more compact compiled code?
				return number.Gradient[index];
			}

			return new DerivativeFunction(this, index);
		}

		private bool DerivativeIsConstant(int index)
		{
			for (int i = 0; i < variables.Count; i++)
			{
				if (number.Hessian[index, i] != 0.0)
				{
					return false;
				}
			}

			return true;
		}

		public override Expression Compile(CodeGenerator generator)
		{
			// Just compile the constant value.
			return ComputeValue(null);
		}

		[Serializable]
		private class DerivativeFunction : Function
		{
			private FixedDualNumberFunction parent;
			private int index;

			public DerivativeFunction(FixedDualNumberFunction parent, int index)
			{
				this.parent = parent;
				this.index = index;
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				return parent.number.Gradient[index];
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				if (!parent.variables.Contains(variable))
				{
					// No dependency.
					return 0.0;
				}

				return parent.number.Hessian[index, parent.variables.IndexOf(variable)];
				//return new SecondDerivativeFunction(parent, index, parent.variables.IndexOf(variable));
			}

			public override Expression Compile(CodeGenerator generator)
			{
				return ComputeValue(null);
			}
		}

		/*[Serializable]
		private class SecondDerivativeFunction : Function
		{
			private FixedDualNumberFunction parent;
			private int index1, index2;

			public SecondDerivativeFunction(FixedDualNumberFunction parent, int index1, int index2)
			{
				this.parent = parent;
				this.index1 = index1;
				this.index2 = index2;
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				return parent.number.Hessian[index1, index2];
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				if (!parent.variables.Contains(variable))
				{
					// No dependency.
					return 0.0;
				}

				// Higher order derivatives aren't supported.
				throw new NotSupportedException();
			}

			public override Expression Compile(CodeGenerator generator)
			{
				return ComputeValue(null);
			}
		}*/
	}
}
