// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Diagnostics;

using FuncLib.Functions.Compilation;

namespace FuncLib.Functions
{
	/// <summary>
	/// Represents a variable and the identity function.
	/// </summary>
	[Serializable]
	[DebuggerStepThrough]
	[DebuggerDisplay("{DebuggerDisplay,nq}")]
	public partial class Variable : Function
	{
		private string name;

		public Variable()
			: this(null)
		{
		}

		public Variable(string name)
		{
			this.name = name;
		}

		protected override double ComputeValue(IEvaluation evaluation)
		{
			return evaluation[this];
		}

		protected override Function ComputeDerivative(Variable variable)
		{
			return variable == this ? 1.0 : 0.0;
		}

		protected override Function ComputePartialValue(IPartialValueEvaluation evaluation)
		{
			if (evaluation.IsAssigned(this))
			{
				// Replace this variable by the specified value.
				return evaluation[this];
			}

			// This variable isn't replaced by a value. Just keep it.
			return this;
		}

		public override Expression Compile(CodeGenerator generator)
		{
			return generator.GetVariableReference(this);
		}

		public string Name
		{
			get
			{
				return name;
			}
		}

		private string DebuggerDisplay
		{
			get
			{
				// Use pattern from here: http://blogs.msdn.com/b/jaredpar/archive/2011/03/18/debuggerdisplay-attribute-best-practices.aspx
				return name ?? "{" + GetType().FullName + "}";
			}
		}
	}
}
