// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Text;

namespace FuncLib.Functions.Compilation
{
	public class ReferenceExpression : Expression
	{
		private Expression innerExpression;
		private int references;
		private string identifier;

		public ReferenceExpression(Expression innerExpression)
		{
			this.innerExpression = innerExpression;

			references = 1;
			identifier = null;
		}

		public void AddReference()
		{
			references++;
		}

		public bool IsInlined()
		{
			//return false; // FIXME
			return references == 1;
		}

		public override string GenerateCode(CodeGenerator generator)
		{
			// Inline if referenced only once.

			if (!IsInlined())
			{
				if (identifier == null)
				{
					identifier = generator.GetNextIdentifier();
				}

				return identifier;
			}
			else
			{
				// inline, add () if not DoubleExpression
				//return innerExpression.GenerateCodeInline(generator);
				return "(" + innerExpression.GenerateCode(generator) + ")";
			}
		}

		public void GenerateCodeStatement(StringBuilder code, CodeGenerator generator)
		{
			if (!IsInlined())
			{
				if (identifier == null)
				{
					identifier = generator.GetNextIdentifier();
				}

				code.AppendLine("double " + identifier + "=" + innerExpression.GenerateCode(generator) + ";");
			}
		}
	}
}
