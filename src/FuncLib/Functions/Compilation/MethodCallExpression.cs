// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;
using System.Text;

namespace FuncLib.Functions.Compilation
{
	public class MethodCallExpression : Expression
	{
		private string methodName;
		private List<ReferenceExpression> argumentExpressions;

		public MethodCallExpression(string methodName, params ReferenceExpression[] argumentExpressions)
		{
			this.methodName = methodName;
			this.argumentExpressions = new List<ReferenceExpression>(argumentExpressions);
		}

		public override string GenerateCode(CodeGenerator generator)
		{
			StringBuilder code = new StringBuilder();
			code.Append(methodName);
			code.Append("(");
			for (int i = 0; i < argumentExpressions.Count; i++)
			{
				if (i > 0)
				{
					code.Append(",");
				}
				code.Append(argumentExpressions[i].GenerateCode(generator));
			}
			code.Append(")");
			return code.ToString();
		}

		//public override string GenerateCodeInline(CodeGenerator generator)
		//{
		//    // No () needed.
		//    return GenerateCode(generator);
		//}
	}
}
