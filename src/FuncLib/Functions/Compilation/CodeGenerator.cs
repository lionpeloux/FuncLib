// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;
using System.Text;

namespace FuncLib.Functions.Compilation
{
	/// <summary>
	/// This class is responsible for generating C# code representation of the functions.
	/// </summary>
	public sealed class CodeGenerator
	{
		private List<Function> functions;
		private List<Variable> variables;
		private List<ReferenceExpression> expressions;
		private Dictionary<Function, ReferenceExpression> references;
		private int nextIdentifier;

		public CodeGenerator(Function[] functions, Variable[] variables)
		{
			this.functions = new List<Function>(functions);
			this.variables = new List<Variable>(variables);
		}

		public Expression GetReference(Function function)
		{
			if (function is Variable)
			{
				return GetVariableReference((Variable)function);
			}

			if (references.ContainsKey(function))
			{
				// Calling AddReference causes the expression NOT to be inlined.
				references[function].AddReference();
				return references[function];
			}

			ReferenceExpression expression = new ReferenceExpression(function.Compile(this));
			expressions.Add(expression);

			references.Add(function, expression);
			return references[function];
		}

		public Expression GetVariableReference(Variable variable)
		{
			if (!variables.Contains(variable))
			{
				throw new VariableNotAssignedException(variable, "A variable not added for compilation encountered.");
			}

			return "x[" + variables.IndexOf(variable) + "]";

			/*ReferenceExpression expression = new ReferenceExpression("x[" + variables.IndexOf(variable) + "]");
			expressions.Add(expression);

			references.Add(variable, expression);
			return references[variable];*/
		}

		public string GetNextIdentifier()
		{
			return "a" + nextIdentifier++;
		}

		public string GenerateCode()
		{
			// List of all expressions as returned by the Function class. Wrapped in ReferenceExpression that's
			// responsible for determining whether to inline or not.
			expressions = new List<ReferenceExpression>();

			// Map the generated expressions to the corresponding function.
			references = new Dictionary<Function, ReferenceExpression>();

			// For local C# variables used with non-inlined expressions.
			nextIdentifier = 0;

			// Implicitly adds all functions to the references dictionary. Using this below to extract final values.
			foreach (Function function in functions)
			{
				GetReference(function);
			}

			StringBuilder code = new StringBuilder();

			code.AppendLine("using System;");
			code.AppendLine("public class GeneratedFunction");
			code.AppendLine("{");
			code.AppendLine("public static void Value(double[] x, double[] y)");
			code.AppendLine("{");

			// All intermediate statements (in the right order).
			foreach (ReferenceExpression expression in expressions)
			{
				expression.GenerateCodeStatement(code, this);
			}

			// Extract final values.
			for (int i = 0; i < functions.Count; i++)
			{
				code.AppendLine("y[" + i + "]=" + references[functions[i]].GenerateCode(this) + ";");
			}

			code.AppendLine("}");
			code.AppendLine("}");

			return code.ToString();
		}
	}
}
