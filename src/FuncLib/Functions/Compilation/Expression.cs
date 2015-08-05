// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Globalization;
using System.Collections.Generic;
using System.Text;
using System.Text.RegularExpressions;

namespace FuncLib.Functions.Compilation
{
	public abstract class Expression
	{
		public abstract string GenerateCode(CodeGenerator generator);

		//public virtual string GenerateCodeInline(CodeGenerator generator)
		//{
		//    // When using the code snippet inlined in other code snippets, the default behavior is to add ().
		//    return "(" + GenerateCode(generator) + ")";
		//}

		public static implicit operator Expression(double value)
		{
			return new DoubleExpression(value);
		}

		public static implicit operator Expression(string expression)
		{
			return new StringExpression(expression);
		}

		public static Expression operator +(Expression leftExpression, Expression rightExpression)
		{
			JoinedExpression expression;
			if (leftExpression is JoinedExpression)
			{
				expression = new JoinedExpression((JoinedExpression)leftExpression, rightExpression);
			}
			else
			{
				expression = new JoinedExpression(leftExpression, rightExpression);
			}

			// Return an expression that's optimized if possible.
			return expression.Optimize();
		}

		private class StringExpression : Expression
		{
			private string expression;

			public StringExpression(string expression)
			{
				this.expression = expression;
			}

			public override string GenerateCode(CodeGenerator generator)
			{
				return expression;
			}

			public string Expression
			{
				get
				{
					return expression;
				}
			}
		}

		private class JoinedExpression : Expression
		{
			private static Regex mathMethod;
			private List<Expression> innerExpressions;

			static JoinedExpression()
			{
				// Used for identifying expressions of the type "Math.X(..., ..., ...)".
				mathMethod = new Regex(@"^(Math\.[A-Za-z]+[[A-Za-z0-9]*)\($");
			}

			public JoinedExpression(params Expression[] innerExpressions)
			{
				this.innerExpressions = new List<Expression>(innerExpressions);
			}

			public JoinedExpression(JoinedExpression leftExpressions, Expression rightExpression)
			{
				innerExpressions = new List<Expression>(leftExpressions.innerExpressions);
				innerExpressions.Add(rightExpression);
			}

			public override string GenerateCode(CodeGenerator generator)
			{
				StringBuilder code = new StringBuilder();
				foreach (Expression expression in innerExpressions)
				{
					code.Append(expression.GenerateCode(generator));
				}
				return code.ToString();
			}

			public Expression Optimize()
			{
				if (innerExpressions.Count >= 3)
				{
					StringExpression e1 = innerExpressions[0] as StringExpression;
					StringExpression e2 = innerExpressions[innerExpressions.Count - 1] as StringExpression;

					if (e1 != null && e2 != null && e2.Expression == ")")
					{
						Match m = mathMethod.Match(e1.Expression);
						if (m.Success)
						{
							List<ReferenceExpression> argumentExpressions = new List<ReferenceExpression>();
							for (int i = 1; i < innerExpressions.Count - 1; i++)
							{
								ReferenceExpression argumentExpression = innerExpressions[i] as ReferenceExpression;
								if (argumentExpression == null)
								{
									// No optimization is possible.
									return this;
								}
								argumentExpressions.Add(argumentExpression);
							}
							return new MethodCallExpression(m.Groups[1].Value, argumentExpressions.ToArray());
						}
					}
				}

				return this;
			}
		}

		private class DoubleExpression : Expression
		{
			private double value;

			public DoubleExpression(double value)
			{
				this.value = value;
			}

			public override string GenerateCode(CodeGenerator generator)
			{
				// Important to add double literal.
				return value.ToString("r", CultureInfo.InvariantCulture) + "d";
			}

			//public override string GenerateCodeInline(CodeGenerator generator)
			//{
			//    if (value >= 0.0)
			//    {
			//        // Skip () for positive numbers.
			//        return GenerateCode(generator);
			//    }
			//
			//    return base.GenerateCodeInline(generator);
			//}
		}
	}
}
