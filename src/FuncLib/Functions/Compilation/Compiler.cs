// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;

namespace FuncLib.Functions.Compilation
{
	/// <summary>
	/// Generates fast IL code representing a function or a group of functions, possibly including derivatives. Use
	/// one of the overloads of the only public static method <see cref="Compile" /> to perform the compilation.
	/// </summary>
	public sealed class Compiler
	{
		private List<Variable> variables;
		private List<Function> compilableFunctions;
		private Dictionary<Function, double> constants;
		private Dictionary<Function, Dictionary<Variable, Function>> derivatives;

		private Compiler(Variable[] variables)
		{
			this.variables = new List<Variable>();

			foreach (Variable variable in variables)
			{
				AddVariable(variable);
			}
		}

		/// <summary>
		/// Compiles a function depending on some variables.
		/// </summary>
		public static CompiledFunction Compile(Function function, params Variable[] variables)
		{
			return Compile(function, variables, 0);
		}

		/// <summary>
		/// Compiles a function depending on some variables, including all partial derivatives up to the order specified.
		/// </summary>
		public static CompiledFunction Compile(Function function, Variable[] variables, int order)
		{
			return Compile(new Function[] { function }, variables, order)[0];
		}

		/// <summary>
		/// Compiles a group of functions depending on some variables.
		/// </summary>
		public static CompiledFunction[] Compile(Function[] functions, Variable[] variables)
		{
			return Compile(functions, variables, 0);
		}

		/// <summary>
		/// Compiles a group of functions depending on some variables, including all partial derivatives up to the order specified.
		/// </summary>
		public static CompiledFunction[] Compile(Function[] functions, Variable[] variables, int order)
		{
			return new Compiler(variables).Compile(functions, order);
		}

		/// <summary>
		/// This non-static method is responsible for handling the actual compilation (code generation is done by CodeGenerator, though)
		/// and keeping track of the derivatives.
		/// </summary>
		private CompiledFunction[] Compile(Function[] functions, int order)
		{
			if (order < 0)
			{
				throw new ArgumentOutOfRangeException("order", "The order of derivatives must be non-negative.");
			}

			if (order > 2)
			{
				// Higher order derivatives would be nice, but would probably make the code below significantly more complex.
				throw new NotSupportedException("Only compilation with first and second partial derivatives is supported.");
			}

			// List of all derivatives encountered.
			derivatives = new Dictionary<Function, Dictionary<Variable, Function>>();

			// List of functions that is known to be constants. They're skipped from the compilation and added later.
			constants = new Dictionary<Function, double>();

			// All functions that is to be compiled, i.e. including derivatives.
			compilableFunctions = new List<Function>();

			foreach (Function function in functions)
			{
				PrepareFunction(function, order);
			}

			CompiledFunctionEvaluator evaluator = null;
			if (compilableFunctions.Count > 0)
			{
				// Generate C# code.
				string generatedCode = new CodeGenerator(compilableFunctions.ToArray(), variables.ToArray()).GenerateCode();

				// And compile this C# code to IL.
				evaluator = new CompiledFunctionEvaluator(generatedCode, variables.ToArray(), compilableFunctions.Count);
			}

			List<CompiledFunction> compiledFunctions = new List<CompiledFunction>();
			foreach (Function function in functions)
			{
				compiledFunctions.Add(BuildCompiledFunction(function, order, evaluator));
			}

			return compiledFunctions.ToArray();
		}

		private void PrepareFunction(Function function, int order)
		{
			AddFunction(function);

			if (order >= 1)
			{
				// Add first order partial derivatives.
				for (int i = 0; i < variables.Count; i++)
				{
					Function firstOrderDerivative = function.Derivative(variables[i]);
					AddFunction(firstOrderDerivative);
					AddDerivative(function, variables[i], firstOrderDerivative);

					if (order >= 2)
					{
						// Add second order partial derivatives.
						for (int j = 0; j <= i; j++)
						{
							Function secondOrderDerivative = firstOrderDerivative.Derivative(variables[j]);
							AddFunction(secondOrderDerivative);
							AddDerivative(firstOrderDerivative, variables[j], secondOrderDerivative);

							// Also add the partial derivative with the order of differentiation interchanged.
							AddDerivative(derivatives[function][variables[j]], variables[i], secondOrderDerivative);
						}
					}
				}
			}
		}

		private CompiledFunction BuildCompiledFunction(Function function, int order, CompiledFunctionEvaluator evaluator)
		{
			if (constants.ContainsKey(function))
			{
				// Just return a CompiledFunction representation of a constant and ignore derivatives (they're zero).
				return new ConstantCompiledFunction(constants[function]);
			}

			// Prepare first order derivatives. Null indicates that this order of derivatives wasn't compiled.
			Dictionary<Variable, Function> firstOrderDerivatives = null;

			if (order >= 1)
			{
				// First order partial derivatives.
				firstOrderDerivatives = new Dictionary<Variable, Function>();

				for (int i = 0; i < variables.Count; i++)
				{
					Function firstOrderDerivative = derivatives[function][variables[i]];
					if (constants.ContainsKey(firstOrderDerivative))
					{
						// Since it's a constant function, we don't have to care about second order partial derivatives. Just add as a constant.
						firstOrderDerivatives.Add(variables[i], constants[firstOrderDerivative]);
					}
					else
					{
						// Do the same for second order partial derivatives.
						Dictionary<Variable, Function> secondOrderDerivatives = null;

						if (order >= 2)
						{
							secondOrderDerivatives = new Dictionary<Variable, Function>();

							for (int j = 0; j < variables.Count; j++)
							{
								Function secondOrderDerivative = derivatives[firstOrderDerivative][variables[j]];
								if (constants.ContainsKey(secondOrderDerivative))
								{
									secondOrderDerivatives.Add(variables[j], constants[secondOrderDerivative]);
								}
								else
								{
									secondOrderDerivatives.Add(variables[j], new InnerCompiledFunction(evaluator, compilableFunctions.IndexOf(secondOrderDerivative), null)); // FIXME Use same instance! [i][j]=[j][i]
								}
							}
						}

						// Now that all second order partial derivatives are taken care of, we can create an
						// object for this particular first order partial derivative.
						firstOrderDerivatives.Add(variables[i], new InnerCompiledFunction(evaluator, compilableFunctions.IndexOf(firstOrderDerivative), secondOrderDerivatives));
					}
				}
			}

			return new InnerCompiledFunction(evaluator, compilableFunctions.IndexOf(function), firstOrderDerivatives);
		}

		private void AddVariable(Variable variable)
		{
			if (variable == null)
			{
				throw new ArgumentNullException("Null reference in list of variables.");
			}

			if (variables.Contains(variable))
			{
				throw new ArgumentException("Duplicated variable.");
			}

			variables.Add(variable);
		}

		private void AddFunction(Function function)
		{
			if (function == null)
			{
				throw new NullReferenceException("Null reference in list of functions.");
			}

			if (compilableFunctions.Contains(function))
			{
				// Already added somehow. Just skip it.
				return;
			}

			if (function.IsZero)
			{
				// The function is zero as indicated by the virtual method IsZero. Add it to the list of known constants.
				constants[function] = 0.0;
				return;
			}

			ConstantFunction function0 = function as ConstantFunction;
			if (function0 != null)
			{
				// The function is a non-zero constant. Add it to the list of known constants.
				constants[function] = function0.Constant;
				return;
			}

			compilableFunctions.Add(function);
		}

		private void AddDerivative(Function function, Variable variable, Function derivative)
		{
			if (!derivatives.ContainsKey(function))
			{
				derivatives.Add(function, new Dictionary<Variable, Function>());
			}

			derivatives[function][variable] = derivative;
		}

		[Serializable]
		private class InnerCompiledFunction : CompiledFunction
		{
			private CompiledFunctionEvaluator evaluator;
			private int functionIndex;

			public InnerCompiledFunction(CompiledFunctionEvaluator evaluator, int functionIndex, Dictionary<Variable, Function> derivatives)
				: base(derivatives)
			{
				this.evaluator = evaluator;
				this.functionIndex = functionIndex;
			}

			public override double Value(params double[] x)
			{
				return evaluator.Value(x)[functionIndex];
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				return evaluator.Value(evaluation)[functionIndex];
			}

			public override string GeneratedCode
			{
				get
				{
					return evaluator.GeneratedCode;
				}
			}
		}

		[Serializable]
		private class ConstantCompiledFunction : CompiledFunction
		{
			private double constant;

			public ConstantCompiledFunction(double constant)
				: base(new Dictionary<Variable, Function>())
			{
				this.constant = constant;
			}

			public override double Value(params double[] values)
			{
				return constant;
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				return constant;
			}

			public override bool IsZero
			{
				get
				{
					return constant == 0.0;
				}
			}

			public override string GeneratedCode
			{
				get
				{
					// Ignore code for the degenerated function.
					return null;
				}
			}
		}
	}
}
