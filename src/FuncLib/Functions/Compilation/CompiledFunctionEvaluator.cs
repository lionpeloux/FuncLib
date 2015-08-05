// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.CodeDom.Compiler;
using System.Reflection;
using Microsoft.CSharp;

namespace FuncLib.Functions.Compilation
{
	[Serializable]
	public class CompiledFunctionEvaluator
	{
		private string generatedCode;
		private Variable[] variables;
		private int functionCount;
		private object lazyLock;
		[NonSerialized]
		private MethodInfo methodInfo;
		[NonSerialized]
		private EvaluationCache<double[]> values;

		public CompiledFunctionEvaluator(string generatedCode, Variable[] variables, int functionCount)
		{
			this.generatedCode = generatedCode;
			this.variables = variables;
			this.functionCount = functionCount;

			lazyLock = new object();

			// Don't use lazy compilation. We want errors to be shown immediately.
			Compile();
		}

		public void Compile()
		{
			// Prepare the compiler.
			CompilerParameters parameters = new CompilerParameters();
			parameters.GenerateInMemory = true;
			parameters.TreatWarningsAsErrors = false;
			parameters.GenerateExecutable = false;
			parameters.CompilerOptions = "/optimize";
			parameters.IncludeDebugInformation = false;
			parameters.ReferencedAssemblies.Add("System.dll");

			// Run the compiler.
			CompilerResults results = new CSharpCodeProvider().CompileAssemblyFromSource(parameters, new string[] { generatedCode });

			if (results.Errors.HasErrors)
			{
				// Something's wrong with the code. Just show the first error.
				throw new CompilerException("Compile error: " + results.Errors[0].ToString(), results.Errors);
			}

			// Locate the method we've just defined (cf. CodeGenerator).
			methodInfo = results.CompiledAssembly.GetModules()[0].GetType("GeneratedFunction").GetMethod("Value");
		}

		public double[] Value(IEvaluation evaluation)
		{
			// Make thread-safe even after deserialization. The lazyLock object is serialized. Do this smarter?
			lock (lazyLock)
			{
				if (values == null)
				{
					values = new EvaluationCache<double[]>();
				}
			}

			double[] y;
			if (!values.TryGetValue(evaluation, out y))
			{
				y = Value(Transform(evaluation));
				values.Add(evaluation, y);
			}

			return y;
		}

		public double[] Value(double[] x)
		{
			lock (lazyLock) // FIXME Used twice! Too slow?
			{
				if (methodInfo == null)
				{
					// Compile again after a deserialization.
					Compile();
				}
			}

			double[] y = new double[functionCount];
			methodInfo.Invoke(null, new object[] { x, y });
			return y;
		}

		private double[] Transform(IEvaluation evaluation)
		{
			double[] x = new double[variables.Length];
			for (int i = 0; i < variables.Length; i++)
			{
				x[i] = evaluation[variables[i]];
			}
			return x;
		}

		public string GeneratedCode
		{
			get
			{
				return generatedCode;
			}
		}
	}
}
