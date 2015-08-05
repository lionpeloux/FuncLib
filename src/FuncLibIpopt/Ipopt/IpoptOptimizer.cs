// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Linq;
using FuncLib.Functions;
using FuncLib.Optimization.Ipopt.Cureos;

namespace FuncLib.Optimization.Ipopt
{
	/// <summary>
	/// Ipopt (Interior Point OPTimizer) non-linear optimizer. Supports sparse structure of gradients and Hessians. The sparsity structure
	/// is determined automatically (to some extend) by the abstract Function framework or by deriving from the Function class manually
	/// and overriding the IsZero property. Consider using <see cref="SparseFunction" /> to express knowledge about the sparsity
	/// structure explicitly.
	/// </summary>
	[Serializable]
	public class IpoptOptimizer : ConstrainedOptimizer
	{

        #region FIELDS

		private IpoptOptionCollection options;
        private bool is_quasinewton_enable;
        
        #endregion

        #region PROPERTIES

        /// <summary>
        /// Any Ipopt option given by its string name and a value of type int, double, or string. Some of the most
        /// common options are exposed as properties in the <see cref="IpoptOptimizer" /> class.
        /// </summary>
        public IpoptOptionCollection Options
        {
            get { return options; }
        }

        /// <summary>
        /// Maximum number of iterations.
        /// </summary>
        public int MaxIterations
        {
            get { return (int)options["max_iter"].Value; }
            set { options["max_iter"].Value = value; }
        }

        /// <summary>
        /// Maximum number of CPU seconds.
        /// </summary>
        public double MaxCpuSeconds
        {
            get { return (double)options["max_cpu_time"].Value; }
            set { options["max_cpu_time"].Value = value; }
        }

        /// <summary>
        /// Ipopt output print level. Setting it to 0 disables output; 5 is the default Ipopt value.
        /// </summary>
        public int PrintLevel
        {
            get { return (int)options["print_level"].Value; }
            set { options["print_level"].Value = value; }
        }

        /// <summary>
        /// Ipopt Quasi-Newton Approximation. True .
        /// </summary>
        public bool EnableQuasiNewton
        {
            get { return is_quasinewton_enable; }
            set
            {
                if (value) options["hessian_approximation"].Value = "limited-memory";
                else options["hessian_approximation"].Value = "exact";

                is_quasinewton_enable = value;
            }
        }
        
        #endregion

        #region CONSTRUCTORS
        
        public IpoptOptimizer()
		{
			options = new IpoptOptionCollection();

			// A few predefined Ipopt options exposed as properties (default Ipopt values, except print_level).
            options.Add(IpoptOptionManager.max_iter, 3000);
            options.Add(IpoptOptionManager.max_cpu_time, 1.0e6);
            options.Add(IpoptOptionManager.print_level, 0);
            options.Add(IpoptOptionManager.hessian_approximation);
            EnableQuasiNewton = false;
		}

        #endregion

        #region METHODS

        /// <summary>
		/// Run the optimizer with the specified initial values. Returns additional Ipopt related convergence status.
		/// </summary>
		public IpoptOptimizerResult RunIpopt(IPoint initialPoint)
		{
			return (IpoptOptimizerResult)Run(initialPoint);
		}

		/// <summary>
		/// Run the optimizer with the specified initial values. Returns additional Ipopt related convergence status.
		/// </summary>
		public IpoptOptimizerResult RunIpopt(params VariableAssignment[] initialAssignments)
		{
			return RunIpopt(new Point(initialAssignments));
		}

		/// <summary>
		/// Prepare the Ipopt optimizer for multiple runs. If used in a multithreading setting, Prepare must be called for each thread.
		/// </summary>
		public override PreparedOptimizer Prepare()
		{
			return new InnerIpoptOptimizer(this);
		}

        #endregion

        #region EVENT 

        /// <summary>
		/// This event is raised once per iteration, providing the user with some information on the state
		/// of the optimization. It also gives the user a way to terminate the optimization prematurely.
		/// </summary>
		public event EventHandler<IpoptIntermediateEventArgs> Intermediate;

        #endregion

        #region PRIVATE CLASS

        [Serializable]
		private class InnerIpoptOptimizer : InnerConstrainedOptimizer
		{
			// Ipopt parameters (following Ipopt naming convention).
            public int n, m, nele_jac, nele_hess;
            public double[] x_L, x_U, g_L, g_U;

			// FuncLib representation of variables and functions, including all required derivatives.
            public Variable[] variables;
            public Function objectiveFunction;
            public Function[] constraintFunctions, objectiveFunctionGradient;
            private Gradient[] constraintFunctionsGradients;
            private HessianStructure hessianStructure;
            private Hessian objectiveFunctionHessian;
            private Hessian[] constraintFunctionsHessians;

			// Other stuff related to this specific optimizer.
            public List<IpoptOption> options;
            public EventHandler<IpoptIntermediateEventArgs> intermediate;
            public int currentIteration;
            public IEvaluation currentEvaluation;

			public InnerIpoptOptimizer(IpoptOptimizer optimizer)
				: base(optimizer)
			{
				// Number of variables and non-linear constraints.
				n = Variables.Count;
				m = FunctionConstraints.Count;

				// Number of nonzeros in the Jacobian of the constraints (incremented below).
				nele_jac = 0;

				// Number of nonzeros in the lower part of the Hessian of the Lagrangian (also incremented below).
				nele_hess = 0;

				// Variables.
				variables = Variables.ToArray();
				// Set the values for the variable bounds.
				x_L = new double[n];
				x_U = new double[n];
				for (int i = 0; i < n; i++)
				{
					VariableConstraint constraint = VariableConstraints.Find(delegate(VariableConstraint c) { return c.Variable == variables[i]; });

					// Use special Ipopt representation of infinity if no constraint on this variable is defined.
					x_L[i] = constraint != null && !double.IsNegativeInfinity(constraint.MinValue) ? constraint.MinValue : IpoptProblem.NegativeInfinity;
					x_U[i] = constraint != null && !double.IsPositiveInfinity(constraint.MaxValue) ? constraint.MaxValue : IpoptProblem.PositiveInfinity;
				}
				// Keeps track of the Hessian entries, whether used by the objective function or one of the constraint functions.
				hessianStructure = new HessianStructure();

				// Prepare objective function, its gradient, and its Hessian. The gradient is assumed dense (otherwise the problem
				// would be rather strange).
				objectiveFunction = ObjectiveFunction;
				objectiveFunctionGradient = new Function[n];
				objectiveFunctionHessian = new Hessian();
				for (int i = 0; i < n; i++)
				{
					objectiveFunctionGradient[i] = objectiveFunction.Derivative(variables[i]);
					for (int j = 0; j <= i; j++)
					{
						Function objectiveFunctionSecondDerivative = objectiveFunctionGradient[i].Derivative(variables[j]);
						if (!objectiveFunctionSecondDerivative.IsZero)
						{
							// Add this entry to the Hessian structure.
							int index = hessianStructure.Add(i, j);

							// Link this specific derivative function the index of that entry.
							objectiveFunctionHessian.Add(index, objectiveFunctionSecondDerivative);
						}
					}
				}
				// Number of nonzeros in the Jacobian of the constraints (incremented below).
				nele_jac = 0;
                
				// Prepare non-linear constraints.
				g_L = new double[m];
				g_U = new double[m];
				constraintFunctions = new Function[m];
				constraintFunctionsGradients = new Gradient[m];
				constraintFunctionsHessians = new Hessian[m];
				for (int i = 0; i < m; i++)
				{
					FunctionConstraint constraint = FunctionConstraints[i];

					// Same as above.
					g_L[i] = !double.IsNegativeInfinity(constraint.MinValue) ? constraint.MinValue : IpoptProblem.NegativeInfinity;
					g_U[i] = !double.IsPositiveInfinity(constraint.MaxValue) ? constraint.MaxValue : IpoptProblem.PositiveInfinity;

					// Constraint function itself, its gradient, and its Hessian.
					constraintFunctions[i] = constraint.Function;
					constraintFunctionsGradients[i] = new Gradient();
					constraintFunctionsHessians[i] = new Hessian();
					for (int j = 0; j < n; j++)
					{
						Function constraintFunctionDerivative = constraintFunctions[i].Derivative(variables[j]);
						if (!constraintFunctionDerivative.IsZero)
						{
							constraintFunctionsGradients[i].Add(j, constraintFunctionDerivative);
							nele_jac++;

							// Fill the lower part of the Hessian, if nonzero.
							for (int k = 0; k <= j; k++)
							{
								Function constraintFunctionSecondDerivative = constraintFunctionDerivative.Derivative(variables[k]);
								if (!constraintFunctionSecondDerivative.IsZero)
								{
									// Add this entry to the Hessian structure, if not already added.
									int index = hessianStructure.Add(j, k);

									// Link this specific derivative function the index of that entry.
									constraintFunctionsHessians[i].Add(index, constraintFunctionSecondDerivative);
								}
							}
						}
					}
				}
				// Number of nonzeros in the lower part of the Hessian of the Lagrangian.
				nele_hess = hessianStructure.Count;

				// Copy and perform checks on Ipopt options.
                options = new List<IpoptOption>();

				// Represent the objective function scaling as an Ipopt option.
				options.Add(IpoptOptionManager.obj_scaling_factor);
                options.Last().Value = ObjectiveFunctionScaling;

                foreach (IpoptOption option in optimizer.Options)
                {
                    if (string.IsNullOrEmpty(option.Name))
                    {
                        throw new OptimizerException("Invalid Ipopt option name.");
                    }

                    if (options.Find(delegate(IpoptOption o) { return o.Name == option.Name; }) != null)
                    {
                        throw new OptimizerException("Duplicated Ipopt option " + option.Name + ".");
                    }

                    // Cloning is needed, since the option object is mutable.
                    options.Add(option.Clone());
                }
				// Copy intermediate callback invocation list.
				intermediate = optimizer.Intermediate;
			}

			public override IOptimizerResult Run(IPoint initialPoint)
			{
				// Doesn't need to call CheckConstraints for Ipopt. That's really cool!
				if (n == 0)
				{
					// No variables to optimize. Just return or Ipopt will fail.
					return new IpoptOptimizerResult(false, false, initialPoint, objectiveFunction.Value(initialPoint), IpoptReturnCode.InvalidProblemDefinition, 0);
				}

				// Reset values from a possible previous run of a prepared optimizer.
				currentIteration = 0;
				currentEvaluation = null;
				// Prepare the initial point.
				double[] x = new double[n];
				for (int i = 0; i < n; i++)
				{
					x[i] = initialPoint[variables[i]];
				}
				double optimalValue;
				IpoptReturnCode returnCode;
				using (IpoptProblem ipopt = new IpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess, eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h))
				{
					DisableWindowCloseEvent();

					// Add options to the Ipopt instance.
					foreach (IpoptOption option in options)
					{
						option.Prepare(ipopt);
					}
					// Requires Ipopt >= 3.9.1. Just uncomment if below this version.
					ipopt.SetIntermediateCallback(Intermediate);
					// This starts solving the problem.
					returnCode = ipopt.SolveProblem(x, out optimalValue, null, null, null, null);
				}
				// Release memory that may be cached in the evaluation object.
				currentEvaluation = null;

				IPoint optimalPoint = CreateOptimalPoint(x);
				bool status;
				bool hasConverged;
				switch (returnCode)
				{
					case IpoptReturnCode.SolveSucceeded:
					case IpoptReturnCode.SolvedToAcceptableLevel:
					case IpoptReturnCode.SearchDirectionBecomesTooSmall:
					case IpoptReturnCode.FeasiblePointFound:
						status = true;
						hasConverged = true;
						break;

					case IpoptReturnCode.MaximumIterationsExceeded:
					case IpoptReturnCode.MaximumCpuTimeExceeded:
						status = true;
						hasConverged = false;
						break;

					default:
						status = false;
						hasConverged = false;
						break;
				}
				return new IpoptOptimizerResult(status, hasConverged, optimalPoint, optimalValue, returnCode, currentIteration);
			}

			[DllImport("Kernel32")]
			private static extern bool SetConsoleCtrlHandler(HandlerRoutine handler, bool add);

			private delegate bool HandlerRoutine(uint ctrlType);

			private static void DisableWindowCloseEvent()
			{
				// Disable annoying "forrtl: error (200): program aborting due to window-CLOSE event" message.
				SetConsoleCtrlHandler(delegate(uint ctrlType) { return true; }, true);
			}

			private bool Intermediate(IpoptAlgorithmMode alg_mod, int iter_count, double obj_value, double inf_pr, double inf_du, double mu, double d_norm, double regularization_size, double alpha_du, double alpha_pr, int ls_trials)
			{
				// Save the total number of iterations.
				currentIteration = iter_count;

				if (intermediate != null)
				{
                    IpoptIntermediateEventArgs e = new IpoptIntermediateEventArgs(alg_mod, iter_count, obj_value, inf_pr, inf_du, mu, d_norm, regularization_size, alpha_du, alpha_pr, ls_trials);
					intermediate(this, e);
					return !e.Cancel;
				}

				return true;
			}

			private bool eval_f(int n, double[] x, bool new_x, out double obj_value)
			{
				if (currentEvaluation == null || new_x)
				{
					// Reuse same evaluation object until instructed to change point. This makes caching (slightly) more efficient.
					currentEvaluation = CreateEvaluation(x);
				}

				try
				{
					obj_value = objectiveFunction.Value(currentEvaluation);
					Validate(obj_value);
				}
				catch (ArithmeticException)
				{
					obj_value = 0.0;
					return false;
				}

				return true;
			}

			private bool eval_grad_f(int n, double[] x, bool new_x, double[] grad_f)
			{
				if (currentEvaluation == null || new_x)
				{
					// Reuse same evaluation object until instructed to change point. This makes caching (slightly) more efficient.
					currentEvaluation = CreateEvaluation(x);
				}

				try
				{
					for (int i = 0; i < n; i++)
					{
						grad_f[i] = objectiveFunctionGradient[i].Value(currentEvaluation);
						Validate(grad_f[i]);
					}
				}
				catch (ArithmeticException)
				{
					return false;
				}

				return true;
			}

			private bool eval_g(int n, double[] x, bool new_x, int m, double[] g)
			{
				if (m == 0)
				{
					// No constraints.
					return true;
				}

				if (currentEvaluation == null || new_x)
				{
					// Reuse same evaluation object until instructed to change point. This makes caching (slightly) more efficient.
					currentEvaluation = CreateEvaluation(x);
				}

				try
				{
					for (int i = 0; i < m; i++)
					{
						g[i] = constraintFunctions[i].Value(currentEvaluation);
						Validate(g[i]);
					}
				}
				catch (ArithmeticException)
				{
					return false;
				}

				return true;
			}

			private bool eval_jac_g(int n, double[] x, bool new_x, int m, int nele_jac, int[] iRow, int[] jCol, double[] values)
			{
				if (m == 0)
				{
					// No constraints.
					return true;
				}

				if (values == null)
				{
					// Set the structure of the Jacobian. May be sparse or dense, depending on the problem.

					for (int i = 0, l = 0; i < m; i++)
					{
						Gradient gradient = constraintFunctionsGradients[i];
						for (int j = 0; j < gradient.Count; j++, l++)
						{
							iRow[l] = i;
							jCol[l] = gradient[j].Index;
						}
					}
				}
				else
				{
					if (currentEvaluation == null || new_x)
					{
						// Reuse same evaluation object until instructed to change point. This makes caching (slightly) more efficient.
						currentEvaluation = CreateEvaluation(x);
					}

					try
					{
						for (int i = 0, l = 0; i < m; i++)
						{
							Gradient gradient = constraintFunctionsGradients[i];
							for (int j = 0; j < gradient.Count; j++, l++)
							{
								values[l] = gradient[j].Function.Value(currentEvaluation);
								Validate(values[l]);
							}
						}
					}
					catch (ArithmeticException)
					{
						return false;
					}
				}

				return true;
			}

			private bool eval_h(int n, double[] x, bool new_x, double obj_factor,
				  int m, double[] lambdas, bool new_lambda,
				  int nele_hess, int[] iRow, int[] jCol,
				  double[] values)
			{
				if (values == null)
				{
					// Set the structure of the Hessian of the Lagrangian. May be sparse or dense, depending on the problem, but
					// fill only the lower part.

					for (int i = 0; i < hessianStructure.Count; i++)
					{
						iRow[i] = hessianStructure[i].Row;
						jCol[i] = hessianStructure[i].Column;
					}
				}
				else
				{
					if (currentEvaluation == null || new_x)
					{
						// Reuse same evaluation object until instructed to change point. This makes caching (slightly) more efficient.
						currentEvaluation = CreateEvaluation(x);
					}

					try
					{
						// Zero out all entries in case they're not used all of them in the objective function.
						for (int i = 0; i < hessianStructure.Count; i++)
						{
							values[i] = 0.0;
						}

						// Scaled Hessian of the objective function.
						foreach (HessianEntry entry in objectiveFunctionHessian)
						{
							double y = entry.Function.Value(currentEvaluation);
							Validate(y);
							values[entry.Index] += obj_factor * y;
						}

						// Add the Hessian of the constraints.
						for (int i = 0; i < m; i++)
						{
							if (lambdas[i] != 0.0)
							{
								foreach (HessianEntry entry in constraintFunctionsHessians[i])
								{
									double y = entry.Function.Value(currentEvaluation);
									Validate(y);
									values[entry.Index] += lambdas[i] * y;
								}
							}
						}
					}
					catch (ArithmeticException)
					{
						return false;
					}
				}

				return true;
			}

			private void Validate(double x)
			{
				if (double.IsNaN(x) || double.IsInfinity(x))
				{
					throw new ArithmeticException();
				}
			}
		}

		[Serializable]
		private class GradientEntry
		{
			private int index;
			private Function function;

			public GradientEntry(int index, Function function)
			{
				this.index = index;
				this.function = function;
			}

			public int Index
			{
				get
				{
					return index;
				}
			}

			public Function Function
			{
				get
				{
					return function;
				}
			}
		}

		[Serializable]
		private class Gradient : List<GradientEntry>
		{
			public void Add(int index, Function function)
			{
				Add(new GradientEntry(index, function));
			}
		}

		[Serializable]
		private class HessianStructureEntry
		{
			private int row, column;

			public HessianStructureEntry(int row, int column)
			{
				this.row = row;
				this.column = column;
			}
			public int Row
			{
                get { return row; }
			}
			public int Column
			{
                get { return column; }
			}
		}

		[Serializable]
		private class HessianStructure : List<HessianStructureEntry>
		{
			public int Add(int row, int column)
			{
				HessianStructureEntry entry = Find(row, column);

				if (entry == null)
				{
					// Doesn't exist. Create new entry.
					entry = new HessianStructureEntry(row, column);
					Add(entry);
				}

				return IndexOf(entry);
			}

			public HessianStructureEntry Find(int row, int column)
			{
				return Find(delegate(HessianStructureEntry e) { return e.Row == row && e.Column == column; });
			}
		}

		[Serializable]
		private class HessianEntry
		{
			private int index;
			private Function function;

			public HessianEntry(int index, Function function)
			{
				this.index = index;
				this.function = function;
			}

			public int Index
			{
				get
				{
					return index;
				}
			}

			public Function Function
			{
				get
				{
					return function;
				}
			}
		}

		[Serializable]
		private class Hessian : List<HessianEntry>
		{
			public void Add(int index, Function function)
			{
				Add(new HessianEntry(index, function));
			}
        }

        #endregion
    }
}
