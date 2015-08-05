// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.ComponentModel;

namespace FuncLib.Optimization.Ipopt
{
	public class IpoptIntermediateEventArgs : CancelEventArgs
	{
		public IpoptAlgorithmMode algorithmMode;
        public int iteration;
        public double objectiveValue;
        public double primalInfeasibility;
        public double dualInfeasibility;
        public double barrierParameter;
        public double primalStep;
        public double regularizationTerm;
        public double dualStepSize;
        public double primalStepSize;
        public int lineSearchStepsCount;

        public double[] x;
        public double[] eval_f;
        public double[] eval_g;
        public double[] eval_grad_f;
        public double[,] eval_jac_g;

        public IpoptIntermediateEventArgs(IpoptAlgorithmMode algorithmMode, int iteration, double objectiveValue)
		{
			this.algorithmMode = algorithmMode;
			this.iteration = iteration;
			this.objectiveValue = objectiveValue;
		}
        public IpoptIntermediateEventArgs(IpoptAlgorithmMode algorithmMode, int iteration, double objectiveValue, double primalInfeasibility, double dualInfeasibility, double barrierParameter, double primalStep, double regularizationTerm, double dualStepSize, double primalStepSize, int lineSearchStepsCount)
        {
            this.algorithmMode = algorithmMode;
			this.iteration = iteration;
			this.objectiveValue = objectiveValue;
            this.primalInfeasibility = primalInfeasibility;
            this.dualInfeasibility = dualInfeasibility;
            this.barrierParameter = barrierParameter;
            this.primalStep = primalStep;
            this.regularizationTerm = regularizationTerm;
            this.dualStepSize = dualStepSize;
            this.primalStepSize = primalStepSize;
            this.lineSearchStepsCount = lineSearchStepsCount;
        }


		/// <summary>
		/// Current Ipopt algorithm mode.
		/// </summary>
		public IpoptAlgorithmMode AlgorithmMode
		{
			get
			{
				return algorithmMode;
			}
		}

		/// <summary>
		/// Current iteration number (iteration).
		/// </summary>
		public int Iteration
		{
			get
			{
				return iteration;
			}
		}

		/// <summary>
		/// The unscaled objective value at the current point (objective).
		/// </summary>
		public double ObjectiveValue
		{
			get
			{
				return objectiveValue;
			}
		}

        /// <summary>
		/// The unscaled constraint violation at the current point (inf_pr).
		/// </summary>
		public double PrimalInfeasibility
		{
			get
			{
				return primalInfeasibility;
			}
		}

        /// <summary>
		/// The scaled dual infeasibility at the current point (inf_du).
		/// </summary>
		public double DualInfeasibility
		{
			get
			{
				return dualInfeasibility;
			}
		}

        /// <summary>
        /// Log_10 of the value of the barrier parameter (mu).
		/// </summary>
		public double BarrierParameter
		{
			get
			{
				return barrierParameter;
			}
		}

        /// <summary>
        /// The infinity norm (max) of the primal step (for the original variables x and the internal slack variables s).
        /// </summary>
        public double PrimalStep
        {
            get
            {
                return primalStep;
            }
        }

        /// <summary>
        /// Log_10 of the value of the regularization term for the Hessian of the Lagrangian in the augmented system.
        /// </summary>
        public double RegularizationTerm
        {
            get
            {
                return regularizationTerm;
            }
        }
        
        /// <summary>
        /// The stepsize for the dual variables.
        /// </summary>
        public double DualStepSize
        {
            get
            {
                return dualStepSize;
            }
        }

        /// <summary>
        /// The stepsize for the primal variables.
        /// </summary>
        public double PriamlStepSize
        {
            get
            {
                return primalStepSize;
            }
        }

        /// <summary>
        /// The number of backtracking line search steps
        /// </summary>
        public int LineSearchStepsCount
        {
            get
            {
                return lineSearchStepsCount;
            }
        }       
	}
}
