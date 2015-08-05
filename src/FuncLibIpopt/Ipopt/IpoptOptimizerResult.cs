// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

using FuncLib.Functions;

namespace FuncLib.Optimization.Ipopt
{
	[Serializable]
	public class IpoptOptimizerResult : OptimizerResult
	{
		private bool hasConverged;
		private IpoptReturnCode returnCode;
		private int iterations;

		public IpoptOptimizerResult(bool status, bool hasConverged, IPoint optimalPoint, double optimalValue, IpoptReturnCode returnCode, int iterations)
			: base(status, optimalPoint, optimalValue)
		{
			this.hasConverged = hasConverged;
			this.returnCode = returnCode;
			this.iterations = iterations;
		}

		/// <summary>
		/// Like <see cref="Status" />, but is only true if the convergence criteria has been met.
		/// </summary>
		public bool HasConverged
		{
			get
			{
				return hasConverged;
			}
		}

		public IpoptReturnCode ReturnCode
		{
			get
			{
				return returnCode;
			}
		}

		/// <summary>
		/// Number of iterations performed.
		/// </summary>
		public int Iterations
		{
			get
			{
				return iterations;
			}
		}
	}
}
