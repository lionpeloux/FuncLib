// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

using FuncLib.Functions;

namespace FuncLib.Optimization
{
	[Serializable]
	public class BfgsOptimizerResult : OptimizerResult
	{
		private BfgsOptimizerConvergenceStatus convergenceStatus;

		public BfgsOptimizerResult(bool status, IPoint optimalPoint, double optimalValue, BfgsOptimizerConvergenceStatus convergenceStatus)
			: base(status, optimalPoint, optimalValue)
		{
			this.convergenceStatus = convergenceStatus;
		}

		public BfgsOptimizerConvergenceStatus ConvergenceStatus
		{
			get
			{
				return convergenceStatus;
			}
		}
	}
}
