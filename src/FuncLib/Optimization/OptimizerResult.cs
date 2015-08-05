// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

using FuncLib.Functions;

namespace FuncLib.Optimization
{
	[Serializable]
	public class OptimizerResult : IOptimizerResult
	{
		private bool status;
		private IPoint optimalPoint;
		private double optimalValue;

		public OptimizerResult(bool status, IPoint optimalPoint, double optimalValue)
		{
			this.status = status;
			this.optimalPoint = optimalPoint;
			this.optimalValue = optimalValue;
		}

		public bool Status
		{
			get
			{
				return status;
			}
		}

		public IPoint OptimalPoint
		{
			get
			{
				return optimalPoint;
			}
		}

		public double OptimalValue
		{
			get
			{
				return optimalValue;
			}
		}
	}
}
