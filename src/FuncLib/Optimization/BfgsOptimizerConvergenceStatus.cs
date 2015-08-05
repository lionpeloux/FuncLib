// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

namespace FuncLib.Optimization
{
	[Serializable]
	public enum BfgsOptimizerConvergenceStatus
	{
		Unknown = 0,
		FunctionDecreaseLessThanEpsF = 1,
		StepLessThanEpsX = 2,
		GradientNormLessThanEpsG = 4,
		MaxIterationsExceeded = 5
	}
}
