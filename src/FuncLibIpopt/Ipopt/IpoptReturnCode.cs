// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

namespace FuncLib.Optimization.Ipopt
{
	/// <summary>
	/// Ipopt optimizer return codes.
	/// </summary>
	public enum IpoptReturnCode
	{
        SolveSucceeded = 0,
        SolvedToAcceptableLevel = 1,
        InfeasibleProblemDetected = 2,
        SearchDirectionBecomesTooSmall = 3,
        DivergingIterates = 4,
        UserRequestedStop = 5,
        FeasiblePointFound = 6,

        MaximumIterationsExceeded = -1,
        RestorationFailed = -2,
        ErrorInStepComputation = -3,
        MaximumCpuTimeExceeded = -4,
        NotEnoughDegreesOfFreedom = -10,
        InvalidProblemDefinition = -11,
        InvalidOption = -12,
        InvalidNumberDetected = -13,

        UnrecoverableException = -100,
        NonIpoptExceptionThrown = -101,
        InsufficientMemory = -102,
        InternalError = -199,

        ProblemNotInitialized = -900
    }
}
