// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

namespace FuncLib.Optimization.Ipopt
{
	/// <summary>
	/// When using <see cref="IpoptIntermediateEventArgs" /> this indicates which mode the algorithm is in.
	/// </summary>
	public enum IpoptAlgorithmMode
	{
		RegularMode = 0,
		RestorationPhaseMode = 1
	}
}
