// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

namespace FuncLib.Functions
{
	public interface IPartialValueEvaluation : IPoint
	{
		bool TryGetPartialValue(Function function, out Function partialValue);

		void AddPartialValue(Function function, Function partialValue);
	}
}
