// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

namespace FuncLib.Functions
{
	/// <summary>
	/// Like <see cref="IPoint" />, but with some mechanism for caching intermediate function values.
	/// </summary>
    public interface IEvaluation : IPoint
    {
		/// <summary>
		/// Try to retrieve a function value from the cache.
		/// </summary>
		bool TryGetValue(Function function, out double value);

		//bool TryGetValue<T>(Function function, out T value);

		/// <summary>
		/// Add a function value to the cache (if available).
		/// </summary>
		void AddValue(Function function, double value);

		//void AddValue<T>(Function function, T value);
    }
}
