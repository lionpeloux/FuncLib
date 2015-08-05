// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;

namespace FuncLib.Functions
{
	/// <summary>
	/// Generic thread-safe class for caching any type depending on <see cref="IEvaluation" />.
	/// </summary>
	/// <typeparam name="T">The type to cache, e.g. double or Matrix.</typeparam>
	public class EvaluationCache<T>
	{
		private static int defaultSize;
		private int size;
		private Dictionary<IEvaluation, T> values;
		private LinkedList<IEvaluation> evaluations;
		private object cacheLock;

		static EvaluationCache()
		{
			// Read from configuration file instead?
			defaultSize = 100;
		}

		/// <summary>
		/// New instance of EvaluationCache using default cache size as specified by <see cref="DefaultSize" />.
		/// </summary>
		public EvaluationCache()
			: this(defaultSize)
		{
		}

		/// <summary>
		/// New instance of EvaluationCache with a custom cache size.
		/// </summary>
		/// <param name="size">The size of the rolling cache.</param>
		public EvaluationCache(int size)
		{
			this.size = size;

			values = new Dictionary<IEvaluation, T>(size);
			evaluations = new LinkedList<IEvaluation>();
			cacheLock = new object();
		}

		public void Add(IEvaluation evaluation, T value)
		{
			if (size == 0)
			{
				// No caching; avoid locking.
				return;
			}

			lock (cacheLock)
			{
				if (evaluations.Count >= size)
				{
					values.Remove(evaluations.First.Value);
					evaluations.RemoveFirst();
				}

				values.Add(evaluation, value);
				evaluations.AddLast(evaluation);
			}
		}

		public bool TryGetValue(IEvaluation evaluation, out T value)
		{
			if (size == 0)
			{
				value = default(T);
				return false;
			}

			lock (cacheLock)
			{
				return values.TryGetValue(evaluation, out value);
			}
		}

		public static int DefaultSize
		{
			get
			{
				return defaultSize;
			}
			set
			{
				defaultSize = value;
			}
		}
	}
}
