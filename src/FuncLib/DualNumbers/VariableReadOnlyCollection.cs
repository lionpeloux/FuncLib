// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;

using FuncLib.Functions;

namespace FuncLib.DualNumbers
{
	public class VariableReadOnlyCollection : ReadOnlyCollection<Variable>
	{
		public VariableReadOnlyCollection(Variable[] variables)
			: base(CreateList(variables))
		{
		}

		public Variable[] ToArray()
		{
			return new List<Variable>(this).ToArray();
		}

		private static List<Variable> CreateList(Variable[] variables)
		{
			List<Variable> list = new List<Variable>();

			foreach (Variable variable in variables)
			{
				if (variable == null)
				{
					throw new NullReferenceException("Null reference in list of variables.");
				}

				if (list.Contains(variable))
				{
					throw new ArgumentException("Duplicated variable.");
				}

				list.Add(variable);
			}

			return list;
		}
	}
}
