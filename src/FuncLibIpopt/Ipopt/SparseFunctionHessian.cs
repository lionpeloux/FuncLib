// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;

using FuncLib.Functions;

namespace FuncLib.Optimization.Ipopt
{
	[Serializable]
	public class SparseFunctionHessian : List<Tuple<Variable, Variable>>
	{
		public SparseFunctionHessian()
		{
		}

		public SparseFunctionHessian(params Tuple<Variable, Variable>[] entries)
			: base(entries)
		{
		}

		private SparseFunctionHessian(SparseFunctionHessian hessian)
			: base(hessian)
		{
		}

		public void Add(Variable variable1, Variable variable2)
		{
			Add(new Tuple<Variable, Variable>(variable1, variable2));
		}

		public bool Contains(Variable variable1, Variable variable2)
		{
			return Find(delegate(Tuple<Variable, Variable> entry) { return entry.Item1 == variable1 && entry.Item2 == variable2 || entry.Item1 == variable2 && entry.Item2 == variable1; }) != null;
		}

		public SparseFunctionHessian Clone()
		{
			return new SparseFunctionHessian(this);
		}
	}

	/*[Serializable]
	public class Tuple<T1, T2>
	{
		private T1 item1;
		private T2 item2;

		public Tuple(T1 item1, T2 item2)
		{
			this.item1 = item1;
			this.item2 = item2;
		}

		public T1 Item1
		{
			get
			{
				return item1;
			}
		}

		public T2 Item2
		{
			get
			{
				return item2;
			}
		}
	}*/
}
