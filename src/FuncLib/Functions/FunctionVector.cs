// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

using FuncLib.Mathematics;

namespace FuncLib.Functions
{
	[Serializable]
	public class FunctionVector
	{
		private FunctionMatrix matrix;

		private FunctionVector(FunctionMatrix matrix)
		{
			this.matrix = matrix;
		}

		public FunctionVector(Function[] entries)
		{
			matrix = new FunctionMatrix(entries);
		}

		public Function[] ToArray()
		{
			int n = Length;

			Function[] entries = new Function[n];
			for (int i = 0; i < n; i++)
			{
				entries[i] = this[i];
			}

			return entries;
		}

		public FunctionVector SetValue(int index, Function f)
		{
			return new FunctionVector(matrix.SetValue(index, 0, f));
		}

		public static FunctionVector Zero(int length)
		{
			return new FunctionVector(FunctionMatrix.Zero(length, 1));
		}

		public static FunctionVector Basis(int length, int index)
		{
			return new FunctionVector(FunctionMatrix.Basis(length, 1, index, 0));
		}

		public static FunctionVector operator +(FunctionVector a, FunctionVector b)
		{
			return new FunctionVector(a.matrix + b.matrix);
		}

		public static FunctionVector operator +(FunctionVector a, Function f)
		{
			return new FunctionVector(a.matrix + f);
		}

		public static FunctionVector operator +(Function f, FunctionVector a)
		{
			return a + f;
		}

		public static FunctionVector operator -(FunctionVector a, FunctionVector b)
		{
			return new FunctionVector(a.matrix - b.matrix);
		}

		public static FunctionVector operator -(FunctionVector a, Function f)
		{
			return a + (-f);
		}

		public static FunctionVector operator -(Function f, FunctionVector a)
		{
			return new FunctionVector(f - a.matrix);
		}

		public static FunctionVector operator -(FunctionVector a)
		{
			return new FunctionVector(-a.matrix);
		}

		public static FunctionVector operator *(FunctionMatrix a, FunctionVector b)
		{
			return new FunctionVector(a * b.matrix);
		}

		public static FunctionVector operator *(FunctionVector a, Function f)
		{
			return new FunctionVector(a.matrix * f);
		}

		public static FunctionVector operator *(Function f, FunctionVector a)
		{
			return a * f;
		}

		public static FunctionVector operator /(FunctionVector a, Function f)
		{
			return new FunctionVector(a.matrix / f);
		}

		public static implicit operator FunctionMatrix(FunctionVector a)
		{
			return a.matrix;
		}

		public static explicit operator FunctionVector(FunctionMatrix a)
		{
			if (a.Columns != 1)
			{
				throw new InvalidCastException("Must be a row matrix.");
			}

			return new FunctionVector(a);
		}

		public static implicit operator FunctionVector(Vector a)
		{
			return new FunctionVector((Matrix)a);
		}

		public Function this[int index]
		{
			get
			{
				return matrix[index, 0];
			}
		}

		public int Length
		{
			get
			{
				return matrix.Rows;
			}
		}
	}
}
