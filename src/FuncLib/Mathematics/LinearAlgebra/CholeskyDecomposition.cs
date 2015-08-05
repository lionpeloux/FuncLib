// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

using FuncLib.Mathematics.LinearAlgebra.AlgLib;

namespace FuncLib.Mathematics.LinearAlgebra
{
	/// <summary>
	/// Cholesky decomposition of a symmetric and positive definite matrix.
	/// </summary>
	public sealed class CholeskyDecomposition // Extend this class like LUDecomposition
	{
		private CholeskyDecomposition()
		{
		}

		public static void InverseDeterminant(Matrix matrix, out Matrix inverse, out double determinant)
		{
			int n = matrix.Rows;

			if (matrix.Columns != n)
			{
				throw new ArgumentException("The matrix isn't a square matrix.");
			}

			double[,] a = matrix.ToArray();

			if (!trfac.spdmatrixcholesky(ref a, n, false))
			{
				throw new ArithmeticException();
			}

			determinant = matdet.spdmatrixcholeskydet(ref a, n);

			int info = 0;
			matinv.matinvreport rep = new matinv.matinvreport();
			matinv.spdmatrixcholeskyinverse(ref a, n, false, ref info, ref rep);

			for (int i = 0; i < n; i++)
			{
				for (int j = i + 1; j < n; j++)
				{
					a[i, j] = a[j, i];
				}
			}

			inverse = new Matrix(a);
		}
	}
}
