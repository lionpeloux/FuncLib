// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

using FuncLib.Mathematics.LinearAlgebra.AlgLib;

namespace FuncLib.Mathematics.LinearAlgebra
{
	/// <summary>
	/// Eigenvalue decomposition of a symmetric matrix.
	/// </summary>
	public sealed class EigenvalueDecomposition
	{
		private Vector eigenvalues;
		private Matrix eigenvectors;

		private EigenvalueDecomposition(Vector eigenvalues, Matrix eigenvectors)
		{
			this.eigenvalues = eigenvalues;
			this.eigenvectors = eigenvectors;
		}

		public static EigenvalueDecomposition Decompose(Matrix matrix)
		{
			Vector eigenvalues;
			Matrix eigenvectors;
			Decompose(matrix, out eigenvalues, out eigenvectors);

			return new EigenvalueDecomposition(eigenvalues, eigenvectors);
		}

		public static void Decompose(Matrix matrix, out Vector eigenvalues, out Matrix eigenvectors)
		{
			int n = matrix.Rows;

			if (matrix.Columns != n)
			{
				throw new ArgumentException("The matrix isn't a square matrix.");
			}

			double[,] a = new double[n, n];
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j <= i; j++)
				{
					if (matrix[i, j] != matrix[j, i])
					{
						throw new ArithmeticException("The matrix must be symmetric.");
					}

					a[i, j] = matrix[i, j];
				}
			}

			double[] d = new double[n];
			double[,] z = new double[n, n];

			if (!evd.smatrixevd(a, n, 1, true, ref d, ref z))
			{
				// The algorithm hasn't converged (rare case).
				throw new ArithmeticException("The algorithm hasn't converged.");
			}

			eigenvalues = new Vector(d);
			eigenvectors = new Matrix(z);
		}

		/// <summary>
		/// Eigenvalues in ascending order.
		/// </summary>
		public Vector Eigenvalues
		{
			get
			{
				return eigenvalues;
			}
		}

		/// <summary>
		/// Eigenvectors stored as columns.
		/// </summary>
		public Matrix Eigenvectors
		{
			get
			{
				return eigenvectors;
			}
		}
	}
}
