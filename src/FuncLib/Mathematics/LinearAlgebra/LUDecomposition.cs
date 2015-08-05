// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

using FuncLib.Mathematics.LinearAlgebra.AlgLib;

namespace FuncLib.Mathematics.LinearAlgebra
{
	/// <summary>
	/// LU decomposition.
	/// </summary>
	public sealed class LUDecomposition
	{
		private int m, n, q;
		private double[,] a;
		private int[] pivots;
		private Matrix lowerTriangular, upperTriangular, permutation;

		private LUDecomposition(Matrix matrix)
		{
			m = matrix.Rows;
			n = matrix.Columns;
			q = Math.Min(m, n);

			a = matrix.ToArray();
			pivots = new int[0];

			trfac.rmatrixlu(ref a, m, n, ref pivots);
		}

		public static LUDecomposition Decompose(Matrix matrix)
		{
			return new LUDecomposition(matrix);
		}

		/*public Vector Solve(Vector rightVector)
		{
			if (m != n)
			{
				throw new ArgumentException("The matrix is not a square matrix.");
			}

			if (rightVector.Length != n)
			{
				throw new ArgumentException("The length of the vector is invalid.");
			}

			double[] b = rightVector.ToArray();
			double[] x = null;

			int info = 0;
			densesolver.densesolverreport rep = new densesolver.densesolverreport();

			densesolver.rmatrixlusolve(ref a, ref pivots, n, ref b, ref info, ref rep, ref x);
			if (info != 1)
			{
				throw new ArithmeticException("The matrix is singular.");
			}

			return new Vector(x);
		}*/

		/// <summary>
		/// Matrix determinant using LU decomposition.
		/// </summary>
		public double Determinant()
		{
			if (m != n)
			{
				throw new ArgumentException("The matrix is not a square matrix.");
			}

			return matdet.rmatrixludet(ref a, ref pivots, n);
		}

		/// <summary>
		/// Matrix inversion using LU decomposition.
		/// </summary>
		public bool TryInverse(out Matrix inverse)
		{
			if (n == m)
			{
				double[,] x = (double[,])a.Clone();
				int info = 0;
				matinv.matinvreport rep = new matinv.matinvreport();

				// Compute the inverse.
				matinv.rmatrixluinverse(ref x, ref pivots, n, ref info, ref rep);
				if (info == 1)
				{
					inverse = new Matrix(x);
					return true;
				}
			}

			inverse = null;
			return false;
		}

		/// <summary>
		/// Matrix inversion using LU decomposition.
		/// </summary>
		public Matrix Inverse()
		{
			Matrix inverse;
			if (!TryInverse(out inverse))
			{
				throw new ArithmeticException("The matrix is singular.");
			}

			return inverse;
		}

		/// <summary>
		/// Matrix determinant using LU decomposition.
		/// </summary>
		public static double Determinant(Matrix matrix)
		{
			return new LUDecomposition(matrix).Determinant();
		}

		/// <summary>
		/// Matrix inversion using LU decomposition.
		/// </summary>
		public static bool TryInverse(Matrix matrix, out Matrix inverse)
		{
			return new LUDecomposition(matrix).TryInverse(out inverse);
		}

		/// <summary>
		/// Matrix inversion using LU decomposition.
		/// </summary>
		public static Matrix Inverse(Matrix matrix)
		{
			return new LUDecomposition(matrix).Inverse();
		}

		/// <summary>
		/// The lower unitriangular matrix.
		/// </summary>
		public Matrix LowerTriangular
		{
			get
			{
				if (lowerTriangular == null)
				{
					double[,] al = new double[m, q];
					for (int j = 0; j < q; j++)
					{
						al[j, j] = 1.0;
						for (int i = j + 1; i < m; i++)
						{
							al[i, j] = a[i, j];
						}
					}

					lowerTriangular = new Matrix(al);
				}

				return lowerTriangular;
			}
		}

		/// <summary>
		/// The upper triangular matrix.
		/// </summary>
		public Matrix UpperTriangular
		{
			get
			{
				if (upperTriangular == null)
				{
					double[,] au = new double[q, n];
					for (int i = 0; i < q; i++)
					{
						for (int j = i; j < n; j++)
						{
							au[i, j] = a[i, j];
						}
					}

					upperTriangular = new Matrix(au);
				}

				return upperTriangular;
			}
		}

		/// <summary>
		/// The permutation matrix.
		/// </summary>
		public Matrix Permutation
		{
			get
			{
				if (permutation == null)
				{
					// Transform to a more natural representation.
					int[] p0 = new int[m];
					for (int i = 0; i < m; i++)
					{
						p0[i] = i;
					}
					for (int i = q - 1; i >= 0; i--)
					{
						int j = pivots[i];
						int r = p0[i];
						p0[i] = p0[j];
						p0[j] = r;
					}

					// Compute the permutation matrix.
					double[,] p = new double[m, m];
					for (int i = 0; i < m; i++)
					{
						p[i, p0[i]] = 1.0;
					}

					permutation = new Matrix(p);
				}

				return permutation;
			}
		}
	}
}
