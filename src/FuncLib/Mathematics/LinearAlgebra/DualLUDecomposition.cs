// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

using FuncLib.DualNumbers;

namespace FuncLib.Mathematics.LinearAlgebra
{
	/// <summary>
	/// Provides the same static methods as the <see cref="LUDecomposition"/> class but in <see cref="DualNumber "/> versions,
	/// i.e. the inverse and the determinant may be computed with gradients and Hessians. Direct computation is used for low
	/// dimensions (less than three).
	/// </summary>
	public static class DualLUDecomposition
	{
		public static DualMatrix Inverse(DualMatrix matrix)
		{
			int m = matrix.Rows;

			if (matrix.Columns != m)
			{
				throw new ArgumentException("The matrix isn't a square matrix.");
			}

			if (m <= 3)
			{
				DualNumber determinant = DeterminantDirect(m, matrix);
				return InverseDirect(m, matrix, determinant);
			}

			int n = GradientLength(m, matrix);

			return Inverse(matrix, m, n, LUDecomposition.Inverse(matrix.GetValues()));
		}

		public static DualNumber Determinant(DualMatrix matrix)
		{
			int m = matrix.Rows;

			if (matrix.Columns != m)
			{
				throw new ArgumentException("The matrix isn't a square matrix.");
			}

			if (m <= 3)
			{
				return DeterminantDirect(m, matrix);
			}

			int n = GradientLength(m, matrix);

			LUDecomposition lu = LUDecomposition.Decompose(matrix.GetValues());

			return Determinant(matrix, m, n, lu.Inverse(), lu.Determinant());
		}

		/// <summary>
		/// This method saves some time if both the inverse and the determinant is needed by only having to compute
		/// the inverse once.
		/// </summary>
		public static void InverseDeterminant(DualMatrix matrix, out DualMatrix inverse, out DualNumber determinant)
		{
			int m = matrix.Rows;

			if (matrix.Columns != m)
			{
				throw new ArgumentException("The matrix isn't a square matrix.");
			}

			if (m <= 3)
			{
				// The naive implementation seems to be faster (for all values of n). Maybe it's
				// faster to perform the LU decomposition directly on the DualNumbers?
				determinant = DeterminantDirect(m, matrix);
				inverse = InverseDirect(m, matrix, determinant);
				return;
			}

			int n = GradientLength(m, matrix);

			LUDecomposition lu = LUDecomposition.Decompose(matrix.GetValues());
			Matrix inverse0 = lu.Inverse();
			double determinant0 = lu.Determinant();

			inverse = Inverse(matrix, m, n, inverse0);
			determinant = Determinant(matrix, m, n, inverse0, determinant0);
		}


		private static DualMatrix Inverse(DualMatrix matrix, int m, int n, Matrix inverse)
		{
			if (n < 0)
			{
				// Return a matrix without derivatives information (using the implicit matrix conversion).
				return inverse;
			}

			double[,][] gradientArrays = new double[m, m][];

			// Compute all gradients in the first loop. These are used to compute the Hessian.

			for (int i0 = 0; i0 < m; i0++)
			{
				for (int j0 = 0; j0 < m; j0++)
				{
					double[] gradientArray = new double[n];

					for (int i = 0; i < n; i++)
					{
						// Formula (36) in The Matrix Cookbook.

						double s = 0.0;
						for (int k = 0; k < m; k++)
						{
							for (int l = 0; l < m; l++)
							{
								s += inverse[i0, k] * matrix[k, l].Gradient[i] * inverse[l, j0];
							}
						}
						gradientArray[i] = -s;
					}

					gradientArrays[i0, j0] = gradientArray;
				}
			}

			// Compute Hessians and DualNumber instances.

			DualNumber[,] a = new DualNumber[m, m];
			for (int i0 = 0; i0 < m; i0++)
			{
				for (int j0 = 0; j0 < m; j0++)
				{
					double[] hessianArray = new double[n * (n + 1) / 2];

					for (int i = 0, h = 0; i < n; i++)
					{
						for (int j = i; j < n; j++, h++)
						{
							double s = 0.0;
							for (int k = 0; k < m; k++)
							{
								for (int l = 0; l < m; l++)
								{
									s += gradientArrays[i0, k][i] * matrix[k, l].Gradient[j] * inverse[l, j0]
										+ inverse[i0, k] * matrix[k, l].Gradient[h] * inverse[l, j0]
										+ inverse[i0, k] * matrix[k, l].Gradient[j] * gradientArrays[l, j0][i];
								}
							}
							hessianArray[h] = -s;
						}
					}

					a[i0, j0] = new DualNumber(inverse[i0, j0], gradientArrays[i0, j0], hessianArray);
				}
			}

			return new DualMatrix(a);
		}

		private static DualNumber Determinant(DualMatrix matrix, int m, int n, Matrix inverse, double determinant)
		{
			if (n < 0)
			{
				return determinant;
			}

			// Compute the gradient in the first loop (not scaled by the determinant yet).
			// These are used to compute the Hessian.

			double[] gradientArray = new double[n];

			for (int i = 0; i < n; i++)
			{
				// Formula (41) in The Matrix Cookbook.

				double s = 0.0;
				for (int k = 0; k < m; k++)
				{
					for (int l = 0; l < m; l++)
					{
						s += inverse[k, l] * matrix[l, k].Gradient[i];
					}
				}
				gradientArray[i] = s;
			}

			// Compute the Hessian.

			double[] hessianArray = new double[n * (n + 1) / 2];

			for (int i = 0, h = 0; i < n; i++)
			{
				for (int j = i; j < n; j++, h++)
				{
					// Formula (42) in The Matrix Cookbook. Also works when taking partial derivatives
					// with respect to two different variables.

					double s = 0.0;
					for (int k = 0; k < m; k++)
					{
						for (int l = 0; l < m; l++)
						{
							double si = 0.0;
							double sj = 0.0;
							for (int q = 0; q < m; q++)
							{
								si += inverse[k, q] * matrix[q, l].Gradient[i];
								sj += inverse[l, q] * matrix[q, k].Gradient[j];
							}
							s += inverse[k, l] * matrix[l, k].Hessian[i, j] - si * sj;
						}
					}
					hessianArray[h] = determinant * (gradientArray[i] * gradientArray[j] + s);
				}

				// Final scaling of the gradient. This index is not used anymore in the loop.
				gradientArray[i] *= determinant;
			}

			// Finally create the DualNumber instance.
			return new DualNumber(determinant, gradientArray, hessianArray);
		}

		private static DualNumber DeterminantDirect(int m, DualMatrix matrix)
		{
			if (m == 1)
			{
				DualNumber a11 = matrix[0, 0];
				return a11;
			}
			else if (m == 2)
			{
				DualNumber a11 = matrix[0, 0];
				DualNumber a12 = matrix[0, 1];
				DualNumber a21 = matrix[1, 0];
				DualNumber a22 = matrix[1, 1];
				return a11 * a22 - a12 * a21;
			}
			else if (m == 3)
			{
				DualNumber a11 = matrix[0, 0];
				DualNumber a12 = matrix[0, 1];
				DualNumber a13 = matrix[0, 2];
				DualNumber a21 = matrix[1, 0];
				DualNumber a22 = matrix[1, 1];
				DualNumber a23 = matrix[1, 2];
				DualNumber a31 = matrix[2, 0];
				DualNumber a32 = matrix[2, 1];
				DualNumber a33 = matrix[2, 2];
				return a11 * a22 * a33 + a12 * a23 * a31 + a13 * a21 * a32 - a11 * a23 * a32 - a12 * a21 * a33 - a13 * a22 * a31;
			}
			else
			{
				throw new NotImplementedException();
			}
		}

		private static DualMatrix InverseDirect(int m, DualMatrix matrix, DualNumber determinant)
		{
			if (m == 1)
			{
				DualNumber a11 = matrix[0, 0];
				return new DualMatrix(new DualNumber[,] { { 1.0 / a11 } });
			}
			else if (m == 2)
			{
				DualNumber a11 = matrix[0, 0];
				DualNumber a12 = matrix[0, 1];
				DualNumber a21 = matrix[1, 0];
				DualNumber a22 = matrix[1, 1];
				return new DualMatrix(new DualNumber[,] { { a22 / determinant, -a12 / determinant }, { -a21 / determinant, a11 / determinant } });
			}
			else if (m == 3)
			{
				DualNumber a11 = matrix[0, 0];
				DualNumber a12 = matrix[0, 1];
				DualNumber a13 = matrix[0, 2];
				DualNumber a21 = matrix[1, 0];
				DualNumber a22 = matrix[1, 1];
				DualNumber a23 = matrix[1, 2];
				DualNumber a31 = matrix[2, 0];
				DualNumber a32 = matrix[2, 1];
				DualNumber a33 = matrix[2, 2];

				// Using http://mathworld.wolfram.com/MatrixInverse.html.
				return new DualMatrix(new DualNumber[,] {
					{ (a22 * a33 - a23 * a32) / determinant, (a13 * a32 - a12 * a33) / determinant, (a12 * a23 - a13 * a22) / determinant },
					{ (a23 * a31 - a21 * a33) / determinant, (a11 * a33 - a13 * a31) / determinant, (a13 * a21 - a11 * a23) / determinant },
					{ (a21 * a32 - a22 * a31) / determinant, (a12 * a31 - a11 * a32) / determinant, (a11 * a22 - a12 * a21) / determinant }
				});
			}
			else
			{
				throw new NotImplementedException();
			}
		}

		private static int GradientLength(int m, DualMatrix matrix)
		{
			int n = -1;

			for (int i = 0; i < m; i++)
			{
				for (int j = 0; j < m; j++)
				{
					Vector gradient = matrix[i, j].Gradient;
					if (gradient != null)
					{
						if (n != -1 && gradient.Length != n)
						{
							throw new ArgumentException("Inconsistent number of derivatives.");
						}

						n = gradient.Length;
					}
				}
			}

			return n;
		}
	}
}
