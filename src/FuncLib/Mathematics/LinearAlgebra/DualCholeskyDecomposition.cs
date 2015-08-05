// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

using FuncLib.DualNumbers;

namespace FuncLib.Mathematics.LinearAlgebra
{
	/// <summary>
	/// Computes the inverse of a <see cref="DualMatrix" /> by keeping track of the gradient and Hessian throughout the Cholesky decomposition.
	/// </summary>
	public sealed class DualCholeskyDecomposition
	{
		private DualCholeskyDecomposition()
		{
		}

		public static void InverseDeterminant(DualMatrix matrix, out DualMatrix inverse, out DualNumber determinant)
		{
			int n = matrix.Rows;
			if (matrix.Columns != n)
			{
				throw new ArgumentException("The matrix is not a square matrix.");
			}

			DualNumber[,] a = matrix.ToArray();

			if (!spdmatrixcholesky(ref a, n, false))
			{
				throw new ArithmeticException();
			}

			determinant = spdmatrixcholeskydet(ref a, n);

			int info = 0;
			spdmatrixcholeskyinverse(ref a, n, false, ref info);

			for (int i = 0; i < n; i++)
			{
				for (int j = i + 1; j < n; j++)
				{
					a[i, j] = a[j, i];
				}
			}

			inverse = new DualMatrix(a);
		}

		public static DualMatrix Inverse(DualMatrix matrix)
		{
			int n = matrix.Rows;
			if (matrix.Columns != n)
			{
				throw new ArgumentException("The matrix is not a square matrix.");
			}

			DualNumber[,] a = matrix.ToArray();
			int info = 0;
			spdmatrixinverse(ref a, n, false, ref info);
			if (info != 1)
			{
				throw new ArithmeticException();
			}

			for (int i = 0; i < n; i++)
			{
				for (int j = i + 1; j < n; j++)
				{
					a[i, j] = a[j, i];
				}
			}

			return new DualMatrix(a);
		}

		/*************************************************************************
		Copyright (c) 2005-2007, Sergey Bochkanov (ALGLIB project).

		>>> SOURCE LICENSE >>>
		This program is free software; you can redistribute it and/or modify
		it under the terms of the GNU General Public License as published by
		the Free Software Foundation (www.fsf.org); either version 2 of the 
		License, or (at your option) any later version.

		This program is distributed in the hope that it will be useful,
		but WITHOUT ANY WARRANTY; without even the implied warranty of
		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
		GNU General Public License for more details.

		A copy of the GNU General Public License is available at
		http://www.fsf.org/licensing/licenses

		>>> END OF LICENSE >>>
		*************************************************************************/

		private static DualNumber spdmatrixdet(DualNumber[,] a, int n, bool isupper)
		{
			DualNumber result = 0;

			a = (DualNumber[,])a.Clone();

			if (!spdmatrixcholesky(ref a, n, isupper))
			{
				result = -1;
			}
			else
			{
				result = spdmatrixcholeskydet(ref a, n);
			}
			return result;
		}

		private static void spdmatrixinverse(ref DualNumber[,] a, int n, bool isupper, ref int info)
		{
			if (n < 1)
			{
				info = -1;
				return;
			}
			info = 1;
			if (spdmatrixcholesky(ref a, n, isupper))
			{
				spdmatrixcholeskyinverse(ref a, n, isupper, ref info);
			}
			else
			{
				info = -3;
			}
		}

		private static DualNumber spdmatrixcholeskydet(ref DualNumber[,] a, int n)
		{
			DualNumber result = 0;
			int i = 0;

			result = 1;
			for (i = 0; i <= n - 1; i++)
			{
				result = result * DualNumber.Sqr(a[i, i]);
			}
			return result;
		}

		private static bool spdmatrixcholesky(ref DualNumber[,] a, int n, bool isupper)
		{
			bool result = new bool();
			DualNumber[] tmp = new DualNumber[0];

			if (n < 1)
			{
				result = false;
				return result;
			}
			tmp = new DualNumber[2 * n];
			result = spdmatrixcholeskyrec(ref a, 0, n, isupper, ref tmp);
			return result;
		}

		private static bool spdmatrixcholeskyrec(ref DualNumber[,] a, int offs, int n, bool isupper, ref DualNumber[] tmp)
		{
			bool result = new bool();
			//int n1 = 0;
			//int n2 = 0;


			//
			// check N
			//
			if (n < 1)
			{
				result = false;
				return result;
			}

			//
			// special cases
			//
			if (n == 1)
			{
				if (a[offs, offs].Value > 0.0)
				{
					a[offs, offs] = DualNumber.Sqrt(a[offs, offs]);
					result = true;
				}
				else
				{
					result = false;
				}
				return result;
			}
			//if (n <= ablas.ablasblocksize(ref a))
			{
				result = spdmatrixcholesky2(ref a, offs, n, isupper, ref tmp);
				return result;
			}
		}

		private static bool spdmatrixcholesky2(ref DualNumber[,] aaa, int offs, int n, bool isupper, ref DualNumber[] tmp)
		{
			bool result = new bool();
			int i = 0;
			int j = 0;
			//int k = 0;
			//int j1 = 0;
			//int j2 = 0;
			DualNumber ajj = 0;
			DualNumber v = 0;
			//double r = 0;
			int i_ = 0;
			int i1_ = 0;

			result = true;
			if (n < 0)
			{
				result = false;
				return result;
			}

			//
			// Quick return if possible
			//
			if (n == 0)
			{
				return result;
			}
			if (isupper)
			{
				throw new NotImplementedException();
			}
			else
			{

				//
				// Compute the Cholesky factorization A = L*L'.
				//
				for (j = 0; j <= n - 1; j++)
				{

					//
					// Compute L(J+1,J+1) and test for non-positive-definiteness.
					//
					v = 0.0;
					for (i_ = offs; i_ <= offs + j - 1; i_++)
					{
						v += aaa[offs + j, i_] * aaa[offs + j, i_];
					}
					ajj = aaa[offs + j, offs + j] - v;
					if (ajj.Value <= 0.0)
					{
						aaa[offs + j, offs + j] = ajj;
						result = false;
						return result;
					}
					ajj = DualNumber.Sqrt(ajj);
					aaa[offs + j, offs + j] = ajj;

					//
					// Compute elements J+1:N of column J.
					//
					if (j < n - 1)
					{
						if (j > 0)
						{
							i1_ = (offs) - (0);
							for (i_ = 0; i_ <= j - 1; i_++)
							{
								tmp[i_] = aaa[offs + j, i_ + i1_];
							}
							rmatrixmv(n - j - 1, j, ref aaa, offs + j + 1, offs, 0, ref tmp, 0, ref tmp, n);
							for (i = 0; i <= n - j - 2; i++)
							{
								aaa[offs + j + 1 + i, offs + j] = (aaa[offs + j + 1 + i, offs + j] - tmp[n + i]) / ajj;
							}
						}
						else
						{
							for (i = 0; i <= n - j - 2; i++)
							{
								aaa[offs + j + 1 + i, offs + j] = aaa[offs + j + 1 + i, offs + j] / ajj;
							}
						}
					}
				}
			}
			return result;
		}

		private static void rmatrixmv(int m, int n, ref DualNumber[,] a, int ia, int ja, int opa, ref DualNumber[] x, int ix, ref DualNumber[] y, int iy)
		{
			int i = 0;
			DualNumber v = 0;
			int i_ = 0;
			int i1_ = 0;

			if (m == 0)
			{
				return;
			}
			if (n == 0)
			{
				for (i = 0; i <= m - 1; i++)
				{
					y[iy + i] = 0;
				}
				return;
			}
			//if (ablasf.rmatrixmvf(m, n, ref a, ia, ja, opa, ref x, ix, ref y, iy))
			//{
			//    return;
			//}
			if (opa == 0)
			{

				//
				// y = A*x
				//
				for (i = 0; i <= m - 1; i++)
				{
					i1_ = (ix) - (ja);
					v = 0.0;
					for (i_ = ja; i_ <= ja + n - 1; i_++)
					{
						v += a[ia + i, i_] * x[i_ + i1_];
					}
					y[iy + i] = v;
				}
				return;
			}
			if (opa == 1)
			{
				throw new NotImplementedException();
			}
		}

		private static void spdmatrixcholeskyinverse(ref DualNumber[,] a, int n, bool isupper, ref int info)
		{
			//int i = 0;
			//int j = 0;
			//int k = 0;
			//double v = 0;
			//double ajj = 0;
			//double aii = 0;
			DualNumber[] tmp = new DualNumber[0];
			//int info2 = 0;
			//matinvreport rep2 = new matinvreport();

			if (n < 1)
			{
				info = -1;
				return;
			}
			info = 1;

			//
			// calculate condition numbers
			//
			//rep.r1 = rcond.spdmatrixcholeskyrcond(ref a, n, isupper);
			//rep.rinf = rep.r1;
			//if ((double)(rep.r1) < (double)(rcond.rcondthreshold()) | (double)(rep.rinf) < (double)(rcond.rcondthreshold()))
			//{
			//    if (isupper)
			//    {
			//        throw new NotImplementedException();
			//    }
			//    else
			//    {
			//        for (i = 0; i <= n - 1; i++)
			//        {
			//            for (j = 0; j <= i; j++)
			//            {
			//                a[i, j] = 0;
			//            }
			//        }
			//    }
			//    rep.r1 = 0;
			//    rep.rinf = 0;
			//    info = -3;
			//    return;
			//}

			//
			// Inverse
			//
			tmp = new DualNumber[n];
			spdmatrixcholeskyinverserec(ref a, 0, n, isupper, ref tmp);
		}

		private static void spdmatrixcholeskyinverserec(ref DualNumber[,] a, int offs, int n, bool isupper, ref DualNumber[] tmp)
		{
			int i = 0;
			int j = 0;
			DualNumber v = 0;
			//int n1 = 0;
			//int n2 = 0;
			int info2 = 0;
			//matinvreport rep2 = new matinvreport();
			int i_ = 0;
			int i1_ = 0;

			if (n < 1)
			{
				return;
			}

			//
			// Base case
			//
			//if (n <= ablas.ablasblocksize(ref a))
			{
				rmatrixtrinverserec(ref a, offs, n, isupper, false, ref tmp, ref info2);
				if (isupper)
				{
					throw new NotImplementedException();
				}
				else
				{

					//
					// Compute the product L' * L
					// NOTE: we never assume that diagonal of L is real
					//
					for (i = 0; i <= n - 1; i++)
					{
						if (i == 0)
						{

							//
							// 1x1 matrix
							//
							a[offs + i, offs + i] = DualNumber.Sqr(a[offs + i, offs + i]);
						}
						else
						{

							//
							// (I+1)x(I+1) matrix,
							//
							// ( A11^H  A21^H )   ( A11      )   ( A11^H*A11+A21^H*A21  A21^H*A22 )
							// (              ) * (          ) = (                                )
							// (        A22^H )   ( A21  A22 )   ( A22^H*A21            A22^H*A22 )
							//
							// A11 is IxI, A22 is 1x1.
							//
							i1_ = (offs) - (0);
							for (i_ = 0; i_ <= i - 1; i_++)
							{
								tmp[i_] = a[offs + i, i_ + i1_];
							}
							for (j = 0; j <= i - 1; j++)
							{
								v = a[offs + i, offs + j];
								i1_ = (0) - (offs);
								for (i_ = offs; i_ <= offs + j; i_++)
								{
									a[offs + j, i_] = a[offs + j, i_] + v * tmp[i_ + i1_];
								}
							}
							v = a[offs + i, offs + i];
							for (i_ = offs; i_ <= offs + i - 1; i_++)
							{
								a[offs + i, i_] = v * a[offs + i, i_];
							}
							a[offs + i, offs + i] = DualNumber.Sqr(a[offs + i, offs + i]);
						}
					}
				}
				return;
			}
		}

		private static void rmatrixtrinverserec(ref DualNumber[,] a, int offs, int n, bool isupper, bool isunit, ref DualNumber[] tmp, ref int info)
		{
			//int n1 = 0;
			//int n2 = 0;
			int i = 0;
			int j = 0;
			DualNumber v = 0;
			DualNumber ajj = 0;
			int i_ = 0;
			int i1_ = 0;

			if (n < 1)
			{
				info = -1;
				return;
			}

			//
			// Base case
			//
			//if (n <= ablas.ablasblocksize(ref a))
			{
				if (isupper)
				{
					throw new NotImplementedException();
				}
				else
				{

					//
					// Compute inverse of lower triangular matrix.
					//
					for (j = n - 1; j >= 0; j--)
					{
						if (!isunit)
						{
							if (a[offs + j, offs + j].Value == 0.0)
							{
								info = -3;
								return;
							}
							a[offs + j, offs + j] = 1 / a[offs + j, offs + j];
							ajj = -a[offs + j, offs + j];
						}
						else
						{
							ajj = -1;
						}
						if (j < n - 1)
						{

							//
							// Compute elements j+1:n of j-th column.
							//
							i1_ = (offs + j + 1) - (j + 1);
							for (i_ = j + 1; i_ <= n - 1; i_++)
							{
								tmp[i_] = a[i_ + i1_, offs + j];
							}
							for (i = j + 1; i <= n - 1; i++)
							{
								if (i > j + 1)
								{
									i1_ = (j + 1) - (offs + j + 1);
									v = 0.0;
									for (i_ = offs + j + 1; i_ <= offs + i - 1; i_++)
									{
										v += a[offs + i, i_] * tmp[i_ + i1_];
									}
								}
								else
								{
									v = 0;
								}
								if (!isunit)
								{
									a[offs + i, offs + j] = v + a[offs + i, offs + i] * tmp[i];
								}
								else
								{
									a[offs + i, offs + j] = v + tmp[i];
								}
							}
							for (i_ = offs + j + 1; i_ <= offs + n - 1; i_++)
							{
								a[i_, offs + j] = ajj * a[i_, offs + j];
							}
						}
					}
				}
				return;
			}
		}
	}
}
