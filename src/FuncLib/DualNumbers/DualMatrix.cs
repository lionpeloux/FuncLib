// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Diagnostics;
using System.Text;

using FuncLib.Mathematics;
using FuncLib.Mathematics.LinearAlgebra;

namespace FuncLib.DualNumbers
{
	[Serializable]
	[DebuggerStepThrough]
	[DebuggerDisplay("{ToString(),nq}")]
	public sealed class DualMatrix
	{
		private DualNumber[,] entries;

		private DualMatrix(int rows, int columns)
		{
			if (rows < 0 || columns < 0)
			{
				throw new ArgumentOutOfRangeException();
			}

			entries = new DualNumber[rows, columns];
		}

		public DualMatrix(DualNumber[,] entries)
		{
			entries = (DualNumber[,])entries.Clone();
		}

		public DualMatrix(Matrix values)
			: this(values, null)
		{
		}

		public DualMatrix(Matrix values, Matrix[] gradients)
			: this(values, gradients, null)
		{
		}

		public DualMatrix(Matrix values, Matrix[] gradients, Matrix[,] hessians)
		{
			int rows = values.Rows;
			int columns = values.Columns;

			entries = new DualNumber[rows, columns];

			int n = 0;

			if (gradients != null)
			{
				n = gradients.Length;

				for (int i = 0; i < n; i++)
				{
					if (gradients[i] == null)
					{
						throw new ArgumentNullException("gradients", "The gradients must be fully specified.");
					}

					if (gradients[i].Rows != rows || gradients[i].Columns != columns)
					{
						throw new ArgumentException("Inconsistent matrix sizes.");
					}
				}
			}

			if (hessians != null)
			{
				if (gradients == null)
				{
					throw new ArgumentException("The gradients must be specified if the Hessians are specified.");
				}

				if (hessians.GetLength(0) != n || hessians.GetLength(1) != n)
				{
					throw new ArgumentException("Inconsistent number of derivatives.");
				}

				for (int i = 0; i < n; i++)
				{
					for (int j = 0; j < n; j++)
					{
						if (hessians[i, j] == null)
						{
							throw new ArgumentNullException("hessians", "The Hessians must be fully specified.");
						}

						if (hessians[i, j].Rows != rows || hessians[i, j].Columns != columns)
						{
							throw new ArgumentException("Inconsistent matrix sizes.");
						}
					}
				}
			}

			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < columns; j++)
				{
					double value = values[i, j];

					Vector gradient = null;
					if (gradients != null)
					{
						double[] a = new double[n];
						for (int k = 0; k < n; k++)
						{
							a[k] = gradients[k][i, j];
						}
						gradient = new Vector(a);
					}

					Matrix hessian = null;
					if (hessians != null)
					{
						double[,] a = new double[n, n];
						for (int k = 0; k < n; k++)
						{
							for (int l = 0; l < n; l++)
							{
								a[k, l] = hessians[k, l][i, j];
							}
						}
						hessian = new Matrix(a);
					}

					entries[i, j] = new DualNumber(value, gradient, hessian);
				}
			}
		}

		public DualMatrix SetEntry(int row, int column, DualNumber entry)
		{
			DualMatrix a = new DualMatrix(entries);
			a.entries[row, column] = entry;

			return a;
		}

		public DualMatrix SetMatrix(int row, int column, DualMatrix matrix)
		{
			int n = matrix.Rows;
			int m = matrix.Columns;

			DualMatrix b = new DualMatrix(entries);
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					b.entries[row + i, column + j] = matrix.entries[i, j];
				}
			}

			return b;
		}

		public DualMatrix GetMatrix(int row, int column, int rows, int columns)
		{
			DualMatrix a = new DualMatrix(rows, columns);
			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < columns; j++)
				{
					a.entries[i, j] = entries[row + i, column + j];
				}
			}

			return a;
		}

		public DualVector GetRow(int row)
		{
			// Use conversion operator defined in DualVector.
			return (DualVector)(GetMatrix(row, 0, Rows, 1));
		}

		public DualVector GetColumn(int column)
		{
			// Use conversion operator defined in DualVector.
			return (DualVector)(GetMatrix(0, column, 1, Columns));
		}

		public Matrix GetValues()
		{
			int n = Rows;
			int m = Columns;

			double[,] a = new double[n, m];
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					a[i, j] = entries[i, j].Value;
				}
			}

			return new Matrix(a);
		}

		public Matrix GetGradients(int index)
		{
			int n = Rows;
			int m = Columns;

			double[,] a = new double[n, m];
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					Vector gradient = entries[i, j].Gradient;
					if (gradient != null)
					{
						a[i, j] = gradient[index];
					}
				}
			}

			return new Matrix(a);
		}

		public Matrix GetHessians(int index1, int index2)
		{
			int n = Rows;
			int m = Columns;

			double[,] a = new double[n, m];
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					Matrix hessian = entries[i, j].Hessian;
					if (hessian != null)
					{
						a[i, j] = hessian[index1, index2];
					}
				}
			}

			return new Matrix(a);
		}

		public DualNumber[,] ToArray()
		{
			return (DualNumber[,])entries.Clone();
		}

		public override string ToString()
		{
			return ToString(null);
		}

		public string ToString(string format)
		{
			StringBuilder sb = new StringBuilder();
			sb.Append("{");
			for (int i = 0; i < Rows; i++)
			{
				if (i > 0)
				{
					sb.Append(", ");
				}
				sb.Append("{");
				for (int j = 0; j < Columns; j++)
				{
					if (j > 0)
					{
						sb.Append(", ");
					}
					sb.Append(entries[i, j].Value.ToString(format));
				}
				sb.Append("}");
			}
			sb.Append("}");
			return sb.ToString();
		}

		public static DualMatrix Zero(int rows, int columns)
		{
			DualMatrix a = new DualMatrix(rows, columns);
			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < columns; j++)
				{
					a.entries[i, j] = 0.0;
				}
			}

			return a;
		}

		public static DualMatrix Identity(int rows, int columns)
		{
			DualMatrix a = new DualMatrix(rows, columns);
			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < columns; j++)
				{
					a.entries[i, j] = i == j ? 1.0 : 0.0;
				}
			}

			return a;
		}

		public static DualMatrix Diagonal(DualNumber[] entries)
		{
			int n = entries.Length;

			DualMatrix a = new DualMatrix(n, n);
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					a.entries[i, j] = i == j ? entries[i] : 0.0;
				}
			}

			return a;
		}

		public static DualMatrix Basis(int rows, int columns, int row, int column)
		{
			if (row < 0 || row >= rows || column < 0 || column >= columns)
			{
				throw new ArgumentOutOfRangeException();
			}

			DualMatrix a = new DualMatrix(rows, columns);
			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < columns; j++)
				{
					a.entries[i, j] = i == row && j == column ? 1.0 : 0.0;
				}
			}

			return a;
		}

		public static DualMatrix Transpose(DualMatrix a)
		{
			int n = a.Rows;
			int m = a.Columns;

			DualMatrix b = new DualMatrix(m, n);
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					b.entries[j, i] = a.entries[i, j];
				}
			}

			return b;
		}

		public static DualMatrix Inverse(DualMatrix a)
		{
			// Use direct computation for low dimenstions; LU decomposition otherwise.
			return DualLUDecomposition.Inverse(a);
		}

		public static DualNumber Determinant(DualMatrix a)
		{
			// Use direct computation for low dimenstions; LU decomposition otherwise.
			return DualLUDecomposition.Determinant(a);
		}

		public static implicit operator DualMatrix(Matrix a)
		{
			return new DualMatrix(a);
		}

		public static DualMatrix operator +(DualMatrix a, DualMatrix b)
		{
			int n = a.Rows;
			int m = a.Rows;

			if (b.Rows != n || b.Columns != m)
			{
				throw new ArgumentException("Matrix sum undefined. Size mismatch.");
			}

			DualMatrix c = new DualMatrix(n, m);
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					c.entries[i, j] = a.entries[i, j] + b.entries[i, j];
				}
			}

			return c;
		}

		public static DualMatrix operator +(DualMatrix a, DualNumber b)
		{
			int n = a.Rows;
			int m = a.Rows;

			DualMatrix c = new DualMatrix(n, m);
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					c.entries[i, j] = a.entries[i, j] + b;
				}
			}

			return c;
		}

		public static DualMatrix operator +(DualNumber a, DualMatrix b)
		{
			return b + a;
		}

		public static DualMatrix operator -(DualMatrix a, DualMatrix b)
		{
			int n = a.Rows;
			int m = a.Rows;

			if (b.Rows != n || b.Columns != m)
			{
				throw new ArgumentException("Matrix difference undefined. Size mismatch.");
			}

			DualMatrix c = new DualMatrix(n, m);
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					c.entries[i, j] = a.entries[i, j] - b.entries[i, j];
				}
			}

			return c;
		}

		public static DualMatrix operator -(DualMatrix a, DualNumber b)
		{
			return a + (-b);
		}

		public static DualMatrix operator -(DualNumber a, DualMatrix b)
		{
			int n = b.Rows;
			int m = b.Columns;

			DualMatrix c = new DualMatrix(n, m);
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					c.entries[i, j] = a - b.entries[i, j];
				}
			}

			return c;
		}

		public static DualMatrix operator -(DualMatrix a)
		{
			// FIXME optimize
			return a * -1.0;
		}

		public static DualMatrix operator *(DualMatrix leftMatrix, DualMatrix rightMatrix)
		{
			int n = leftMatrix.Columns;
			if (n != rightMatrix.Rows)
			{
				throw new ArgumentException("Matrix product undefined. Size mismatch.");
			}

			DualMatrix resultMatrix = new DualMatrix(leftMatrix.Rows, rightMatrix.Columns);
			for (int i = 0; i < resultMatrix.Rows; i++)
			{
				for (int j = 0; j < resultMatrix.Columns; j++)
				{
					DualNumber s = 0.0;
					for (int k = 0; k < n; k++)
					{
						s += leftMatrix.entries[i, k] * rightMatrix.entries[k, j];
					}

					resultMatrix.entries[i, j] = s;
				}
			}

			return resultMatrix;
		}

		public static DualMatrix operator *(DualMatrix a, DualNumber b)
		{
			int n = a.Rows;
			int m = a.Columns;

			DualMatrix c = new DualMatrix(n, m);
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					c.entries[i, j] = a.entries[i, j] * b;
				}
			}

			return c;
		}

		public static DualMatrix operator *(DualNumber a, DualMatrix b)
		{
			return b * a;
		}

		public static DualMatrix operator /(DualMatrix a, DualNumber b)
		{
			return a * (1.0 / b);
		}

		public DualNumber this[int row, int column]
		{
			get
			{
				return entries[row, column];
			}
		}

		public int Rows
		{
			get
			{
				return entries.GetLength(0);
			}
		}

		public int Columns
		{
			get
			{
				return entries.GetLength(1);
			}
		}
	}
}
