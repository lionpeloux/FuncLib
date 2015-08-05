// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Diagnostics;
using System.Globalization;
using System.Text;

using FuncLib.Mathematics.LinearAlgebra;

namespace FuncLib.Mathematics
{
	[Serializable]
	[DebuggerStepThrough]
	public sealed class Matrix
	{
		private int rows, columns;
		private double[] values;

		private Matrix(int rows, int columns)
		{
			this.rows = rows;
			this.columns = columns;

			values = new double[rows * columns];
		}

		private Matrix(Matrix matrix)
		{
			rows = matrix.rows;
			columns = matrix.columns;
			values = (double[])matrix.values.Clone();
		}

		public Matrix(double[] values, int rows)
		{
			if (rows < 0)
			{
				throw new ArgumentOutOfRangeException();
			}

			if (values.Length > 0)
			{
				this.rows = rows;
				this.columns = values.Length / rows;
				if (rows * columns != values.Length)
				{
					throw new ArgumentException("Array length is not a multiple of the number of rows.");
				}
			}
			this.values = (double[])values.Clone();
		}

		public Matrix(double[,] values)
			: this(values.GetLength(0), values.GetLength(1))
		{
			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < columns; j++)
				{
					SetValueInternal(i, j, values[i, j]);
				}
			}
		}

		public Matrix(int rows, int columns, double value)
			: this(rows, columns)
		{
			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < columns; j++)
				{
					SetValueInternal(i, j, value);
				}
			}
		}

		private void SetValueInternal(int row, int column, double value)
		{
			values[row + column * rows] = value;
		}

		internal double GetValueInternal(int row, int column)
		{
			return values[row + column * rows];
		}

		public Matrix SetValue(int row, int column, double value)
		{
			if (row < 0 || row >= rows || column < 0 || column >= columns)
			{
				throw new IndexOutOfRangeException();
			}

			Matrix resultMatrix = new Matrix(this);
			resultMatrix.SetValueInternal(row, column, value);
			return resultMatrix;
		}

		public Matrix SetMatrix(int row, int column, Matrix matrix)
		{
			if (row < 0 || row + matrix.rows > rows || column < 0 || column + matrix.columns > columns)
			{
				throw new IndexOutOfRangeException();
			}

			Matrix resultMatrix = new Matrix(this);
			for (int i = 0; i < matrix.rows; i++)
			{
				for (int j = 0; j < matrix.columns; j++)
				{
					resultMatrix.SetValueInternal(row + i, column + j, matrix.GetValueInternal(i, j));
				}
			}

			return resultMatrix;
		}

		public Matrix GetMatrix(int row, int column, int rows, int columns)
		{
			if (row < 0 || row + rows > this.rows || column < 0 || column + columns > this.columns)
			{
				throw new IndexOutOfRangeException();
			}

			Matrix resultMatrix = new Matrix(rows, columns);
			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < columns; j++)
				{
					resultMatrix.SetValueInternal(i, j, GetValueInternal(row + i, column + j));
				}
			}

			return resultMatrix;
		}

		public Vector GetRow(int row)
		{
			// return new Vector(GetMatrix(row, 0, rows, 1))

			if (row < 0 || row >= rows)
			{
				throw new IndexOutOfRangeException();
			}

			Matrix resultMatrix = new Matrix(columns, 1);
			for (int i = 0; i < columns; i++)
			{
				resultMatrix.SetValueInternal(i, 0, GetValueInternal(row, i));
			}

			return new Vector(resultMatrix);
		}

		public Vector GetColumn(int column)
		{
			// return new Vector(GetMatrix(0, column, 1, columns))

			if (column < 0 || column >= columns)
			{
				throw new IndexOutOfRangeException();
			}

			Matrix resultMatrix = new Matrix(rows, 1);
			for (int i = 0; i < rows; i++)
			{
				resultMatrix.SetValueInternal(i, 0, GetValueInternal(i, column));
			}

			return new Vector(resultMatrix);
		}

		public double[,] ToArray()
		{
			double[,] values = new double[rows, columns];
			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < columns; j++)
				{
					values[i, j] = GetValueInternal(i, j);
				}
			}

			return values;
		}

		public double[] ToLinearArray()
		{
			return (double[])values.Clone();
		}

		public double[][] ToJaggedArray()
		{
			double[][] values = new double[rows][];
			for (int i = 0; i < rows; i++)
			{
				values[i] = new double[columns];
				for (int j = 0; j < columns; j++)
				{
					values[i][j] = GetValueInternal(i, j);
				}
			}
			return values;
		}

		public override string ToString()
		{
			return ToString(null);
		}

		public string ToString(string format)
		{
			StringBuilder sb = new StringBuilder();
			sb.Append("{");
			for (int i = 0; i < rows; i++)
			{
				if (i > 0)
					sb.Append(", ");
				sb.Append("{");
				for (int j = 0; j < columns; j++)
				{
					if (j > 0)
						sb.Append(", ");
					sb.Append(GetValueInternal(i, j).ToString(format, NumberFormatInfo.InvariantInfo));
				}
				sb.Append("}");
			}
			sb.Append("}");
			return sb.ToString();
		}

		public override int GetHashCode()
		{
			// http://stackoverflow.com/questions/892618/create-a-hashcode-of-two-numbers

			int hash = 23;
			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < columns; j++)
				{
					hash = hash * 31 + GetValueInternal(i, j).GetHashCode();
				}
			}

			return hash;
		}

		public static Matrix Zero(int rows, int columns)
		{
			Matrix resultMatrix = new Matrix(rows, columns);
			resultMatrix.FillZero();
			return resultMatrix;
		}

		public static Matrix Identity(int rows, int columns)
		{
			Matrix resultMatrix = new Matrix(rows, columns);
			resultMatrix.FillZero();
			for (int i = 0; i < rows && i < columns; i++)
			{
				resultMatrix.SetValueInternal(i, i, 1.0);
			}

			return resultMatrix;
		}

		public static Matrix Diagonal(double[] values)
		{
			Matrix resultMatrix = new Matrix(values.Length, values.Length);
			resultMatrix.FillZero();
			for (int i = 0; i < values.Length; i++)
			{
				resultMatrix.SetValueInternal(i, i, values[i]);
			}

			return resultMatrix;
		}

		public static Matrix Basis(int rows, int columns, int row, int column)
		{
			if (row < 0 || row >= rows || column < 0 || column >= columns)
			{
				throw new IndexOutOfRangeException();
			}

			Matrix resultMatrix = new Matrix(rows, columns);
			resultMatrix.FillZero();
			resultMatrix.SetValueInternal(row, column, 1.0);
			return resultMatrix;
		}

		public static Matrix Transpose(Matrix matrix)
		{
			Matrix resultMatrix = new Matrix(matrix.columns, matrix.rows);
			for (int i = 0; i < matrix.rows; i++)
			{
				for (int j = 0; j < matrix.columns; j++)
				{
					resultMatrix.SetValueInternal(j, i, matrix.GetValueInternal(i, j));
				}
			}

			return resultMatrix;
		}

		/// <summary>
		/// General-purpose matrix inversion through LU decomposition. No exceptions are thrown if the matrix isn't invertible.
		/// </summary>
		public static Matrix Inverse(Matrix a)
		{
			int n = a.Rows;
			if (a.Columns != n)
			{
				throw new ArgumentException("The matrix is not a square matrix.");
			}

			Matrix b;
			if (!LUDecomposition.TryInverse(a, out b))
			{
				b = new Matrix(n, n, double.NaN);
			}

			return b;
		}

		/// <summary>
		/// General-purpose matrix determinant through LU decomposition.
		/// </summary>
		public static double Determinant(Matrix a)
		{
			return LUDecomposition.Determinant(a);
		}

		public static Matrix operator +(Matrix leftMatrix, Matrix rightMatrix)
		{
			Matrix resultMatrix = new Matrix(leftMatrix.rows, leftMatrix.columns);

			if (resultMatrix.rows != rightMatrix.rows || resultMatrix.columns != rightMatrix.columns)
			{
				throw new ArgumentException("Matrix sum undefined. Size mismatch.");
			}

			for (int i = 0; i < resultMatrix.rows; i++)
			{
				for (int j = 0; j < resultMatrix.columns; j++)
				{
					resultMatrix.SetValueInternal(i, j, leftMatrix.GetValueInternal(i, j) + rightMatrix.GetValueInternal(i, j));
				}
			}

			return resultMatrix;
		}

		public static Matrix operator +(Matrix matrix, double scalar)
		{
			Matrix resultMatrix = new Matrix(matrix.rows, matrix.columns);
			for (int i = 0; i < resultMatrix.rows; i++)
			{
				for (int j = 0; j < resultMatrix.columns; j++)
				{
					resultMatrix.SetValueInternal(i, j, matrix.GetValueInternal(i, j) + scalar);
				}
			}

			return resultMatrix;
		}

		public static Matrix operator +(double scalar, Matrix matrix)
		{
			return matrix + scalar;
		}

		public static Matrix operator -(Matrix leftMatrix, Matrix rightMatrix)
		{
			Matrix resultMatrix = new Matrix(leftMatrix.rows, leftMatrix.columns);

			if (resultMatrix.rows != rightMatrix.rows || resultMatrix.columns != rightMatrix.columns)
			{
				throw new ArgumentException("Matrix difference undefined. Size mismatch.");
			}

			for (int i = 0; i < resultMatrix.rows; i++)
			{
				for (int j = 0; j < resultMatrix.columns; j++)
				{
					resultMatrix.SetValueInternal(i, j, leftMatrix.GetValueInternal(i, j) - rightMatrix.GetValueInternal(i, j));
				}
			}

			return resultMatrix;
		}

		public static Matrix operator -(Matrix matrix, double scalar)
		{
			return matrix + (-scalar);
		}

		public static Matrix operator -(double scalar, Matrix matrix)
		{
			Matrix resultMatrix = new Matrix(matrix.rows, matrix.columns);
			for (int i = 0; i < resultMatrix.rows; i++)
			{
				for (int j = 0; j < resultMatrix.columns; j++)
				{
					resultMatrix.SetValueInternal(i, j, scalar - matrix.GetValueInternal(i, j));
				}
			}

			return resultMatrix;
		}

		public static Matrix operator -(Matrix matrix)
		{
			return matrix * -1.0;
		}

		public static Matrix operator *(Matrix leftMatrix, Matrix rightMatrix)
		{
			int n = leftMatrix.columns;
			if (n != rightMatrix.rows)
			{
				throw new ArgumentException("Matrix product undefined. Size mismatch.");
			}

			Matrix resultMatrix = new Matrix(leftMatrix.rows, rightMatrix.columns);
			for (int i = 0; i < resultMatrix.rows; i++)
			{
				for (int j = 0; j < resultMatrix.columns; j++)
				{
					double s = 0.0;
					for (int k = 0; k < n; k++)
					{
						s += leftMatrix.GetValueInternal(i, k) * rightMatrix.GetValueInternal(k, j);
					}

					resultMatrix.SetValueInternal(i, j, s);
				}
			}

			return resultMatrix;
		}

		public static Matrix operator *(Matrix matrix, double scalar)
		{
			Matrix resultMatrix = new Matrix(matrix.rows, matrix.columns);
			for (int i = 0; i < resultMatrix.rows; i++)
			{
				for (int j = 0; j < resultMatrix.columns; j++)
				{
					resultMatrix.SetValueInternal(i, j, matrix.GetValueInternal(i, j) * scalar);
				}
			}

			return resultMatrix;
		}

		public static Matrix operator *(double scalar, Matrix matrix)
		{
			return matrix * scalar;
		}

		public static Matrix operator /(Matrix matrix, double scalar)
		{
			return matrix * (1.0 / scalar);
		}

		public double this[int row, int column]
		{
			get
			{
				if (row < 0 || row >= rows || column < 0 || column >= columns)
				{
					throw new IndexOutOfRangeException();
				}

				return GetValueInternal(row, column);
			}
		}

		public int Rows
		{
			get
			{
				return rows;
			}
		}

		public int Columns
		{
			get
			{
				return columns;
			}
		}

		private void FillZero()
		{
			// Not needed. Arrays of value-types are guaranteed to be filled with 0.
		}
	}
}
