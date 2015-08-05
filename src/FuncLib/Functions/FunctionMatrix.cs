// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;

using FuncLib.Functions.Compilation;
using FuncLib.Mathematics;

namespace FuncLib.Functions
{
	/// <summary>
	/// A matrix class with <see cref="Function" /> entries. Performance is not as critical as for the <see cref="Matrix" /> class due to the significant
	/// overhead of object representation of the entries, so this class is just a wrapper around a Function[,] array.
	/// </summary>
	[Serializable]
	public class FunctionMatrix
	{
		private Function[,] entries;
		//private int rows, columns;

		public FunctionMatrix(Function[,] entries)
		{
			// Copy to keep immutable.
			this.entries = (Function[,])entries.Clone();

			//rows = entries.GetLength(0);
			//columns = entries.GetLength(1);
		}

		/// <summary>
		/// New row matrix.
		/// </summary>
		public FunctionMatrix(Function[] entries)
			: this(entries.Length, 1)
		{
			for (int i = 0; i < entries.Length; i++)
			{
				this.entries[i, 0] = entries[i];
			}
		}

		private FunctionMatrix(int rows, int columns)
		{
			entries = new Function[rows, columns];
		}

		private FunctionMatrix(FunctionMatrix a)
			: this(a.entries)
		{
		}

		public Matrix Value(IPoint point)
		{
			double[,] values = new double[Rows, Columns];
			for (int i = 0; i < Rows; i++)
			{
				for (int j = 0; j < Columns; j++)
				{
					values[i, j] = entries[i, j].Value(point);
				}
			}

			return new Matrix(values);
		}

		public Matrix Value(params VariableAssignment[] assignments)
		{
			return Value(new Point(assignments));
		}

		public FunctionMatrix Derivative(Variable variable)
		{
			Function[,] derivatives = new Function[Rows, Columns];
			for (int i = 0; i < Rows; i++)
			{
				for (int j = 0; j < Columns; j++)
				{
					derivatives[i, j] = entries[i, j].Derivative(variable);
				}
			}

			return new FunctionMatrix(derivatives);
		}

		/*public MatrixFunction Derivative(Variable variable, int order)
		{
			if (order < 0)
			{
				throw new ArgumentOutOfRangeException("order");
			}

			MatrixFunction a = this;
			for (int i = 0; i < order; i++)
			{
				a = a.Derivative(variable);
			}
			return a;
		}

		public MatrixFunction Derivative(params Variable[] variables)
		{
			MatrixFunction a = this;
			for (int i = 0; i < variables.Length; i++)
			{
				a = a.Derivative(variables[i]);
			}
			return a;
		}*/

		public Function[,] ToArray()
		{
			return (Function[,])entries.Clone();
		}

		public FunctionMatrix SetValue(int row, int column, Function f)
		{
			if (row < 0 || row >= Rows || column < 0 || column >= Columns)
			{
				throw new IndexOutOfRangeException();
			}

			FunctionMatrix a = new FunctionMatrix(this);
			a.entries[row, column] = f;
			return a;
		}

		public static FunctionMatrix Zero(int rows, int columns)
		{
			FunctionMatrix a = new FunctionMatrix(rows, columns);
			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < columns; j++)
				{
					a.entries[i, j] = 0.0;
				}
			}

			return a;
		}

		public static FunctionMatrix Identity(int rows, int columns)
		{
			FunctionMatrix a = FunctionMatrix.Zero(rows, columns); //new FunctionMatrix(rows, columns);
			for (int i = 0; i < rows && i < columns; i++)
			{
				a.entries[i, i] = 1.0;
			}

			return a;
		}

		public static FunctionMatrix Diagonal(Function[] values)
		{
			int n = values.Length;

			FunctionMatrix a = FunctionMatrix.Zero(n, n); //new FunctionMatrix(values.Length, values.Length);
			for (int i = 0; i < n; i++)
			{
				a.entries[i, i] = values[i];
			}

			return a;
		}

		public static FunctionMatrix Basis(int rows, int columns, int row, int column)
		{
			if (row < 0 || row >= rows || column < 0 || column >= columns)
			{
				throw new IndexOutOfRangeException();
			}

			FunctionMatrix a = FunctionMatrix.Zero(rows, columns); //new FunctionMatrix(rows, columns);
			a.entries[row, column] = 1.0;
			return a;
		}

		public static FunctionMatrix Transpose(FunctionMatrix a)
		{
			FunctionMatrix b = new FunctionMatrix(a.Columns, a.Rows);
			for (int i = 0; i < a.Rows; i++)
			{
				for (int j = 0; j < a.Columns; j++)
				{
					b.entries[j, i] = a.entries[i, j];
				}
			}

			return b;
		}

		public static Function Trace(FunctionMatrix a)
		{
			Function s = 0.0;
			for (int i = 0; i < a.Rows && i < a.Columns; i++)
			{
				s += a[i, i];
			}

			return s;
		}

		public static FunctionMatrix Inverse(FunctionMatrix a)
		{
			FunctionMatrix b;
			Function d;
			InverseDeterminant(a, out b, out d);
			return b;
		}

		public static Function Determinant(FunctionMatrix a)
		{
			FunctionMatrix b;
			Function d;
			InverseDeterminant(a, out b, out d);
			return d;
		}

		/// <summary>
		/// Compute the inverse and the determinant simultaneously. This avoids computing similar expressions twice.
		/// </summary>
		/// <param name="a">The matrix to compute inverse and determinant of.</param>
		/// <param name="b">The inverse.</param>
		/// <param name="d">The determinant.</param>
		public static void InverseDeterminant(FunctionMatrix a, out FunctionMatrix b, out Function d)
		{
			// Avoid using complex and slow MatrixFunction objects small matrices with formulas for inverses.

			if (a.Rows == 1 && a.Columns == 1)
			{
				Function a11 = a[0, 0];

				d = a11;
				b = new FunctionMatrix(new Function[,] { { 1.0 / a11 } });

				return;
			}

			if (a.Rows == 2 && a.Columns == 2) // FIXME Enable only at final optimization
			{
				Function a11 = a[0, 0];
				Function a12 = a[0, 1];
				Function a21 = a[1, 0];
				Function a22 = a[1, 1];

				d = a11 * a22 - a12 * a21;
				b = new FunctionMatrix(new Function[,] { { a22 / d, -a12 / d }, { -a21 / d, a11 / d } });

				return;
			}

			if (a.Rows == 3 && a.Columns == 3)
			{
				Function a11 = a[0, 0];
				Function a12 = a[0, 1];
				Function a13 = a[0, 2];
				Function a21 = a[1, 0];
				Function a22 = a[1, 1];
				Function a23 = a[1, 2];
				Function a31 = a[2, 0];
				Function a32 = a[2, 1];
				Function a33 = a[2, 2];

				// Using http://mathworld.wolfram.com/MatrixInverse.html.
				d = a11 * a22 * a33 + a12 * a23 * a31 + a13 * a21 * a32 - a11 * a23 * a32 - a12 * a21 * a33 - a13 * a22 * a31;
				b = new FunctionMatrix(new Function[,] {
					{ (a22 * a33 - a23 * a32) / d, (a13 * a32 - a12 * a33) / d, (a12 * a23 - a13 * a22) / d },
					{ (a23 * a31 - a21 * a33) / d, (a11 * a33 - a13 * a31) / d, (a13 * a21 - a11 * a23) / d },
					{ (a21 * a32 - a22 * a31) / d, (a12 * a31 - a11 * a32) / d, (a11 * a22 - a12 * a21) / d }
				});

				return;
			}

			b = new InverseMatrixFunction(a);
			d = new DeterminantFunction(a, b);
		}

		public static FunctionMatrix operator +(FunctionMatrix a, FunctionMatrix b)
		{
			int n = a.Rows;
			int m = a.Columns;
			if (b.Rows != n || b.Columns != m)
			{
				throw new ArgumentException("Matrix sum undefined. Size mismatch.");
			}

			FunctionMatrix c = new FunctionMatrix(n, m);
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					c.entries[i, j] = a.entries[i, j] + b.entries[i, j];
				}
			}

			return c;
		}

		public static FunctionMatrix operator +(FunctionMatrix a, Function f)
		{
			int n = a.Rows;
			int m = a.Columns;

			FunctionMatrix b = new FunctionMatrix(n, m);
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					b.entries[i, j] = a.entries[i, j] + f;
				}
			}

			return b;
		}

		public static FunctionMatrix operator +(Function f, FunctionMatrix a)
		{
			return a + f;
		}

		public static FunctionMatrix operator -(FunctionMatrix a, FunctionMatrix b)
		{
			int n = a.Rows;
			int m = a.Columns;
			if (b.Rows != n || b.Columns != m)
			{
				throw new ArgumentException("Matrix sum undefined. Size mismatch.");
			}

			FunctionMatrix c = new FunctionMatrix(n, m);
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					c.entries[i, j] = a.entries[i, j] - b.entries[i, j];
				}
			}

			return c;
		}

		public static FunctionMatrix operator -(FunctionMatrix a, Function f)
		{
			int n = a.Rows;
			int m = a.Columns;

			FunctionMatrix b = new FunctionMatrix(n, m);
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					b.entries[i, j] = a.entries[i, j] - f;
				}
			}

			return b;
		}

		public static FunctionMatrix operator -(Function f, FunctionMatrix a)
		{
			int n = a.Rows;
			int m = a.Columns;

			FunctionMatrix b = new FunctionMatrix(n, m);
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					b.entries[i, j] = f - a.entries[i, j];
				}
			}

			return b;
		}

		public static FunctionMatrix operator -(FunctionMatrix a)
		{
			int n = a.Rows;
			int m = a.Columns;

			FunctionMatrix b = new FunctionMatrix(n, m);
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					b.entries[i, j] = -a.entries[i, j];
				}
			}

			return b;
		}

		public static FunctionMatrix operator *(FunctionMatrix a, FunctionMatrix b)
		{
			int n = a.Columns;
			if (n != b.Rows)
			{
				throw new ArgumentException("Matrix product undefined. Size mismatch.");
			}

			Function[,] c = new Function[a.Rows, b.Columns];
			for (int i = 0; i < a.Rows; i++)
			{
				for (int j = 0; j < b.Columns; j++)
				{
					Function s = 0.0;
					for (int k = 0; k < n; k++)
					{
						s += a[i, k] * b[k, j];
					}

					c[i, j] = s;
				}
			}

			return new FunctionMatrix(c);
		}

		public static FunctionMatrix operator *(FunctionMatrix a, Function f)
		{
			int n = a.Rows;
			int m = a.Columns;

			FunctionMatrix b = new FunctionMatrix(n, m);
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					b.entries[i, j] = a.entries[i, j] * f;
				}
			}

			return b;
		}

		public static FunctionMatrix operator *(Function f, FunctionMatrix a)
		{
			return a * f;
		}

		public static FunctionMatrix operator /(FunctionMatrix a, Function f)
		{
			int n = a.Rows;
			int m = a.Columns;

			FunctionMatrix b = new FunctionMatrix(n, m);
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					b.entries[i, j] = a.entries[i, j] / f;
				}
			}

			return b;
		}

		public static implicit operator FunctionMatrix(Matrix a)
		{
			int n = a.Rows;
			int m = a.Columns;

			FunctionMatrix b = new FunctionMatrix(n, m);
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					b.entries[i, j] = a[i, j];
				}
			}

			return b;
		}

		public Function this[int row, int column]
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

		[Serializable]
		private class InverseMatrixFunction : FunctionMatrix
		{
			private FunctionMatrix a;
			private EvaluationCache<Matrix> inverses;
			private Dictionary<Variable, FunctionMatrix> derivatives;

			public InverseMatrixFunction(FunctionMatrix a)
				: base(BuildEntries(a))
			{
				this.a = a;

				// Add a reference to this object to all entries create above. Rather ugly construction but it
				// works perfectly for now.
				for (int i = 0; i < Rows; i++)
				{
					for (int j = 0; j < Columns; j++)
					{
						((ValueFunction)entries[i, j]).UpdateParent(this);
					}
				}

				inverses = new EvaluationCache<Matrix>();
				derivatives = new Dictionary<Variable, FunctionMatrix>();
			}

			private static Function[,] BuildEntries(FunctionMatrix a)
			{
				int n = a.Rows;
				if (a.Columns != n)
				{
					throw new ArgumentException();
				}

				Function[,] entries = new Function[n, n];
				for (int i = 0; i < n; i++)
				{
					for (int j = 0; j < n; j++)
					{
						entries[i, j] = new ValueFunction(i, j);
					}
				}

				return entries;
			}

			private Matrix ComputeValue(IEvaluation evaluation)
			{
				Matrix inverse;
				if (!inverses.TryGetValue(evaluation, out inverse))
				{
					inverse = Matrix.Inverse(a.Value(evaluation));
					inverses.Add(evaluation, inverse);
				}

				return inverse;
			}

			private FunctionMatrix ComputeDerivative(Variable variable)
			{
				// It's important for efficient evaluation to use the same object for all entries.

				FunctionMatrix derivative;
				if (!derivatives.TryGetValue(variable, out derivative))
				{
					// (53) in The Matrix Cookbook.
					derivative = -this * a.Derivative(variable) * this;
					derivatives.Add(variable, derivative);
				}

				return derivative;
			}

			[Serializable]
			private class ValueFunction : Function
			{
				private InverseMatrixFunction parent;
				private int row, column;

				public ValueFunction(int row, int column)
				{
					this.row = row;
					this.column = column;
				}

				public void UpdateParent(InverseMatrixFunction parent)
				{
					this.parent = parent;
				}

				protected override double ComputeValue(IEvaluation evaluation)
				{
					return parent.ComputeValue(evaluation)[row, column];
				}

				protected override Function ComputeDerivative(Variable variable)
				{
					return parent.ComputeDerivative(variable)[row, column];
				}

				public override Expression Compile(CodeGenerator generator)
				{
					// public override MatrixExpression Compile() = Her ville multi-evaluarings compile være brugbar!!
					throw new NotImplementedException();
				}
			}
		}

		[Serializable]
		private class DeterminantFunction : Function
		{
			private FunctionMatrix a, ainv;

			public DeterminantFunction(FunctionMatrix a, FunctionMatrix ainv)
			{
				if (a.Rows != a.Columns)
				{
					throw new ArgumentException();
				}

				this.a = a;
				this.ainv = ainv;
			}

			protected override double ComputeValue(IEvaluation evaluation)
			{
				return Matrix.Determinant(a.Value(evaluation));
			}

			protected override Function ComputeDerivative(Variable variable)
			{
				// (41) in The Matrix Cookbook.
				return this * FunctionMatrix.Trace(ainv * a.Derivative(variable));
			}

			public override Expression Compile(CodeGenerator generator)
			{
				throw new NotImplementedException();
			}
		}
	}
}
