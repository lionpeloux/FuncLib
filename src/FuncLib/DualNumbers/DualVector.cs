// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Diagnostics;

using FuncLib.Mathematics;

namespace FuncLib.DualNumbers
{
	[Serializable]
	[DebuggerStepThrough]
	[DebuggerDisplay("{ToString(),nq}")]
	public sealed class DualVector
	{
		private DualMatrix inner;

		private DualVector(DualMatrix inner)
		{
			if (inner.Columns != 1)
			{
				throw new ArgumentException();
			}

			this.inner = inner;
		}

		public DualVector(DualNumber[] entries)
		{
			int n = entries.Length;

			DualNumber[,] a = new DualNumber[n, 1];
			for (int i = 0; i < n; i++)
			{
				a[i, 0] = entries[i];
			}

			inner = new DualMatrix(a);
		}

		public DualVector(Vector values)
			: this(values, null)
		{
		}

		public DualVector(Vector values, Vector[] gradients)
			: this(values, gradients, null)
		{
		}

		public DualVector(Vector values, Vector[] gradients, Vector[,] hessians)
		{
			int n = 0;

			Matrix[] matrixGradients = null;
			if (gradients != null)
			{
				n = gradients.Length;

				matrixGradients = new Matrix[n];
				for (int i = 0; i < n; i++)
				{
					if (gradients[i] == null)
					{
						throw new ArgumentNullException("gradients", "The gradients must be fully specified.");
					}

					matrixGradients[i] = (Matrix)gradients[i];
				}
			}

			Matrix[,] matrixHessians = null;
			if (hessians != null)
			{
				if (hessians.GetLength(0) != n || hessians.GetLength(1) != n)
				{
					throw new ArgumentException("Inconsistent number of derivatives.");
				}

				matrixHessians = new Matrix[n, n];
				for (int i = 0; i < n; i++)
				{
					for (int j = 0; j < n; j++)
					{
						if (hessians[i, j] == null)
						{
							throw new ArgumentNullException("hessians", "The Hessians must be fully specified.");
						}

						matrixHessians[i, j] = (Matrix)hessians[i, j];
					}
				}
			}

			inner = new DualMatrix((Matrix)values, matrixGradients, matrixHessians);
		}

		public DualVector SetEntry(int index, DualNumber entry)
		{
			return new DualVector(inner.SetEntry(index, 0, entry));
		}

		public DualVector SetVector(int index, DualVector vector)
		{
			return new DualVector(inner.SetMatrix(index, 0, vector.inner));
		}

		public DualVector GetVector(int index, int length)
		{
			return new DualVector(inner.GetMatrix(index, 0, length, 1));
		}

		public Vector GetValues()
		{
			return new Vector(inner.GetValues());
		}

		public Vector GetGradients(int index)
		{
			return new Vector(inner.GetGradients(index));
		}

		public Vector GetHessians(int index1, int index2)
		{
			return new Vector(inner.GetHessians(index1, index2));
		}

		public DualNumber[] ToArray()
		{
			int n = Length;

			DualNumber[] a = new DualNumber[n];
			for (int i = 0; i < n; i++)
			{
				a[i] = inner[i, 0];
			}

			return a;
		}

		public override string ToString()
		{
			return inner.ToString();
		}

		public string ToString(string format)
		{
			return inner.ToString(format);
		}

		public static DualVector Zero(int length)
		{
			return new DualVector(DualMatrix.Zero(length, 1));
		}

		public static DualVector Basis(int length, int index)
		{
			return new DualVector(DualMatrix.Basis(length, 1, index, 0));
		}

		public static implicit operator DualVector(Vector a)
		{
			return new DualVector(a);
		}

		public static implicit operator DualMatrix(DualVector a)
		{
			return a.inner;
		}

		public static explicit operator DualVector(DualMatrix a)
		{
			if (a.Columns == 1)
			{
				return new DualVector(a);
			}

			if (a.Rows == 1)
			{
				// Requires copying of the references to DualNumber entries.
				return new DualVector(DualMatrix.Transpose(a));
			}

			throw new InvalidCastException("The matrix has no vector representation.");
		}

		public static DualVector operator +(DualVector a, DualVector b)
		{
			return new DualVector(a.inner + b.inner);
		}

		public static DualVector operator +(DualVector a, DualNumber b)
		{
			return new DualVector(a.inner + b);
		}

		public static DualVector operator +(DualNumber a, DualVector b)
		{
			return b + a;
		}

		public static DualVector operator -(DualVector a, DualVector b)
		{
			return new DualVector(a.inner - b.inner);
		}

		public static DualVector operator -(DualVector a, DualNumber b)
		{
			return a + (-b);
		}

		public static DualVector operator -(DualNumber a, DualVector b)
		{
			return new DualVector(a - b.inner);
		}

		public static DualVector operator -(DualVector a)
		{
			return a * -1.0;
		}

		public static DualVector operator *(DualMatrix a, DualVector b)
		{
			return new DualVector(a * b.inner);
		}

		public static DualVector operator *(DualVector a, DualNumber b)
		{
			return new DualVector(a.inner * b);
		}

		public static DualVector operator *(DualNumber a, DualVector b)
		{
			return b * a;
		}

		public static DualVector operator /(DualVector vector, DualNumber scalar)
		{
			return vector * (1.0 / scalar);
		}

		public DualNumber this[int index]
		{
			get
			{
				return inner[index, 0];
			}
		}

		public int Length
		{
			get
			{
				return inner.Rows;
			}
		}
	}
}
