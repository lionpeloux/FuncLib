// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Diagnostics;

namespace FuncLib.Mathematics
{
	[Serializable]
	[DebuggerStepThrough]
	public sealed class Vector
	{
		private Matrix innerMatrix;

		internal Vector(Matrix innerMatrix)
		{
			if (innerMatrix.Columns != 1)
			{
				throw new ArgumentException();
			}

			this.innerMatrix = innerMatrix;
		}

		public Vector(double[] values)
		{
			innerMatrix = new Matrix(values, values.Length);
		}

		public Vector(int length, double value)
		{
			innerMatrix = new Matrix(length, 1, value);
		}

		internal double GetValueInternal(int index)
		{
			return innerMatrix.GetValueInternal(index, 0);
		}

		public Vector SetValue(int index, double value)
		{
			return new Vector(innerMatrix.SetValue(index, 0, value));
		}

		public Vector SetVector(int index, Vector vector)
		{
			return new Vector(innerMatrix.SetMatrix(index, 0, vector.innerMatrix));
		}

		public Vector GetVector(int index, int length)
		{
			return new Vector(innerMatrix.GetMatrix(index, 0, length, 1));
		}

		public double[] ToArray()
		{
			return innerMatrix.ToLinearArray();
		}

		public override string ToString()
		{
			return ToString(null);
		}

		public string ToString(string format)
		{
			string s = Matrix.Transpose(innerMatrix).ToString(format);
			return s.Substring(1, s.Length - 2);
		}

		public override int GetHashCode()
		{
			return innerMatrix.GetHashCode();
		}

		public static Vector Zero(int length)
		{
			return new Vector(Matrix.Zero(length, 1));
		}

		public static Vector Basis(int length, int index)
		{
			return new Vector(Matrix.Basis(length, 1, index, 0));
		}

		public static double Dot(Vector leftVector, Vector rightVector)
		{
			int n = leftVector.Length;
			if (n != rightVector.Length)
			{
				throw new ArgumentException("Dot product undefined. Size mismatch.");
			}

			double s = 0.0;
			for (int i = 0; i < n; i++)
			{
				s += leftVector.GetValueInternal(i) * rightVector.GetValueInternal(i);
			}

			return s;
		}

		public static Vector operator +(Vector leftVector, Vector rightVector)
		{
			return new Vector(leftVector.innerMatrix + rightVector.innerMatrix);
		}

		public static Vector operator +(Vector vector, double scalar)
		{
			return new Vector(vector.innerMatrix + scalar);
		}

		public static Vector operator +(double scalar, Vector vector)
		{
			return vector + scalar;
		}

		public static Vector operator -(Vector leftVector, Vector rightVector)
		{
			return new Vector(leftVector.innerMatrix - rightVector.innerMatrix);
		}

		public static Vector operator -(Vector vector, double scalar)
		{
			return vector + (-scalar);
		}

		public static Vector operator -(double scalar, Vector vector)
		{
			return new Vector(scalar - vector.innerMatrix);
		}

		public static Vector operator -(Vector vector)
		{
			return vector * -1.0;
		}

		public static Vector operator *(Matrix leftMatrix, Vector rightVector)
		{
			return new Vector(leftMatrix * rightVector.innerMatrix);
		}

		public static Vector operator *(Vector vector, double scalar)
		{
			return new Vector(vector.innerMatrix * scalar);
		}

		public static Vector operator *(double scalar, Vector vector)
		{
			return vector * scalar;
		}

		public static Vector operator /(Vector vector, double scalar)
		{
			return vector * (1.0 / scalar);
		}

		public static implicit operator Matrix(Vector vector)
		{
			return vector.innerMatrix;
		}

		public static explicit operator Vector(Matrix matrix)
		{
			if (matrix.Columns != 1)
			{
				throw new InvalidCastException("Must be a row matrix.");
			}

			return new Vector(matrix);
		}

		public double this[int index]
		{
			get
			{
				return innerMatrix[index, 0];
			}
		}

		public int Length
		{
			get
			{
				return innerMatrix.Rows;
			}
		}
	}
}
