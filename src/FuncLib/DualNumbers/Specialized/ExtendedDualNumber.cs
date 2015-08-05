// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Diagnostics;

using FuncLib.Mathematics;

namespace FuncLib.DualNumbers.Specialized
{
	/// <summary>
	/// This class represents a function value and its first, second, and third partial derivatives at a given point. Several standard
	/// mathematical operations are defined. Compare with <see cref="DualNumber" />. The class is immutable.
	/// </summary>
	[Serializable]
	[DebuggerStepThrough]
	[DebuggerDisplay("{Value}")]
	public sealed class ExtendedDualNumber
	{
		private static ExtendedDualNumber zero;
		private double value;
		private int n;
		private double[] gradientArray, hessianArray, thirdArray;
		[NonSerialized]
		private Vector gradient;
		[NonSerialized]
		private Matrix hessian;
		[NonSerialized]
		private Tensor third;

		static ExtendedDualNumber()
		{
			zero = new ExtendedDualNumber(0.0);
		}

		public ExtendedDualNumber(double value)
			: this(value, null)
		{
		}

		public ExtendedDualNumber(double value, Vector gradient)
			: this(value, gradient, null)
		{
		}

		public ExtendedDualNumber(double value, Vector gradient, Matrix hessian)
		{
			this.value = value;
			this.gradient = gradient;
			this.hessian = hessian;

			if (gradient != null)
			{
				n = gradient.Length;
				gradientArray = gradient.ToArray();
			}

			if (hessian != null)
			{
				if (gradient == null)
				{
					throw new ArgumentException("The gradient must be specified if the Hessian is specified.");
				}

				if (hessian.Rows != n || hessian.Columns != n)
				{
					throw new ArgumentException("Inconsistent number of derivatives.");
				}

				// Since the Hessian is symmetric we only need to store the upper triangular part of it. Use a
				// one dimensional array until the matrix is requested (if ever). Doing it this way is almost
				// a factor of 10 faster than using naive matrix operations.

				hessianArray = new double[HessianLength(n)];

				for (int i = 0, l = 0; i < n; i++)
				{
					for (int j = i; j < n; j++, l++)
					{
						if (hessian[i, j] != hessian[j, i] && !(double.IsNaN(hessian[i, j]) && double.IsNaN(hessian[j, i])) && !(double.IsPositiveInfinity(hessian[i, j]) && double.IsPositiveInfinity(hessian[j, i])) && !(double.IsNegativeInfinity(hessian[i, j]) && double.IsNegativeInfinity(hessian[j, i])))
						{
							throw new ArgumentException("The Hessian must be symmetric.");
						}

						hessianArray[l] = hessian[i, j];
					}
				}
			}
		}

		public ExtendedDualNumber(double value, double[] gradientArray, double[] hessianArray, double[] thirdArray)
		{
			this.value = value;

			if (gradientArray != null)
			{
				n = gradientArray.Length;
				this.gradientArray = (double[])gradientArray.Clone();
			}

			if (hessianArray != null)
			{
				if (gradientArray == null)
				{
					throw new ArgumentException("The gradient must be specified if the Hessian is specified.");
				}

				if (hessianArray.Length != HessianLength(n))
				{
					throw new ArgumentException("Inconsistent number of derivatives.");
				}

				this.hessianArray = (double[])hessianArray.Clone();
			}

			if (thirdArray != null)
			{
				if (gradientArray == null || hessianArray == null)
				{
					throw new ArgumentException("The gradient and the Hessian must be specified if the third derivatives are specified.");
				}

				if (thirdArray.Length != ThirdLength(n))
				{
					throw new ArgumentException("Inconsistent number of derivatives.");
				}

				this.thirdArray = (double[])thirdArray.Clone();
			}
		}

		public ExtendedDualNumber(ExtendedDualNumber f, double g, double g1, double g11, double g111)
		{
			value = g;

			if (f.gradientArray != null)
			{
				n = f.n;
				gradientArray = new double[n];

				if (g1 != 0.0)
				{
					for (int i = 0; i < n; i++)
					{
						gradientArray[i] += g1 * f.gradientArray[i];
					}
				}

				if (f.hessianArray != null || g11 != 0.0 || g111 != 0.0)
				{
					hessianArray = new double[HessianLength(n)];

					if (g1 != 0.0 && f.hessianArray != null)
					{
						for (int i = 0, l = 0; i < n; i++)
						{
							for (int j = i; j < n; j++, l++)
							{
								hessianArray[l] += g1 * f.hessianArray[l];
							}
						}
					}

					if (g11 != 0.0)
					{
						for (int i = 0, l = 0; i < n; i++)
						{
							for (int j = i; j < n; j++, l++)
							{
								hessianArray[l] += g11 * f.gradientArray[i] * f.gradientArray[j];
							}
						}
					}

					if (f.thirdArray != null || g111 != 0.0)
					{
						thirdArray = new double[ThirdLength(n)];

						// See how the counters work in the other constructor.

						if (g1 != 0.0 && f.thirdArray != null)
						{
							for (int i = 0, m = 0; i < n; i++)
							{
								for (int j = i; j < n; j++)
								{
									for (int k = j; k < n; k++, m++)
									{
										thirdArray[m] += g1 * f.thirdArray[m];
									}
								}
							}
						}

						if (g11 != 0.0 && f.hessianArray != null)
						{
							for (int i = 0, l1 = 0, l2 = 0, m = 0; i < n; i++)
							{
								for (int j = i, l3 = l2; j < n; j++, l1++, l2 -= n - j)
								{
									for (int k = j; k < n; k++, l2++, l3++, m++)
									{
										thirdArray[m] += g11 * (f.gradientArray[i] * f.hessianArray[l3] + f.gradientArray[j] * f.hessianArray[l2] + f.gradientArray[k] * f.hessianArray[l1]);
									}
								}
							}
						}

						if (g111 != 0.0)
						{
							for (int i = 0, m = 0; i < n; i++)
							{
								for (int j = i; j < n; j++)
								{
									for (int k = j; k < n; k++, m++)
									{
										thirdArray[m] += g111 * f.gradientArray[i] * f.gradientArray[j] * f.gradientArray[k];
									}
								}
							}
						}
					}
				}
			}
		}

		public ExtendedDualNumber(ExtendedDualNumber f1, ExtendedDualNumber f2, double g, double g1, double g2, double g11, double g12, double g22, double g111, double g112, double g122, double g222)
		{
			value = g;

			if (f1.gradientArray != null || f2.gradientArray != null)
			{
				if (f1.gradientArray != null && f2.gradientArray != null && f1.n != f2.n)
				{
					throw new ArgumentException("Inconsistent number of derivatives.");
				}

				// One of the counters may be zero if the corresponding ExtendedDualNumber is a constant.
				n = Math.Max(f1.n, f2.n);
				gradientArray = new double[n];

				if (g1 != 0.0 && f1.gradientArray != null)
				{
					for (int i = 0; i < n; i++)
					{
						gradientArray[i] += g1 * f1.gradientArray[i];
					}
				}

				if (g2 != 0.0 && f2.gradientArray != null)
				{
					for (int i = 0; i < n; i++)
					{
						gradientArray[i] += g2 * f2.gradientArray[i];
					}
				}

				if (f1.hessianArray != null || f2.hessianArray != null || g11 != 0.0 || g12 != 0.0 || g22 != 0.0 || g111 != 0.0 || g112 != 0.0 || g122 != 0.0 || g222 != 0.0)
				{
					hessianArray = new double[HessianLength(n)];

					if (g1 != 0.0 && f1.hessianArray != null)
					{
						for (int i = 0, l = 0; i < n; i++)
						{
							for (int j = i; j < n; j++, l++)
							{
								hessianArray[l] += g1 * f1.hessianArray[l];
							}
						}
					}

					if (g2 != 0.0 && f2.hessianArray != null)
					{
						for (int i = 0, l = 0; i < n; i++)
						{
							for (int j = i; j < n; j++, l++)
							{
								hessianArray[l] += g2 * f2.hessianArray[l];
							}
						}
					}

					if (g11 != 0.0 && f1.gradientArray != null)
					{
						for (int i = 0, l = 0; i < n; i++)
						{
							for (int j = i; j < n; j++, l++)
							{
								hessianArray[l] += g11 * f1.gradientArray[i] * f1.gradientArray[j];
							}
						}
					}

					if (g22 != 0.0 && f2.gradientArray != null)
					{
						for (int i = 0, l = 0; i < n; i++)
						{
							for (int j = i; j < n; j++, l++)
							{
								hessianArray[l] += g22 * f2.gradientArray[i] * f2.gradientArray[j];
							}
						}
					}

					if (g12 != 0.0 && f1.gradientArray != null && f2.gradientArray != null)
					{
						for (int i = 0, l = 0; i < n; i++)
						{
							for (int j = i; j < n; j++, l++)
							{
								hessianArray[l] += g12 * (f1.gradientArray[i] * f2.gradientArray[j] + f2.gradientArray[i] * f1.gradientArray[j]);
							}
						}
					}

					if (f1.hessianArray != null || f2.hessianArray != null || g111 != 0.0 || g112 != 0.0 || g122 != 0.0 || g222 != 0.0)
					{
						thirdArray = new double[ThirdLength(n)];

						// The counters below are constructed so that:
						//
						// hessianArray[l1] is the (i,j)th   entry of the Hessian,
						// hessianArray[l2] is the (i,k)th   entry of the Hessian,
						// hessianArray[l3] is the (j,k)th   entry of the Hessian,
						// thirdArray[m]    is the (i,j,k)th entry of the third derivatives tensor.
						//
						// It's made this way to eliminate redundant entries. In case of the Hessian, only
						// the upper triangular part is stored.

						if (g1 != 0.0 && f1.thirdArray != null)
						{
							for (int i = 0, m = 0; i < n; i++)
							{
								for (int j = i; j < n; j++)
								{
									for (int k = j; k < n; k++, m++)
									{
										thirdArray[m] += g1 * f1.thirdArray[m];
									}
								}
							}
						}

						if (g2 != 0.0 && f2.thirdArray != null)
						{
							for (int i = 0, m = 0; i < n; i++)
							{
								for (int j = i; j < n; j++)
								{
									for (int k = j; k < n; k++, m++)
									{
										thirdArray[m] += g2 * f2.thirdArray[m];
									}
								}
							}
						}

						if (g11 != 0.0 && f1.hessianArray != null)
						{
							for (int i = 0, l1 = 0, l2 = 0, m = 0; i < n; i++)
							{
								for (int j = i, l3 = l2; j < n; j++, l1++, l2 -= n - j)
								{
									for (int k = j; k < n; k++, l2++, l3++, m++)
									{
										thirdArray[m] += g11 * (f1.gradientArray[i] * f1.hessianArray[l3] + f1.gradientArray[j] * f1.hessianArray[l2] + f1.gradientArray[k] * f1.hessianArray[l1]);
									}
								}
							}
						}

						if (g22 != 0.0 && f2.hessianArray != null)
						{
							for (int i = 0, l1 = 0, l2 = 0, m = 0; i < n; i++)
							{
								for (int j = i, l3 = l2; j < n; j++, l1++, l2 -= n - j)
								{
									for (int k = j; k < n; k++, l2++, l3++, m++)
									{
										thirdArray[m] += g22 * (f2.gradientArray[i] * f2.hessianArray[l3] + f2.gradientArray[j] * f2.hessianArray[l2] + f2.gradientArray[k] * f2.hessianArray[l1]);
									}
								}
							}
						}

						if (g12 != 0.0 && f1.gradientArray != null && f2.hessianArray != null)
						{
							for (int i = 0, l1 = 0, l2 = 0, m = 0; i < n; i++)
							{
								for (int j = i, l3 = l2; j < n; j++, l1++, l2 -= n - j)
								{
									for (int k = j; k < n; k++, l2++, l3++, m++)
									{
										thirdArray[m] += g12 * (f1.gradientArray[i] * f2.hessianArray[l3] + f1.gradientArray[j] * f2.hessianArray[l2] + f1.gradientArray[k] * f2.hessianArray[l1]);
									}
								}
							}
						}

						if (g12 != 0.0 && f2.gradientArray != null && f1.hessianArray != null)
						{
							for (int i = 0, l1 = 0, l2 = 0, m = 0; i < n; i++)
							{
								for (int j = i, l3 = l2; j < n; j++, l1++, l2 -= n - j)
								{
									for (int k = j; k < n; k++, l2++, l3++, m++)
									{
										thirdArray[m] += g12 * (f2.gradientArray[i] * f1.hessianArray[l3] + f2.gradientArray[j] * f1.hessianArray[l2] + f2.gradientArray[k] * f1.hessianArray[l1]);
									}
								}
							}
						}

						if (g111 != 0.0 && f1.gradientArray != null)
						{
							for (int i = 0, m = 0; i < n; i++)
							{
								for (int j = i; j < n; j++)
								{
									for (int k = j; k < n; k++, m++)
									{
										thirdArray[m] += g111 * f1.gradientArray[i] * f1.gradientArray[j] * f1.gradientArray[k];
									}
								}
							}
						}

						if (g222 != 0.0 && f2.gradientArray != null)
						{
							for (int i = 0, m = 0; i < n; i++)
							{
								for (int j = i; j < n; j++)
								{
									for (int k = j; k < n; k++, m++)
									{
										thirdArray[m] += g222 * f2.gradientArray[i] * f2.gradientArray[j] * f2.gradientArray[k];
									}
								}
							}
						}

						if (g112 != 0.0 && f1.gradientArray != null && f2.gradientArray != null)
						{
							for (int i = 0, m = 0; i < n; i++)
							{
								for (int j = i; j < n; j++)
								{
									for (int k = j; k < n; k++, m++)
									{
										thirdArray[m] += g112 * (f1.gradientArray[i] * (f1.gradientArray[j] * f2.gradientArray[k] + f2.gradientArray[j] * f1.gradientArray[k]) + f2.gradientArray[i] * f1.gradientArray[j] * f1.gradientArray[k]);
									}
								}
							}
						}

						if (g122 != 0.0 && f1.gradientArray != null && f2.gradientArray != null)
						{
							for (int i = 0, m = 0; i < n; i++)
							{
								for (int j = i; j < n; j++)
								{
									for (int k = j; k < n; k++, m++)
									{
										thirdArray[m] += g122 * (f2.gradientArray[i] * (f2.gradientArray[j] * f1.gradientArray[k] + f1.gradientArray[j] * f2.gradientArray[k]) + f1.gradientArray[i] * f2.gradientArray[j] * f2.gradientArray[k]);
									}
								}
							}
						}
					}
				}
			}
		}

		private int HessianLength(int n)
		{
			return n * (n + 1) / 2;
		}

		private int ThirdLength(int n)
		{
			return n * (n + 1) * (n + 2) / 6;
		}

		public static implicit operator ExtendedDualNumber(double f)
		{
			if (f == 0.0)
			{
				return zero;
			}

			return new ExtendedDualNumber(f);
		}

		public static ExtendedDualNumber Basis(double value, int length, int index)
		{
			return new ExtendedDualNumber(value, Vector.Basis(length, index), Matrix.Zero(length, length));
		}

		public static ExtendedDualNumber Exp(ExtendedDualNumber f)
		{
			double g = Math.Exp(f.value);
			double g1 = g;
			double g11 = g;
			double g111 = g;

			return new ExtendedDualNumber(f, g, g1, g11, g111);
		}

		public static ExtendedDualNumber Log(ExtendedDualNumber f)
		{
			double g = Math.Log(f.value);
			double g1 = 1.0 / f.value;
			double g11 = -g1 / f.value;
			double g111 = -2.0 * g11 / f.value;

			return new ExtendedDualNumber(f, g, g1, g11, g111);
		}

		public static ExtendedDualNumber Sqr(ExtendedDualNumber f)
		{
			return new ExtendedDualNumber(f, f.value * f.value, 2.0 * f.value, 2.0, 0.0);
		}

		public static ExtendedDualNumber Sqrt(ExtendedDualNumber f)
		{
			double g = Math.Sqrt(f.value);
			double g1 = 0.5 / g;
			double g11 = -0.5 * g1 / f.value;
			double g111 = -3.0 * g11 / (2.0 * f.value);

			return new ExtendedDualNumber(f, g, g1, g11, g111);
		}

		public static ExtendedDualNumber Pow(ExtendedDualNumber f1, ExtendedDualNumber f2)
		{
			/*Variable f1_ = new Variable();
			Variable f2_ = new Variable();
			Function g_ = Function.Pow(f1_, f2_);*/

			// FIXME Optimize this!

			double g = Math.Pow(f1.value, f2.value);

			double g1 = f2.value * Math.Pow(f1.value, f2.value - 1.0);
			double g2 = Math.Log(f1.value) * Math.Pow(f1.value, f2.value);

			double g11 = f2.value * (f2.value - 1.0) * Math.Pow(f1.value, f2.value - 2.0);
			double g22 = Math.Log(f1.value) * Math.Log(f1.value) * Math.Pow(f1.value, f2.value);
			double g12 = Math.Pow(f1.value, f2.value - 1.0) * (1.0 + Math.Log(f1.value) * f2.value);

			double g111 = f2.value * (f2.value - 1.0) * (f2.value - 2.0) * Math.Pow(f1.value, f2.value - 3.0);
			double g222 = Math.Log(f1.value) * Math.Log(f1.value) * Math.Log(f1.value) * Math.Pow(f1.value, f2.value);
			double g112 = (2.0 * f2.value - 1.0 + f2.value * (f2.value - 1.0) * Math.Log(f1.value)) * Math.Pow(f1.value, f2.value - 2.0);
			double g122 = 2.0 * Math.Log(f1.value) / f1.value * Math.Pow(f1.value, f2.value) + Math.Log(f1.value) * Math.Log(f1.value) * f2.value * Math.Pow(f1.value, f2.value - 1.0);

			return new ExtendedDualNumber(f1, f2, g, g1, g2, g11, g12, g22, g111, g112, g122, g222);
		}

		public static ExtendedDualNumber Cos(ExtendedDualNumber f)
		{
			double g = Math.Cos(f.value);
			double g111 = Math.Sin(f.value);
			double g1 = -g111;
			double g11 = -g;

			return new ExtendedDualNumber(f, g, g1, g11, g111);
		}

		public static ExtendedDualNumber Sin(ExtendedDualNumber f)
		{
			double g = Math.Sin(f.value);
			double g1 = Math.Cos(f.value);
			double g11 = -g;
			double g111 = -g1;

			return new ExtendedDualNumber(f, g, g1, g11, g111);
		}

		public static ExtendedDualNumber operator +(ExtendedDualNumber f1, ExtendedDualNumber f2)
		{
			return new ExtendedDualNumber(f1, f2, f1.value + f2.value, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
		}

		public static ExtendedDualNumber operator -(ExtendedDualNumber f1, ExtendedDualNumber f2)
		{
			return new ExtendedDualNumber(f1, f2, f1.value - f2.value, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
		}

		public static ExtendedDualNumber operator -(ExtendedDualNumber f)
		{
			return new ExtendedDualNumber(f, -f.value, -1.0, 0.0, 0.0);
		}

		public static ExtendedDualNumber operator *(ExtendedDualNumber f1, ExtendedDualNumber f2)
		{
			return new ExtendedDualNumber(f1, f2, f1.value * f2.value, f2.value, f1.value, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
		}

		public static ExtendedDualNumber operator /(ExtendedDualNumber f1, ExtendedDualNumber f2)
		{
			double g = f1.value / f2.value;
			double g1 = 1.0 / f2.value;
			double g2 = -g / f2.value;
			double g11 = 0.0;
			double g12 = -g1 / f2.value;
			double g22 = -2.0 * g2 / f2.value;
			double g111 = 0.0;
			double g112 = 0.0;
			double g122 = 2.0 / (f2.value * f2.value * f2.value); // FIXME optimize!
			double g222 = -6.0 * f1.value / (f2.value * f2.value * f2.value * f2.value);

			/*Variable f1_ = new Variable();
			Variable f2_ = new Variable();
			Function g_ = f1_/f2_;*/

			return new ExtendedDualNumber(f1, f2, g, g1, g2, g11, g12, g22, g111, g112, g122, g222);
		}

		public double Value
		{
			get
			{
				return value;
			}
		}

		public Vector Gradient
		{
			get
			{
				if (gradient == null && gradientArray != null)
				{
					gradient = new Vector(gradientArray);
				}

				return gradient;
			}
		}

		public Matrix Hessian
		{
			get
			{
				if (hessian == null && hessianArray != null)
				{
					double[,] a = new double[n, n];

					for (int i = 0, k = 0; i < n; i++)
					{
						for (int j = i; j < n; j++, k++)
						{
							a[i, j] = a[j, i] = hessianArray[k];
						}
					}

					hessian = new Matrix(a);
				}

				return hessian;
			}
		}

		public Tensor Third
		{
			get
			{
				if (third == null && thirdArray != null)
				{
					double[, ,] a = new double[n, n, n];

					for (int i = 0, m = 0; i < n; i++)
					{
						for (int j = i; j < n; j++)
						{
							for (int k = j; k < n; k++, m++)
							{
								a[i, j, k] = a[i, k, j] = a[j, i, k] = a[j, k, i] = a[k, i, j] = a[k, j, i] = thirdArray[m];
							}
						}
					}

					third = new Tensor(a);
				}

				return third;
			}
		}

		/*[Serializable]
		private class ThirdTensor : Tensor
		{
			private double[, ,] entries;
			private int n, n0;

			public ThirdTensor(double[, ,] entries)
			{
				this.entries = (double[, ,])entries.Clone();
			}

			//public ThirdTensor(int n, double[] thirdArray)
			//{
			//this.thirdArray=thirdArray; // copy ref.
			//}

			public override double this[int i, int j, int k]
			{
				get
				{
					if (i < 0 || i >= n0 || j < 0 || j >= n || k < 0 || k >= n)
					{
						throw new IndexOutOfRangeException();
					}
					// calculate index = i*n+...
					return entries[i, j, k];
				}
			}
		}*/
	}
}
