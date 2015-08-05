// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Diagnostics;

namespace FuncLib.Mathematics
{
	[Serializable]
	[DebuggerStepThrough]
	public struct Complex
	{
		private double re, im;

		public Complex(double re, double im)
		{
			this.re = re;
			this.im = im;
		}

		public override string ToString()
		{
			if (double.IsNaN(re) && double.IsNaN(im))
			{
				return double.NaN.ToString();
			}
			else if (im == 0.0)
			{
				return re.ToString();
			}
			else if (im > 0.0)
			{
				return string.Format("{0}+i{1}", re, im);
			}
			else
			{
				return string.Format("{0}-i{1}", re, -im);
			}
		}

		public override int GetHashCode()
		{
			return re.GetHashCode() * re.GetHashCode() + im.GetHashCode() * im.GetHashCode();
		}

		public override bool Equals(object other)
		{
			if (other is Complex)
			{
				return ((Complex)other) == this;
			}
			else
			{
				return false;
			}
		}

		public static Complex operator +(Complex z1, Complex z2)
		{
			return new Complex(z1.re + z2.re, z1.im + z2.im);
		}

		public static Complex operator +(double a, Complex z)
		{
			return new Complex(a + z.re, z.im);
		}

		public static Complex operator +(Complex z, double a)
		{
			return new Complex(z.re + a, z.im);
		}

		public static Complex operator -(Complex z)
		{
			return new Complex(-z.re, -z.im);
		}

		public static Complex operator -(Complex z1, Complex z2)
		{
			return new Complex(z1.re - z2.re, z1.im - z2.im);
		}

		public static Complex operator -(double a, Complex z)
		{
			return new Complex(a - z.re, -z.im);
		}

		public static Complex operator -(Complex z, double a)
		{
			return new Complex(z.re - a, z.im);
		}

		public static Complex operator *(Complex z1, Complex z2)
		{
			return new Complex(z1.re * z2.re - z1.im * z2.im, z1.im * z2.re + z1.re * z2.im);
		}

		public static Complex operator *(double a, Complex z)
		{
			return new Complex(a * z.re, a * z.im);
		}

		public static Complex operator *(Complex z, double a)
		{
			return new Complex(z.re * a, z.im * a);
		}

		public static Complex operator /(Complex z1, Complex z2)
		{
			double c = z2.re * z2.re + z2.im * z2.im;
			return new Complex((z1.re * z2.re + z1.im * z2.im) / c, (z1.im * z2.re - z1.re * z2.im) / c);
		}

		public static Complex operator /(double a, Complex z)
		{
			double c = z.re * z.re + z.im * z.im;
			return new Complex(a * z.re / c, -a * z.im / c);
		}

		public static Complex operator /(Complex z, double a)
		{
			return new Complex(z.re / a, z.im / a);
		}

		public static bool operator ==(Complex z1, Complex z2)
		{
			return z1.re == z2.re && z1.im == z2.im;
		}

		public static bool operator !=(Complex z1, Complex z2)
		{
			return z1.re != z2.re || z1.im != z2.im;
		}

		public static implicit operator Complex(double a)
		{
			return new Complex(a, 0.0);
		}

		public static explicit operator double(Complex z)
		{
			return z.re;
		}

		public static double Re(Complex z)
		{
			return z.re;
		}

		public static double Im(Complex z)
		{
			return z.im;
		}

		public static Complex Conjugate(Complex z)
		{
			return new Complex(z.re, -z.im);
		}

		public static double Abs(Complex z)
		{
			return Math.Sqrt(z.re * z.re + z.im * z.im);
		}

		public static Complex Exp(Complex z)
		{
			double c = Math.Exp(z.re);
			if (z.im == 0.0)
			{
				return new Complex(c, 0.0);
			}
			else
			{
				return new Complex(c * Math.Cos(z.im), c * Math.Sin(z.im));
			}
		}

		public static Complex Sqr(Complex z)
		{
			return new Complex(z.re * z.re - z.im * z.im, 2.0 * z.re * z.im);
		}

		public static Complex Sqrt(Complex z)
		{
			// Using the identity described at
			// http://en.wikipedia.org/wiki/Square_root#Square_roots_of_negative_and_complex_numbers.

			if (z.im != 0.0)
			{
				double c = Complex.Abs(z) + z.re;
				return new Complex(Math.Sqrt(0.5 * c), z.im / Math.Sqrt(2.0 * c));
			}
			else if (z.re < 0.0)
			{
				return new Complex(0.0, Math.Sqrt(-z.re));
			}
			else
			{
				return new Complex(Math.Sqrt(z.re), 0.0);
			}
		}

		public static double Arg(Complex z)
		{
			if (z.re > 0.0 || z.im != 0.0)
			{
				return 2.0 * Math.Atan2(z.im, Abs(z) + z.re);
			}
			else if (z.re < 0.0 && z.im == 0.0)
			{
				return Math.PI;
			}
			else
			{
				return double.NaN;
			}
		}

		public static Complex Log(Complex z)
		{
			return new Complex(0.5 * Math.Log(z.re * z.re + z.im * z.im), Arg(z));
		}

		public static Complex I
		{
			get
			{
				return new Complex(0.0, 1.0);
			}
		}
	}
}
