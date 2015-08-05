// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Diagnostics;
using System.Globalization;

using FuncLib.Mathematics;

namespace FuncLib.DualNumbers.Specialized
{
	/// <summary>
	/// This class represents a function value and its first derivatives at a given point. Several standard
	/// mathematical operations are defined. Compare with <see cref="DualNumber" />. The class is immutable.
	/// </summary>
	[Serializable]
	[DebuggerStepThrough]
	public sealed class ReducedDualNumber
	{
		private static ReducedDualNumber zero;
		private double value;
		private Vector gradient; // FIXME gradientArray instead?
		private double multiplier; // gradientMultipler?

        static ReducedDualNumber()
        {
			// Used very often (e.g. in matrix initialization).
            zero = new ReducedDualNumber(0.0);
        }

        public ReducedDualNumber(double value)
			: this(value, null)
		{
		}

		public ReducedDualNumber(double value, Vector gradient)
            : this(value, gradient, 1.0)
		{
		}

        /// <summary>
        /// Using this constructor we may reuse the gradient vector instance for all unary operations.
        /// </summary>
        private ReducedDualNumber(double value, Vector gradient, double multiplier)
		{
			this.value = value;
			this.gradient = gradient;
			this.multiplier = multiplier;
		}

        /// <summary>
        /// Constructor for binary operations.
        /// </summary>
        private ReducedDualNumber(double value, Vector gradient1, double multiplier1, Vector gradient2, double multiplier2)
		{
			this.value = value;

            if (gradient1 != null && gradient2 != null)
            {
                // This is the only situation where we have to create a new gradient vector instance.
                gradient = multiplier1 * gradient1 + multiplier2 * gradient2;
                multiplier = 1.0;
            }
            else if (gradient1 != null)
            {
                gradient = gradient1;
                multiplier = multiplier1;
            }
            else if (gradient2 != null)
            {
                gradient = gradient2;
                multiplier = multiplier2;
            }
        }

		public override string ToString()
		{
			return ToString(null);
		}

		public string ToString(string format)
		{
			return value.ToString(format, CultureInfo.InvariantCulture);
		}

        /// <summary>
        /// Use the DualScalar returned by this method to represent variables as the starting points of the
        /// calculation. The gradient of a variable with respect to a set of variables is always a unit vector.
        /// </summary>
        public static ReducedDualNumber Basis(double value, int length, int index)
        {
            return new ReducedDualNumber(value, Vector.Basis(length, index));
        }

		public static ReducedDualNumber Exp(ReducedDualNumber s)
		{
			double u = Math.Exp(s.value);
			return new ReducedDualNumber(u, s.gradient, s.multiplier * u);
		}

		public static ReducedDualNumber Log(ReducedDualNumber s)
		{
			return new ReducedDualNumber(Math.Log(s.value), s.gradient, s.multiplier / s.value);
		}

		public static ReducedDualNumber Sqr(ReducedDualNumber s)
		{
			return new ReducedDualNumber(s.value * s.value, s.gradient, 2.0 * s.multiplier * s.value);
		}

		public static ReducedDualNumber Sqrt(ReducedDualNumber s)
		{
			double u = Math.Sqrt(s.value);
			return new ReducedDualNumber(u, s.gradient, 0.5 * s.multiplier / u);
		}

		public static ReducedDualNumber Pow(ReducedDualNumber s, double t)
		{
			return new ReducedDualNumber(Math.Pow(s.value, t), s.gradient, t * Math.Pow(s.value, t - 1.0) * s.multiplier);
		}

		public static ReducedDualNumber Pow(ReducedDualNumber s, ReducedDualNumber t)
		{
			double u = Math.Pow(s.value, t.value);
			return new ReducedDualNumber(u, s.gradient, t.value * Math.Pow(s.value, t.value - 1.0) * s.multiplier, t.gradient, u * Math.Log(s.value) * t.multiplier);
		}

		public static ReducedDualNumber operator +(ReducedDualNumber s, ReducedDualNumber t)
		{
			return new ReducedDualNumber(s.value + t.value, s.gradient, s.multiplier, t.gradient, t.multiplier);
		}

		public static ReducedDualNumber operator +(ReducedDualNumber s, double t)
		{
			return new ReducedDualNumber(s.value + t, s.gradient, s.multiplier);
		}

		public static ReducedDualNumber operator +(double s, ReducedDualNumber t)
		{
			return t + s;
		}

		public static ReducedDualNumber operator -(ReducedDualNumber s, ReducedDualNumber t)
		{
			return new ReducedDualNumber(s.value - t.value, s.gradient, s.multiplier, t.gradient, -t.multiplier);
		}

		public static ReducedDualNumber operator -(ReducedDualNumber s, double t)
		{
			return new ReducedDualNumber(s.value - t, s.gradient, s.multiplier);
		}

		public static ReducedDualNumber operator -(double s, ReducedDualNumber t)
		{
			return new ReducedDualNumber(s - t.value, t.gradient, -t.multiplier);
		}

		public static ReducedDualNumber operator -(ReducedDualNumber s)
		{
			return new ReducedDualNumber(-s.value, s.gradient, -s.multiplier);
		}

		public static ReducedDualNumber operator *(ReducedDualNumber s, ReducedDualNumber t)
		{
			return new ReducedDualNumber(s.value * t.value, s.gradient, t.value * s.multiplier, t.gradient, s.value * t.multiplier);
		}

		public static ReducedDualNumber operator *(ReducedDualNumber s, double t)
		{
			return new ReducedDualNumber(s.value * t, s.gradient, t * s.multiplier);
		}

		public static ReducedDualNumber operator *(double s, ReducedDualNumber t)
		{
			return t * s;
		}

		public static ReducedDualNumber operator /(ReducedDualNumber s, ReducedDualNumber t)
		{
			return new ReducedDualNumber(s.value / t.value, s.gradient, s.multiplier / t.value, t.gradient, -s.value * t.multiplier / (t.value * t.value));
		}

		public static ReducedDualNumber operator /(ReducedDualNumber s, double t)
		{
			return s * (1.0 / t);
		}

		public static ReducedDualNumber operator /(double s, ReducedDualNumber t)
		{
			return new ReducedDualNumber(s / t.value, t.gradient, -s * t.multiplier / (t.value * t.value));
		}

		public static implicit operator ReducedDualNumber(double s)
		{
            if (s == 0.0)
            {
                return zero;
            }
            
			return new ReducedDualNumber(s);
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
                if (multiplier != 1.0 && gradient != null)
                {
                    gradient *= multiplier;
                    multiplier = 1.0;
                }

                return gradient;
			}
		}
	}
}
