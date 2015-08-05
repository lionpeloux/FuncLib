// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Text;

namespace FuncLib.Functions
{
	[Serializable]
	[DebuggerStepThrough]
	[DebuggerDisplay("{ToString(),nq}")]
	public class Point : IPoint, IEquatable<Point>
	{
		private static Point empty;
		private Dictionary<Variable, double> assignments;
		private int hashCode;

		static Point()
		{
			empty = new Point();
		}

		private Point()
		{
			assignments = new Dictionary<Variable, double>();
			hashCode = 23;
		}

		public Point(params VariableAssignment[] assignments)
			: this()
		{
			// Copy to keep immutable (and faster using a hash table).
			foreach (VariableAssignment assignment in assignments)
			{
				Add(assignment.Variable, assignment.Value);
			}
		}

		public Point(IDictionary<Variable, double> assignments)
			: this()
		{
			foreach (KeyValuePair<Variable, double> assignment in assignments)
			{
				Add(assignment.Key, assignment.Value);
			}
		}

		private void Add(Variable variable, double value)
		{
			if (variable == null)
			{
				throw new ArgumentNullException("variable");
			}

			// Overwrite if already assigned.
			assignments[variable] = value;

			unchecked
			{
				// http://stackoverflow.com/questions/892618/create-a-hashcode-of-two-numbers
				hashCode = hashCode * 31 + variable.GetHashCode();
				hashCode = hashCode * 31 + value.GetHashCode();
			}
		}

		public bool IsAssigned(Variable variable)
		{
			return assignments.ContainsKey(variable);
		}

		public IDictionary<Variable, double> ToDictionary()
		{
			return new Dictionary<Variable, double>(assignments);
		}

		public bool Equals(Point other)
		{
			// Follows Item 6 in Effective C# in first three tests.

			if (object.ReferenceEquals(other, null))
			{
				return false;
			}

			if (object.ReferenceEquals(this, other))
			{
				return true;
			}

			if (GetType() != other.GetType())
			{
				return false;
			}

			if (hashCode != other.hashCode)
			{
				// Can't be equal if they have different hash codes.
				return false;
			}

			if (assignments.Count != other.assignments.Count)
			{
				// Can't be equal if they contain a different number of variables.
				return false;
			}

			foreach (KeyValuePair<Variable, double> assignment in assignments)
			{
				Variable variable = assignment.Key;
				double value = assignment.Value;

				double otherValue;
				if (!other.assignments.TryGetValue(variable, out otherValue))
				{
					// The other point doesn't contain this variable so they can't be equal.
					return false;
				}

				if (value != otherValue)
				{
					return false;
				}
			}

			// No differences found.
			return true;
		}

		public override bool Equals(object other)
		{
			return Equals(other as Point);
		}

		public override int GetHashCode()
		{
			return hashCode;
		}

		public override string ToString()
		{
			return ToString("r");
		}

		public string ToString(string format)
		{
			StringBuilder sb = new StringBuilder();

			foreach (KeyValuePair<Variable, double> assignment in assignments)
			{
				Variable variable = assignment.Key;
				double value = assignment.Value;

				if (sb.Length > 0)
				{
					//sb.Append(", ");
                    sb.AppendLine();
				}

				if (!string.IsNullOrEmpty(variable.Name))
				{
					sb.Append(variable.Name);
					sb.Append(" = ");
				}

				sb.Append(value.ToString(format, CultureInfo.InvariantCulture));
			}

			return sb.ToString();
		}

		public double this[Variable variable]
		{
			get
			{
				double value;
				if (!assignments.TryGetValue(variable, out value))
				{
					throw new VariableNotAssignedException(variable, "Value of a non-assigned variable is required.");
				}

				return value;
			}
		}

		public static bool operator ==(Point point1, Point point2)
		{
			if (object.ReferenceEquals(point1, point2))
			{
				return true;
			}

			if (object.ReferenceEquals(point1, null))
			{
				return false;
			}

			return point1.Equals(point2);
		}

		public static bool operator !=(Point point1, Point point2)
		{
			return !(point1 == point2);
		}

		/// <summary>
		/// Point with no variables for evaluation of constant functions.
		/// </summary>
		public static Point Empty
		{
			get
			{
				return empty;
			}
		}
	}
}
