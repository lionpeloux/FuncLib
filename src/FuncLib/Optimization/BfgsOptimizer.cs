// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;

using FuncLib.Functions;

namespace FuncLib.Optimization
{
	/// <summary>
	/// Implements the L-BFGS-B and L-BFGS optimization methods following C. Zhu, R.H. Byrd, and J. Nocedal.
	/// Using Sergey Bochkanov's translation of Fortran code from the ALGLIB project. Only the function value
	/// itself and it's first partial derivatives are needed.
	/// </summary>
	[Serializable]
	public sealed class BfgsOptimizer : VariableConstrainedOptimizer
	{
		private double epsG, epsF, epsX;
		private int maxIterations;

		public BfgsOptimizer()
		{
			epsG = 1.0e-6;
			epsF = 1.0e-8;
			epsX = 1.0e-8;
			maxIterations = 250;
		}

		/// <summary>
		/// Run the optimizer with the specified initial values. Returns additional BFGS related convergence status.
		/// </summary>
		public BfgsOptimizerResult RunBfgs(IPoint initialPoint)
		{
			return (BfgsOptimizerResult)Run(initialPoint);
		}

		/// <summary>
		/// Run the optimizer with the specified initial values. Returns additional BFGS related convergence status.
		/// </summary>
		public BfgsOptimizerResult RunBfgs(params VariableAssignment[] initialAssignments)
		{
			return RunBfgs(new Point(initialAssignments));
		}

		public override PreparedOptimizer Prepare()
		{
			return new InnerBfgsOptimizer(this);
		}

		/// <summary>
		/// A positive number which defines a precision of search. The subroutine finishes its work if the
		/// condition ||G|| &lt; EpsG is satisfied, where ||.|| is the Euclidian norm and G is the gradient
		/// projection onto a feasible set.
		/// </summary>
		public double EpsG
		{
			get
			{
				return epsG;
			}
			set
			{
				epsG = value;
			}
		}

		/// <summary>
		/// A positive number which defines a precision of search. The subroutine finishes its work if on
		/// iteration number k+1 the condition |F(k+1)-F(k)| &lt;= EpsF*max{|F(k)|, |F(k+1)|, 1} is satisfied.
		/// </summary>
		public double EpsF
		{
			get
			{
				return epsF;
			}
			set
			{
				epsF = value;
			}
		}

		/// <summary>
		/// A positive number which defines a precision of search. The subroutine finishes its work if on
		/// iteration number k+1 the condition |X(k+1)-X(k)| &lt;= EpsX is satisfied.
		/// </summary>
		public double EpsX
		{
			get
			{
				return epsX;
			}
			set
			{
				epsX = value;
			}
		}

		/// <summary>
		/// Maximum number of iterations. If set to zero, the number of iterations is unlimited.
		/// </summary>
		public int MaxIterations
		{
			get
			{
				return maxIterations;
			}
			set
			{
				maxIterations = value;
			}
		}

		[Serializable]
		private class InnerBfgsOptimizer : InnerVariableConstrainedOptimizer
		{
			private Function function;
			private Function[] gradient;
			private int n, m;

			private double epsG, epsF, epsX;
			private int maxIterations;

			private int[] nbd;
			private double[] l, u;

			public InnerBfgsOptimizer(BfgsOptimizer optimizer)
				: base(optimizer)
			{
				function = ObjectiveFunction * ObjectiveFunctionScaling;

				// Use the maximally recommended number of corrections in the BFGS scheme.
				n = Variables.Count;
				m = Math.Min(7, n);

				gradient = new Function[n];
				for (int i = 0; i < n; i++)
				{
					gradient[i] = function.Derivative(Variables[i]);
				}

				epsG = optimizer.EpsG;
				epsF = optimizer.EpsF;
				epsX = optimizer.EpsX;
				maxIterations = optimizer.MaxIterations;

				nbd = new int[n + 1];
				l = new double[n + 1];
				u = new double[n + 1];
				foreach (VariableConstraint constraint in VariableConstraints)
				{
					int i = Variables.IndexOf(constraint.Variable) + 1;
					double minValue = constraint.MinValue;
					double maxValue = constraint.MaxValue;

					if (!double.IsNegativeInfinity(minValue))
					{
						// Has lower boundary.
						nbd[i] = 1;
						l[i] = minValue;
					}

					if (!double.IsPositiveInfinity(maxValue))
					{
						// Has upper boundary.
						nbd[i] = 3;
						u[i] = maxValue;
					}

					if (!double.IsNegativeInfinity(minValue) && !double.IsPositiveInfinity(maxValue))
					{
						// Has both lower and upper boundaries (values are already defined).
						nbd[i] = 2;
					}
				}
			}

			public override IOptimizerResult Run(IPoint initialPoint)
			{
				CheckConstraints(initialPoint);

				double[] x = new double[n + 1];
				for (int i = 0; i < n; i++)
				{
					x[i + 1] = initialPoint[Variables[i]];
				}

				int info = 0;
				if (VariableConstraints.Count > 0)
				{
					// Use L-BFGS-B.
					lbfgsbminimize(n, m, ref x, epsG, epsF, epsX, maxIterations, ref nbd, ref l, ref u, ref info);
				}
				else
				{
					// Use L-BFGS.
					lbfgsminimize(n, m, ref x, epsG, epsF, epsX, maxIterations, ref info);
				}

				IPoint optimalPoint = CreateOptimalPoint(Transform(x));
				double optimalValue = ObjectiveFunction.Value(CreateEvaluation(Transform(x)));

				return new BfgsOptimizerResult(info > 0, optimalPoint, optimalValue, (BfgsOptimizerConvergenceStatus)info);
			}

			private double[] Transform(double[] x)
			{
				// Transform to usual zero-based indexing and create evaluation object.
				double[] values = new double[n];
				for (int i = 0; i < n; i++)
				{
					values[i] = x[i + 1];
				}
				return values;
			}

			/*
			This members must be defined by you:
			static void funcgrad(ref double[] x,
				ref double f,
				ref double[] g)
			*/

			private void funcgrad(ref double[] x,
				ref double f,
				ref double[] g)
			{
				IEvaluation evaluation = CreateEvaluation(Transform(x));

				f = function.Value(evaluation);
				for (int i = 0; i < n; i++)
				{
					g[i + 1] = gradient[i].Value(evaluation);
				}
			}

			private const double MachineEpsilon = 5e-16;

			private double Sqr(double x)
			{
				return x * x;
			}

			/*************************************************************************
			NEOS, November 1994. (Latest revision June 1996.)
			Optimization Technology Center.
			Argonne National Laboratory and Northwestern University.

			Written by Ciyou Zhu in collaboration with
			R.H. Byrd, P. Lu-Chen and J. Nocedal.

			Contributors:
				* Sergey Bochkanov (ALGLIB project). Translation from FORTRAN to
				  pseudocode.
		      
			This software is freely available, but we  expect  that  all  publications
			describing  work using this software, or all commercial products using it,
			quote at least one of the references given below:
				* R. H. Byrd, P. Lu and J. Nocedal.  A Limited  Memory  Algorithm  for
				  Bound Constrained Optimization, (1995), SIAM Journal  on  Scientific
				  and Statistical Computing , 16, 5, pp. 1190-1208.
				* C. Zhu, R.H. Byrd and J. Nocedal. L-BFGS-B: Algorithm 778: L-BFGS-B,
				  FORTRAN routines for  large  scale  bound  constrained  optimization
				  (1997), ACM Transactions on Mathematical Software,  Vol 23,  Num. 4,
				  pp. 550 - 560.
			*************************************************************************/

			/*************************************************************************
			The  subroutine  minimizes  the  function  F(x) of N arguments with simple
			constraints using a quasi-Newton method (LBFGS scheme) which is  optimized
			to use a minimum amount of memory.

			The subroutine generates the approximation of an inverse Hessian matrix by
			using information about the last M steps of the algorithm (instead  of N).
			It lessens a required amount of memory from a value  of  order  N^2  to  a
			value of order 2*N*M.

			This subroutine uses the FuncGrad subroutine which calculates the value of
			the function F and gradient G in point X. The programmer should define the
			FuncGrad subroutine by himself.  It should be noted  that  the  subroutine
			doesn't need to waste  time for memory allocation of array G, because  the
			memory is allocated in calling the  subroutine.  Setting  a  dimension  of
			array G each time when calling a subroutine will excessively slow down  an
			algorithm.

			The programmer could also redefine the LBFGSNewIteration subroutine  which
			is called on each new step. The current point X, the function value F  and
			the gradient G are passed  into  this  subroutine.  It  is  reasonable  to
			redefine the subroutine for better debugging, for  example,  to  visualize
			the solution process.

			Input parameters:
				N       -   problem dimension. N>0
				M       -   number of  corrections  in  the  BFGS  scheme  of  Hessian
							approximation  update.  Recommended value:  3<=M<=7.   The
							smaller value causes worse convergence,  the  bigger  will
							not  cause  a  considerably  better  convergence, but will
							cause a fall in the performance. M<=N.
				X       -   initial solution approximation.
							Array whose index ranges from 1 to N.
				EpsG    -   positive number which defines a precision of  search.  The
							subroutine finishes its work if the condition ||G|| < EpsG
							is satisfied, where ||.|| means Euclidian norm, G - gradient
							projection onto a feasible set, X - current approximation.
				EpsF    -   positive number which defines a precision of  search.  The
							subroutine  finishes  its  work if on iteration number k+1
							the condition |F(k+1)-F(k)| <= EpsF*max{|F(k)|, |F(k+1)|, 1}
							is satisfied.
				EpsX    -   positive number which defines a precision of  search.  The
							subroutine  finishes  its  work if on iteration number k+1
							the condition |X(k+1)-X(k)| <= EpsX is satisfied.
				MaxIts  -   maximum number of iterations.
							If MaxIts=0, the number of iterations is unlimited.
				NBD     -   constraint type. If NBD(i) is equal to:
							* 0, X(i) has no constraints,
							* 1, X(i) has only lower boundary,
							* 2, X(i) has both lower and upper boundaries,
							* 3, X(i) has only upper boundary,
							Array whose index ranges from 1 to N.
				L       -   lower boundaries of X(i) variables.
							Array whose index ranges from 1 to N.
				U       -   upper boundaries of X(i) variables.
							Array whose index ranges from 1 to N.

			Output parameters:
				X       -   solution approximation.
			Array whose index ranges from 1 to N.
				Info    -   a return code:
								* -2 unknown internal error,
								* -1 wrong parameters were specified,
								* 0 interrupted by user,
								* 1 relative function decreasing is less or equal to EpsF,
								* 2 step is less or equal to EpsX,
								* 4 gradient norm is less or equal to EpsG,
								* 5 number of iterations exceeds MaxIts.

			FuncGrad routine description. User-defined.
			Input parameters:
				X   -   array whose index ranges from 1 to N.
			Output parameters:
				F   -   function value at X.
				G   -   function gradient.
						Array whose index ranges from 1 to N.
			The memory for array G has already been allocated in the calling subroutine,
			and it isn't necessary to allocate it in the FuncGrad subroutine.

				NEOS, November 1994. (Latest revision June 1996.)
				Optimization Technology Center.
				Argonne National Laboratory and Northwestern University.

				Written by Ciyou Zhu in collaboration with
				R.H. Byrd, P. Lu-Chen and J. Nocedal.
			*************************************************************************/

			private void lbfgsbminimize(int n,
				int m,
				ref double[] x,
				double epsg,
				double epsf,
				double epsx,
				int maxits,
				ref int[] nbd,
				ref double[] l,
				ref double[] u,
				ref int info)
			{
				double f = 0;
				double[] g = new double[0];
				double[] xold = new double[0];
				double[] xdiff = new double[0];
				double[,] ws = new double[0, 0];
				double[,] wy = new double[0, 0];
				double[,] sy = new double[0, 0];
				double[,] ss = new double[0, 0];
				double[,] yy = new double[0, 0];
				double[,] wt = new double[0, 0];
				double[,] wn = new double[0, 0];
				double[,] snd = new double[0, 0];
				double[] z = new double[0];
				double[] r = new double[0];
				double[] d = new double[0];
				double[] t = new double[0];
				double[] wa = new double[0];
				double[] sg = new double[0];
				double[] sgo = new double[0];
				double[] yg = new double[0];
				double[] ygo = new double[0];
				int[] index = new int[0];
				int[] iwhere = new int[0];
				int[] indx2 = new int[0];
				int csave = 0;
				bool[] lsave = new bool[0];
				int[] isave = new int[0];
				double[] dsave = new double[0];
				int task = 0;
				bool prjctd = new bool();
				bool cnstnd = new bool();
				bool boxed = new bool();
				bool updatd = new bool();
				bool wrk = new bool();
				int i = 0;
				int k = 0;
				int nintol = 0;
				int iback = 0;
				int nskip = 0;
				int head = 0;
				int col = 0;
				int iter = 0;
				int itail = 0;
				int iupdat = 0;
				int nint = 0;
				int nfgv = 0;
				int internalinfo = 0;
				int ifun = 0;
				int iword = 0;
				int nfree = 0;
				int nact = 0;
				int ileave = 0;
				int nenter = 0;
				double theta = 0;
				double fold = 0;
				double dr = 0;
				double rr = 0;
				double dnrm = 0;
				double xstep = 0;
				double sbgnrm = 0;
				double ddum = 0;
				double dtd = 0;
				double gd = 0;
				double gdold = 0;
				double stp = 0;
				double stpmx = 0;
				double tf = 0;
				double[] workvec = new double[0];
				double[] workvec2 = new double[0];
				double[] dsave13 = new double[0];
				double[] wa0 = new double[0];
				double[] wa1 = new double[0];
				double[] wa2 = new double[0];
				double[] wa3 = new double[0];
				double[,] workmat = new double[0, 0];
				int[] isave2 = new int[0];
				int i_ = 0;
				int i1_ = 0;

				workvec = new double[m + 1];
				workvec2 = new double[2 * m + 1];
				workmat = new double[m + 1, m + 1];
				isave2 = new int[2 + 1];
				dsave13 = new double[13 + 1];
				wa0 = new double[2 * m + 1];
				wa1 = new double[2 * m + 1];
				wa2 = new double[2 * m + 1];
				wa3 = new double[2 * m + 1];
				g = new double[n + 1];
				xold = new double[n + 1];
				xdiff = new double[n + 1];
				ws = new double[n + 1, m + 1];
				wy = new double[n + 1, m + 1];
				sy = new double[m + 1, m + 1];
				ss = new double[m + 1, m + 1];
				yy = new double[m + 1, m + 1];
				wt = new double[m + 1, m + 1];
				wn = new double[2 * m + 1, 2 * m + 1];
				snd = new double[2 * m + 1, 2 * m + 1];
				z = new double[n + 1];
				r = new double[n + 1];
				d = new double[n + 1];
				t = new double[n + 1];
				wa = new double[8 * m + 1];
				sg = new double[m + 1];
				sgo = new double[m + 1];
				yg = new double[m + 1];
				ygo = new double[m + 1];
				index = new int[n + 1];
				iwhere = new int[n + 1];
				indx2 = new int[n + 1];
				lsave = new bool[4 + 1];
				isave = new int[23 + 1];
				dsave = new double[29 + 1];
				col = 0;
				head = 1;
				theta = 1;
				iupdat = 0;
				updatd = false;
				iter = 0;
				nfgv = 0;
				nint = 0;
				nintol = 0;
				nskip = 0;
				nfree = n;
				internalinfo = 0;
				lbfgsberrclb(n, m, epsf, ref l, ref u, ref nbd, ref task, ref internalinfo, ref k);
				if (task == 2 | maxits < 0 | epsg < 0 | epsx < 0)
				{
					info = -1;
					return;
				}
				lbfgsbactive(n, ref l, ref u, ref nbd, ref x, ref iwhere, ref prjctd, ref cnstnd, ref boxed);
				for (i_ = 1; i_ <= n; i_++)
				{
					xold[i_] = x[i_];
				}
				funcgrad(ref x, ref f, ref g);
				nfgv = 1;
				lbfgsbprojgr(n, ref l, ref u, ref nbd, ref x, ref g, ref sbgnrm);
				if (sbgnrm <= epsg)
				{
					info = 4;
					return;
				}
				while (true)
				{
					iword = -1;
					if (!cnstnd & col > 0)
					{
						for (i_ = 1; i_ <= n; i_++)
						{
							z[i_] = x[i_];
						}
						wrk = updatd;
						nint = 0;
					}
					else
					{
						for (i_ = 1; i_ <= 2 * m; i_++)
						{
							wa0[i_] = wa[i_];
						}
						i1_ = (2 * m + 1) - (1);
						for (i_ = 1; i_ <= 2 * m; i_++)
						{
							wa1[i_] = wa[i_ + i1_];
						}
						i1_ = (4 * m + 1) - (1);
						for (i_ = 1; i_ <= 2 * m; i_++)
						{
							wa2[i_] = wa[i_ + i1_];
						}
						i1_ = (6 * m + 1) - (1);
						for (i_ = 1; i_ <= 2 * m; i_++)
						{
							wa3[i_] = wa[i_ + i1_];
						}
						lbfgsbcauchy(n, ref x, ref l, ref u, ref nbd, ref g, ref indx2, ref iwhere, ref t, ref d, ref z, m, ref wy, ref ws, ref sy, ref wt, theta, col, head, ref wa0, ref wa1, ref wa2, ref wa3, ref nint, ref sg, ref yg, sbgnrm, ref internalinfo, ref workvec);
						for (i_ = 1; i_ <= 2 * m; i_++)
						{
							wa[i_] = wa0[i_];
						}
						i1_ = (1) - (2 * m + 1);
						for (i_ = 2 * m + 1; i_ <= 4 * m; i_++)
						{
							wa[i_] = wa1[i_ + i1_];
						}
						i1_ = (1) - (4 * m + 1);
						for (i_ = 4 * m + 1; i_ <= 6 * m; i_++)
						{
							wa[i_] = wa2[i_ + i1_];
						}
						i1_ = (1) - (6 * m + 1);
						for (i_ = 6 * m + 1; i_ <= 8 * m; i_++)
						{
							wa[i_] = wa3[i_ + i1_];
						}
						if (internalinfo != 0)
						{
							internalinfo = 0;
							col = 0;
							head = 1;
							theta = 1;
							iupdat = 0;
							updatd = false;
							continue;
						}
						nintol = nintol + nint;
						lbfgsbfreev(n, ref nfree, ref index, ref nenter, ref ileave, ref indx2, ref iwhere, ref wrk, updatd, cnstnd, iter);
						nact = n - nfree;
					}
					if (nfree != 0 & col != 0)
					{
						if (wrk)
						{
							lbfgsbformk(n, nfree, ref index, nenter, ileave, ref indx2, iupdat, updatd, ref wn, ref snd, m, ref ws, ref wy, ref sy, theta, col, head, ref internalinfo, ref workvec, ref workmat);
						}
						if (internalinfo != 0)
						{
							internalinfo = 0;
							col = 0;
							head = 1;
							theta = 1;
							iupdat = 0;
							updatd = false;
							continue;
						}
						lbfgsbcmprlb(n, m, ref x, ref g, ref ws, ref wy, ref sy, ref wt, ref z, ref r, ref wa, ref index, theta, col, head, nfree, cnstnd, ref internalinfo, ref workvec, ref workvec2);
						if (internalinfo == 0)
						{
							lbfgsbsubsm(n, m, nfree, ref index, ref l, ref u, ref nbd, ref z, ref r, ref ws, ref wy, theta, col, head, ref iword, ref wa, ref wn, ref internalinfo);
						}
						if (internalinfo != 0)
						{
							internalinfo = 0;
							col = 0;
							head = 1;
							theta = 1;
							iupdat = 0;
							updatd = false;
							continue;
						}
					}
					for (i = 1; i <= n; i++)
					{
						d[i] = z[i] - x[i];
					}
					task = 0;
					while (true)
					{
						lbfgsblnsrlb(n, ref l, ref u, ref nbd, ref x, f, ref fold, ref gd, ref gdold, ref g, ref d, ref r, ref t, ref z, ref stp, ref dnrm, ref dtd, ref xstep, ref stpmx, iter, ref ifun, ref iback, ref nfgv, ref internalinfo, ref task, boxed, cnstnd, ref csave, ref isave2, ref dsave13);
						if (internalinfo != 0 | iback >= 20 | task != 1)
						{
							break;
						}
						funcgrad(ref x, ref f, ref g);
					}
					if (internalinfo != 0)
					{
						for (i_ = 1; i_ <= n; i_++)
						{
							x[i_] = t[i_];
						}
						for (i_ = 1; i_ <= n; i_++)
						{
							g[i_] = r[i_];
						}
						f = fold;
						if (col == 0)
						{
							if (internalinfo == 0)
							{
								internalinfo = -9;
								nfgv = nfgv - 1;
								ifun = ifun - 1;
								iback = iback - 1;
							}
							task = 2;
							iter = iter + 1;
							info = -2;
							return;
						}
						else
						{
							if (internalinfo == 0)
							{
								nfgv = nfgv - 1;
							}
							internalinfo = 0;
							col = 0;
							head = 1;
							theta = 1;
							iupdat = 0;
							updatd = false;
							continue;
						}
					}
					iter = iter + 1;
					lbfgsbnewiteration(ref x, f, ref g);
					lbfgsbprojgr(n, ref l, ref u, ref nbd, ref x, ref g, ref sbgnrm);
					if (sbgnrm <= epsg)
					{
						info = 4;
						return;
					}
					for (i_ = 1; i_ <= n; i_++)
					{
						xdiff[i_] = xold[i_];
					}
					for (i_ = 1; i_ <= n; i_++)
					{
						xdiff[i_] = xdiff[i_] - x[i_];
					}
					tf = 0.0;
					for (i_ = 1; i_ <= n; i_++)
					{
						tf += xdiff[i_] * xdiff[i_];
					}
					tf = Math.Sqrt(tf);
					if (tf <= epsx)
					{
						info = 2;
						return;
					}
					ddum = Math.Max(Math.Abs(fold), Math.Max(Math.Abs(f), 1));
					if (fold - f <= epsf * ddum)
					{
						info = 1;
						return;
					}
					if (iter > maxits & maxits > 0)
					{
						info = 5;
						return;
					}
					if (additionallbfgsbstoppingcriterion(iter, ref x, f, ref g))
					{
						info = 0;
						return;
					}
					for (i_ = 1; i_ <= n; i_++)
					{
						xold[i_] = x[i_];
					}
					for (i = 1; i <= n; i++)
					{
						r[i] = g[i] - r[i];
					}
					rr = 0.0;
					for (i_ = 1; i_ <= n; i_++)
					{
						rr += r[i_] * r[i_];
					}
					if (stp == 1)
					{
						dr = gd - gdold;
						ddum = -gdold;
					}
					else
					{
						dr = (gd - gdold) * stp;
						for (i_ = 1; i_ <= n; i_++)
						{
							d[i_] = stp * d[i_];
						}
						ddum = -(gdold * stp);
					}
					if (dr <= MachineEpsilon * ddum)
					{
						nskip = nskip + 1;
						updatd = false;
					}
					else
					{
						updatd = true;
						iupdat = iupdat + 1;
						lbfgsbmatupd(n, m, ref ws, ref wy, ref sy, ref ss, ref d, ref r, ref itail, iupdat, ref col, ref head, ref theta, rr, dr, stp, dtd);
						lbfgsbformt(m, ref wt, ref sy, ref ss, col, theta, ref internalinfo);
						if (internalinfo != 0)
						{
							internalinfo = 0;
							col = 0;
							head = 1;
							theta = 1;
							iupdat = 0;
							updatd = false;
							continue;
						}
					}
				}
			}

			private void lbfgsbactive(int n,
				ref double[] l,
				ref double[] u,
				ref int[] nbd,
				ref double[] x,
				ref int[] iwhere,
				ref bool prjctd,
				ref bool cnstnd,
				ref bool boxed)
			{
				int nbdd = 0;
				int i = 0;

				nbdd = 0;
				prjctd = false;
				cnstnd = false;
				boxed = true;
				for (i = 1; i <= n; i++)
				{
					if (nbd[i] > 0)
					{
						if (nbd[i] <= 2 & x[i] <= l[i])
						{
							if (x[i] < l[i])
							{
								prjctd = true;
								x[i] = l[i];
							}
							nbdd = nbdd + 1;
						}
						else
						{
							if (nbd[i] >= 2 & x[i] >= u[i])
							{
								if (x[i] > u[i])
								{
									prjctd = true;
									x[i] = u[i];
								}
								nbdd = nbdd + 1;
							}
						}
					}
				}
				for (i = 1; i <= n; i++)
				{
					if (nbd[i] != 2)
					{
						boxed = false;
					}
					if (nbd[i] == 0)
					{
						iwhere[i] = -1;
					}
					else
					{
						cnstnd = true;
						if (nbd[i] == 2 & u[i] - l[i] <= 0)
						{
							iwhere[i] = 3;
						}
						else
						{
							iwhere[i] = 0;
						}
					}
				}
			}

			private void lbfgsbbmv(int m,
				ref double[,] sy,
				ref double[,] wt,
				int col,
				ref double[] v,
				ref double[] p,
				ref int info,
				ref double[] workvec)
			{
				int i = 0;
				int k = 0;
				int i2 = 0;
				double s = 0;
				int i_ = 0;
				int i1_ = 0;

				if (col == 0)
				{
					return;
				}
				p[col + 1] = v[col + 1];
				for (i = 2; i <= col; i++)
				{
					i2 = col + i;
					s = 0.0;
					for (k = 1; k <= i - 1; k++)
					{
						s = s + sy[i, k] * v[k] / sy[k, k];
					}
					p[i2] = v[i2] + s;
				}
				i1_ = (col + 1) - (1);
				for (i_ = 1; i_ <= col; i_++)
				{
					workvec[i_] = p[i_ + i1_];
				}
				lbfgsbdtrsl(ref wt, col, ref workvec, 11, ref info);
				i1_ = (1) - (col + 1);
				for (i_ = col + 1; i_ <= col + col; i_++)
				{
					p[i_] = workvec[i_ + i1_];
				}
				if (info != 0)
				{
					return;
				}
				for (i = 1; i <= col; i++)
				{
					p[i] = v[i] / Math.Sqrt(sy[i, i]);
				}
				i1_ = (col + 1) - (1);
				for (i_ = 1; i_ <= col; i_++)
				{
					workvec[i_] = p[i_ + i1_];
				}
				lbfgsbdtrsl(ref wt, col, ref workvec, 1, ref info);
				i1_ = (1) - (col + 1);
				for (i_ = col + 1; i_ <= col + col; i_++)
				{
					p[i_] = workvec[i_ + i1_];
				}
				if (info != 0)
				{
					return;
				}
				for (i = 1; i <= col; i++)
				{
					p[i] = -(p[i] / Math.Sqrt(sy[i, i]));
				}
				for (i = 1; i <= col; i++)
				{
					s = 0;
					for (k = i + 1; k <= col; k++)
					{
						s = s + sy[k, i] * p[col + k] / sy[i, i];
					}
					p[i] = p[i] + s;
				}
			}

			private void lbfgsbcauchy(int n,
				ref double[] x,
				ref double[] l,
				ref double[] u,
				ref int[] nbd,
				ref double[] g,
				ref int[] iorder,
				ref int[] iwhere,
				ref double[] t,
				ref double[] d,
				ref double[] xcp,
				int m,
				ref double[,] wy,
				ref double[,] ws,
				ref double[,] sy,
				ref double[,] wt,
				double theta,
				int col,
				int head,
				ref double[] p,
				ref double[] c,
				ref double[] wbp,
				ref double[] v,
				ref int nint,
				ref double[] sg,
				ref double[] yg,
				double sbgnrm,
				ref int info,
				ref double[] workvec)
			{
				bool xlower = new bool();
				bool xupper = new bool();
				bool bnded = new bool();
				int i = 0;
				int j = 0;
				int col2 = 0;
				int nfree = 0;
				int nbreak = 0;
				int pointr = 0;
				int ibp = 0;
				int nleft = 0;
				int ibkmin = 0;
				int iter = 0;
				double f1 = 0;
				double f2 = 0;
				double dt = 0;
				double dtm = 0;
				double tsum = 0;
				double dibp = 0;
				double zibp = 0;
				double dibp2 = 0;
				double bkmin = 0;
				double tu = 0;
				double tl = 0;
				double wmc = 0;
				double wmp = 0;
				double wmw = 0;
				double tj = 0;
				double tj0 = 0;
				double neggi = 0;
				double f2org = 0;
				double tmpv = 0;
				int i_ = 0;

				if (sbgnrm <= 0)
				{
					for (i_ = 1; i_ <= n; i_++)
					{
						xcp[i_] = x[i_];
					}
					return;
				}
				bnded = true;
				nfree = n + 1;
				nbreak = 0;
				ibkmin = 0;
				bkmin = 0;
				col2 = 2 * col;
				f1 = 0;
				for (i = 1; i <= col2; i++)
				{
					p[i] = 0;
				}
				for (i = 1; i <= n; i++)
				{
					neggi = -g[i];
					if (iwhere[i] != 3 & iwhere[i] != -1)
					{
						tl = 0;
						tu = 0;
						if (nbd[i] <= 2)
						{
							tl = x[i] - l[i];
						}
						if (nbd[i] >= 2)
						{
							tu = u[i] - x[i];
						}
						xlower = nbd[i] <= 2 & tl <= 0;
						xupper = nbd[i] >= 2 & tu <= 0;
						iwhere[i] = 0;
						if (xlower)
						{
							if (neggi <= 0)
							{
								iwhere[i] = 1;
							}
						}
						else
						{
							if (xupper)
							{
								if (neggi >= 0)
								{
									iwhere[i] = 2;
								}
							}
							else
							{
								if (Math.Abs(neggi) <= 0)
								{
									iwhere[i] = -3;
								}
							}
						}
					}
					pointr = head;
					if (iwhere[i] != 0 & iwhere[i] != -1)
					{
						d[i] = 0;
					}
					else
					{
						d[i] = neggi;
						f1 = f1 - neggi * neggi;
						for (j = 1; j <= col; j++)
						{
							p[j] = p[j] + wy[i, pointr] * neggi;
							p[col + j] = p[col + j] + ws[i, pointr] * neggi;
							pointr = pointr % m + 1;
						}
						if (nbd[i] <= 2 & nbd[i] != 0 & neggi < 0)
						{
							nbreak = nbreak + 1;
							iorder[nbreak] = i;
							t[nbreak] = tl / -neggi;
							if (nbreak == 1 | t[nbreak] < bkmin)
							{
								bkmin = t[nbreak];
								ibkmin = nbreak;
							}
						}
						else
						{
							if (nbd[i] >= 2 & neggi > 0)
							{
								nbreak = nbreak + 1;
								iorder[nbreak] = i;
								t[nbreak] = tu / neggi;
								if (nbreak == 1 | t[nbreak] < bkmin)
								{
									bkmin = t[nbreak];
									ibkmin = nbreak;
								}
							}
							else
							{
								nfree = nfree - 1;
								iorder[nfree] = i;
								if (Math.Abs(neggi) > 0)
								{
									bnded = false;
								}
							}
						}
					}
				}
				if (theta != 1)
				{
					for (i_ = col + 1; i_ <= col + col; i_++)
					{
						p[i_] = theta * p[i_];
					}
				}
				for (i_ = 1; i_ <= n; i_++)
				{
					xcp[i_] = x[i_];
				}
				if (nbreak == 0 & nfree == n + 1)
				{
					return;
				}
				for (j = 1; j <= col2; j++)
				{
					c[j] = 0;
				}
				f2 = -(theta * f1);
				f2org = f2;
				if (col > 0)
				{
					lbfgsbbmv(m, ref sy, ref wt, col, ref p, ref v, ref info, ref workvec);
					if (info != 0)
					{
						return;
					}
					tmpv = 0.0;
					for (i_ = 1; i_ <= col2; i_++)
					{
						tmpv += v[i_] * p[i_];
					}
					f2 = f2 - tmpv;
				}
				dtm = -(f1 / f2);
				tsum = 0;
				nint = 1;
				if (nbreak != 0)
				{
					nleft = nbreak;
					iter = 1;
					tj = 0;
					while (true)
					{
						tj0 = tj;
						if (iter == 1)
						{
							tj = bkmin;
							ibp = iorder[ibkmin];
						}
						else
						{
							if (iter == 2)
							{
								if (ibkmin != nbreak)
								{
									t[ibkmin] = t[nbreak];
									iorder[ibkmin] = iorder[nbreak];
								}
							}
							lbfgsbhpsolb(nleft, ref t, ref iorder, iter - 2);
							tj = t[nleft];
							ibp = iorder[nleft];
						}
						dt = tj - tj0;
						if (dtm < dt)
						{
							break;
						}
						tsum = tsum + dt;
						nleft = nleft - 1;
						iter = iter + 1;
						dibp = d[ibp];
						d[ibp] = 0;
						if (dibp > 0)
						{
							zibp = u[ibp] - x[ibp];
							xcp[ibp] = u[ibp];
							iwhere[ibp] = 2;
						}
						else
						{
							zibp = l[ibp] - x[ibp];
							xcp[ibp] = l[ibp];
							iwhere[ibp] = 1;
						}
						if (nleft == 0 & nbreak == n)
						{
							dtm = dt;
							if (col > 0)
							{
								for (i_ = 1; i_ <= col2; i_++)
								{
									c[i_] = c[i_] + dtm * p[i_];
								}
							}
							return;
						}
						nint = nint + 1;
						dibp2 = Sqr(dibp);
						f1 = f1 + dt * f2 + dibp2 - theta * dibp * zibp;
						f2 = f2 - theta * dibp2;
						if (col > 0)
						{
							for (i_ = 1; i_ <= col2; i_++)
							{
								c[i_] = c[i_] + dt * p[i_];
							}
							pointr = head;
							for (j = 1; j <= col; j++)
							{
								wbp[j] = wy[ibp, pointr];
								wbp[col + j] = theta * ws[ibp, pointr];
								pointr = pointr % m + 1;
							}
							lbfgsbbmv(m, ref sy, ref wt, col, ref wbp, ref v, ref info, ref workvec);
							if (info != 0)
							{
								return;
							}
							wmc = 0.0;
							for (i_ = 1; i_ <= col2; i_++)
							{
								wmc += c[i_] * v[i_];
							}
							wmp = 0.0;
							for (i_ = 1; i_ <= col2; i_++)
							{
								wmp += p[i_] * v[i_];
							}
							wmw = 0.0;
							for (i_ = 1; i_ <= col2; i_++)
							{
								wmw += wbp[i_] * v[i_];
							}
							for (i_ = 1; i_ <= col2; i_++)
							{
								p[i_] = p[i_] - dibp * wbp[i_];
							}
							f1 = f1 + dibp * wmc;
							f2 = f2 + 2.0 * dibp * wmp - dibp2 * wmw;
						}
						f2 = Math.Max(MachineEpsilon * f2org, f2);
						if (nleft > 0)
						{
							dtm = -(f1 / f2);
							continue;
						}
						else
						{
							if (bnded)
							{
								f1 = 0;
								f2 = 0;
								dtm = 0;
							}
							else
							{
								dtm = -(f1 / f2);
							}
						}
						break;
					}
				}
				if (dtm <= 0)
				{
					dtm = 0;
				}
				tsum = tsum + dtm;
				for (i_ = 1; i_ <= n; i_++)
				{
					xcp[i_] = xcp[i_] + tsum * d[i_];
				}
				if (col > 0)
				{
					for (i_ = 1; i_ <= col2; i_++)
					{
						c[i_] = c[i_] + dtm * p[i_];
					}
				}
			}

			private void lbfgsbcmprlb(int n,
				int m,
				ref double[] x,
				ref double[] g,
				ref double[,] ws,
				ref double[,] wy,
				ref double[,] sy,
				ref double[,] wt,
				ref double[] z,
				ref double[] r,
				ref double[] wa,
				ref int[] index,
				double theta,
				int col,
				int head,
				int nfree,
				bool cnstnd,
				ref int info,
				ref double[] workvec,
				ref double[] workvec2)
			{
				int i = 0;
				int j = 0;
				int k = 0;
				int pointr = 0;
				double a1 = 0;
				double a2 = 0;
				int i_ = 0;
				int i1_ = 0;

				if (!cnstnd & col > 0)
				{
					for (i = 1; i <= n; i++)
					{
						r[i] = -g[i];
					}
				}
				else
				{
					for (i = 1; i <= nfree; i++)
					{
						k = index[i];
						r[i] = -(theta * (z[k] - x[k])) - g[k];
					}
					i1_ = (2 * m + 1) - (1);
					for (i_ = 1; i_ <= 2 * m; i_++)
					{
						workvec2[i_] = wa[i_ + i1_];
					}
					lbfgsbbmv(m, ref sy, ref wt, col, ref workvec2, ref wa, ref info, ref workvec);
					i1_ = (1) - (2 * m + 1);
					for (i_ = 2 * m + 1; i_ <= 4 * m; i_++)
					{
						wa[i_] = workvec2[i_ + i1_];
					}
					if (info != 0)
					{
						info = -8;
						return;
					}
					pointr = head;
					for (j = 1; j <= col; j++)
					{
						a1 = wa[j];
						a2 = theta * wa[col + j];
						for (i = 1; i <= nfree; i++)
						{
							k = index[i];
							r[i] = r[i] + wy[k, pointr] * a1 + ws[k, pointr] * a2;
						}
						pointr = pointr % m + 1;
					}
				}
			}

			private void lbfgsberrclb(int n,
				int m,
				double factr,
				ref double[] l,
				ref double[] u,
				ref int[] nbd,
				ref int task,
				ref int info,
				ref int k)
			{
				int i = 0;

				if (n <= 0)
				{
					task = 2;
				}
				if (m <= 0)
				{
					task = 2;
				}
				if (m > n)
				{
					task = 2;
				}
				if (factr < 0)
				{
					task = 2;
				}
				for (i = 1; i <= n; i++)
				{
					if (nbd[i] < 0 | nbd[i] > 3)
					{
						task = 2;
						info = -6;
						k = i;
					}
					if (nbd[i] == 2)
					{
						if (l[i] > u[i])
						{
							task = 2;
							info = -7;
							k = i;
						}
					}
				}
			}

			private void lbfgsbformk(int n,
				int nsub,
				ref int[] ind,
				int nenter,
				int ileave,
				ref int[] indx2,
				int iupdat,
				bool updatd,
				ref double[,] wn,
				ref double[,] wn1,
				int m,
				ref double[,] ws,
				ref double[,] wy,
				ref double[,] sy,
				double theta,
				int col,
				int head,
				ref int info,
				ref double[] workvec,
				ref double[,] workmat)
			{
				int m2 = 0;
				int ipntr = 0;
				int jpntr = 0;
				int iy = 0;
				int iis = 0;
				int jy = 0;
				int js = 0;
				int is1 = 0;
				int js1 = 0;
				int k1 = 0;
				int i = 0;
				int k = 0;
				int col2 = 0;
				int pbegin = 0;
				int pend = 0;
				int dbegin = 0;
				int dend = 0;
				int upcl = 0;
				double temp1 = 0;
				double temp2 = 0;
				double temp3 = 0;
				double temp4 = 0;
				double v = 0;
				int j = 0;
				int i_ = 0;
				int i1_ = 0;

				if (updatd)
				{
					if (iupdat > m)
					{
						for (jy = 1; jy <= m - 1; jy++)
						{
							js = m + jy;
							i1_ = (jy + 1) - (jy);
							for (i_ = jy; i_ <= m - 1; i_++)
							{
								wn1[i_, jy] = wn1[i_ + i1_, jy + 1];
							}
							i1_ = (js + 1) - (js);
							for (i_ = js; i_ <= js + m - jy - 1; i_++)
							{
								wn1[i_, js] = wn1[i_ + i1_, js + 1];
							}
							i1_ = (m + 2) - (m + 1);
							for (i_ = m + 1; i_ <= m + m - 1; i_++)
							{
								wn1[i_, jy] = wn1[i_ + i1_, jy + 1];
							}
						}
					}
					pbegin = 1;
					pend = nsub;
					dbegin = nsub + 1;
					dend = n;
					iy = col;
					iis = m + col;
					ipntr = head + col - 1;
					if (ipntr > m)
					{
						ipntr = ipntr - m;
					}
					jpntr = head;
					for (jy = 1; jy <= col; jy++)
					{
						js = m + jy;
						temp1 = 0;
						temp2 = 0;
						temp3 = 0;
						for (k = pbegin; k <= pend; k++)
						{
							k1 = ind[k];
							temp1 = temp1 + wy[k1, ipntr] * wy[k1, jpntr];
						}
						for (k = dbegin; k <= dend; k++)
						{
							k1 = ind[k];
							temp2 = temp2 + ws[k1, ipntr] * ws[k1, jpntr];
							temp3 = temp3 + ws[k1, ipntr] * wy[k1, jpntr];
						}
						wn1[iy, jy] = temp1;
						wn1[iis, js] = temp2;
						wn1[iis, jy] = temp3;
						jpntr = jpntr % m + 1;
					}
					jy = col;
					jpntr = head + col - 1;
					if (jpntr > m)
					{
						jpntr = jpntr - m;
					}
					ipntr = head;
					for (i = 1; i <= col; i++)
					{
						iis = m + i;
						temp3 = 0;
						for (k = pbegin; k <= pend; k++)
						{
							k1 = ind[k];
							temp3 = temp3 + ws[k1, ipntr] * wy[k1, jpntr];
						}
						ipntr = ipntr % m + 1;
						wn1[iis, jy] = temp3;
					}
					upcl = col - 1;
				}
				else
				{
					upcl = col;
				}
				ipntr = head;
				for (iy = 1; iy <= upcl; iy++)
				{
					iis = m + iy;
					jpntr = head;
					for (jy = 1; jy <= iy; jy++)
					{
						js = m + jy;
						temp1 = 0;
						temp2 = 0;
						temp3 = 0;
						temp4 = 0;
						for (k = 1; k <= nenter; k++)
						{
							k1 = indx2[k];
							temp1 = temp1 + wy[k1, ipntr] * wy[k1, jpntr];
							temp2 = temp2 + ws[k1, ipntr] * ws[k1, jpntr];
						}
						for (k = ileave; k <= n; k++)
						{
							k1 = indx2[k];
							temp3 = temp3 + wy[k1, ipntr] * wy[k1, jpntr];
							temp4 = temp4 + ws[k1, ipntr] * ws[k1, jpntr];
						}
						wn1[iy, jy] = wn1[iy, jy] + temp1 - temp3;
						wn1[iis, js] = wn1[iis, js] - temp2 + temp4;
						jpntr = jpntr % m + 1;
					}
					ipntr = ipntr % m + 1;
				}
				ipntr = head;
				for (iis = m + 1; iis <= m + upcl; iis++)
				{
					jpntr = head;
					for (jy = 1; jy <= upcl; jy++)
					{
						temp1 = 0;
						temp3 = 0;
						for (k = 1; k <= nenter; k++)
						{
							k1 = indx2[k];
							temp1 = temp1 + ws[k1, ipntr] * wy[k1, jpntr];
						}
						for (k = ileave; k <= n; k++)
						{
							k1 = indx2[k];
							temp3 = temp3 + ws[k1, ipntr] * wy[k1, jpntr];
						}
						if (iis <= jy + m)
						{
							wn1[iis, jy] = wn1[iis, jy] + temp1 - temp3;
						}
						else
						{
							wn1[iis, jy] = wn1[iis, jy] - temp1 + temp3;
						}
						jpntr = jpntr % m + 1;
					}
					ipntr = ipntr % m + 1;
				}
				m2 = 2 * m;
				for (iy = 1; iy <= col; iy++)
				{
					iis = col + iy;
					is1 = m + iy;
					for (jy = 1; jy <= iy; jy++)
					{
						js = col + jy;
						js1 = m + jy;
						wn[jy, iy] = wn1[iy, jy] / theta;
						wn[js, iis] = wn1[is1, js1] * theta;
					}
					for (jy = 1; jy <= iy - 1; jy++)
					{
						wn[jy, iis] = -wn1[is1, jy];
					}
					for (jy = iy; jy <= col; jy++)
					{
						wn[jy, iis] = wn1[is1, jy];
					}
					wn[iy, iy] = wn[iy, iy] + sy[iy, iy];
				}
				info = 0;
				if (!lbfgsbdpofa(ref wn, col))
				{
					info = -1;
					return;
				}
				col2 = 2 * col;
				for (js = col + 1; js <= col2; js++)
				{
					for (i_ = 1; i_ <= col; i_++)
					{
						workvec[i_] = wn[i_, js];
					}
					lbfgsbdtrsl(ref wn, col, ref workvec, 11, ref info);
					for (i_ = 1; i_ <= col; i_++)
					{
						wn[i_, js] = workvec[i_];
					}
				}
				for (iis = col + 1; iis <= col2; iis++)
				{
					for (js = iis; js <= col2; js++)
					{
						v = 0.0;
						for (i_ = 1; i_ <= col; i_++)
						{
							v += wn[i_, iis] * wn[i_, js];
						}
						wn[iis, js] = wn[iis, js] + v;
					}
				}
				for (j = 1; j <= col; j++)
				{
					i1_ = (col + 1) - (1);
					for (i_ = 1; i_ <= col; i_++)
					{
						workmat[j, i_] = wn[col + j, i_ + i1_];
					}
				}
				info = 0;
				if (!lbfgsbdpofa(ref workmat, col))
				{
					info = -2;
					return;
				}
				for (j = 1; j <= col; j++)
				{
					i1_ = (1) - (col + 1);
					for (i_ = col + 1; i_ <= col + col; i_++)
					{
						wn[col + j, i_] = workmat[j, i_ + i1_];
					}
				}
			}

			private void lbfgsbformt(int m,
				ref double[,] wt,
				ref double[,] sy,
				ref double[,] ss,
				int col,
				double theta,
				ref int info)
			{
				int i = 0;
				int j = 0;
				int k = 0;
				int k1 = 0;
				double ddum = 0;

				for (j = 1; j <= col; j++)
				{
					wt[1, j] = theta * ss[1, j];
				}
				for (i = 2; i <= col; i++)
				{
					for (j = i; j <= col; j++)
					{
						k1 = Math.Min(i, j) - 1;
						ddum = 0;
						for (k = 1; k <= k1; k++)
						{
							ddum = ddum + sy[i, k] * sy[j, k] / sy[k, k];
						}
						wt[i, j] = ddum + theta * ss[i, j];
					}
				}
				info = 0;
				if (!lbfgsbdpofa(ref wt, col))
				{
					info = -3;
				}
			}

			private void lbfgsbfreev(int n,
				ref int nfree,
				ref int[] index,
				ref int nenter,
				ref int ileave,
				ref int[] indx2,
				ref int[] iwhere,
				ref bool wrk,
				bool updatd,
				bool cnstnd,
				int iter)
			{
				int iact = 0;
				int i = 0;
				int k = 0;

				nenter = 0;
				ileave = n + 1;
				if (iter > 0 & cnstnd)
				{
					for (i = 1; i <= nfree; i++)
					{
						k = index[i];
						if (iwhere[k] > 0)
						{
							ileave = ileave - 1;
							indx2[ileave] = k;
						}
					}
					for (i = 1 + nfree; i <= n; i++)
					{
						k = index[i];
						if (iwhere[k] <= 0)
						{
							nenter = nenter + 1;
							indx2[nenter] = k;
						}
					}
				}
				wrk = ileave < n + 1 | nenter > 0 | updatd;
				nfree = 0;
				iact = n + 1;
				for (i = 1; i <= n; i++)
				{
					if (iwhere[i] <= 0)
					{
						nfree = nfree + 1;
						index[nfree] = i;
					}
					else
					{
						iact = iact - 1;
						index[iact] = i;
					}
				}
			}

			private void lbfgsbhpsolb(int n,
				ref double[] t,
				ref int[] iorder,
				int iheap)
			{
				int i = 0;
				int j = 0;
				int k = 0;
				int indxin = 0;
				int indxou = 0;
				double ddum = 0;
				double dout = 0;

				if (iheap == 0)
				{
					for (k = 2; k <= n; k++)
					{
						ddum = t[k];
						indxin = iorder[k];
						i = k;
						while (true)
						{
							if (i > 1)
							{
								j = i / 2;
								if (ddum < t[j])
								{
									t[i] = t[j];
									iorder[i] = iorder[j];
									i = j;
									continue;
								}
							}
							break;
						}
						t[i] = ddum;
						iorder[i] = indxin;
					}
				}
				if (n > 1)
				{
					i = 1;
					dout = t[1];
					indxou = iorder[1];
					ddum = t[n];
					indxin = iorder[n];
					while (true)
					{
						j = i + i;
						if (j <= n - 1)
						{
							if (t[j + 1] < t[j])
							{
								j = j + 1;
							}
							if (t[j] < ddum)
							{
								t[i] = t[j];
								iorder[i] = iorder[j];
								i = j;
								continue;
							}
						}
						break;
					}
					t[i] = ddum;
					iorder[i] = indxin;
					t[n] = dout;
					iorder[n] = indxou;
				}
			}

			private void lbfgsblnsrlb(int n,
				ref double[] l,
				ref double[] u,
				ref int[] nbd,
				ref double[] x,
				double f,
				ref double fold,
				ref double gd,
				ref double gdold,
				ref double[] g,
				ref double[] d,
				ref double[] r,
				ref double[] t,
				ref double[] z,
				ref double stp,
				ref double dnrm,
				ref double dtd,
				ref double xstep,
				ref double stpmx,
				int iter,
				ref int ifun,
				ref int iback,
				ref int nfgv,
				ref int info,
				ref int task,
				bool boxed,
				bool cnstnd,
				ref int csave,
				ref int[] isave,
				ref double[] dsave)
			{
				int i = 0;
				double a1 = 0;
				double a2 = 0;
				double v = 0;
				double ftol = 0;
				double gtol = 0;
				double xtol = 0;
				double big = 0;
				int addinfo = 0;
				int i_ = 0;

				addinfo = 0;
				big = 1.0E10;
				ftol = 1.0E-3;
				gtol = 0.9E0;
				xtol = 0.1E0;
				if (task != 1)
				{
					v = 0.0;
					for (i_ = 1; i_ <= n; i_++)
					{
						v += d[i_] * d[i_];
					}
					dtd = v;
					dnrm = Math.Sqrt(dtd);
					stpmx = big;
					if (cnstnd)
					{
						if (iter == 0)
						{
							stpmx = 1;
						}
						else
						{
							for (i = 1; i <= n; i++)
							{
								a1 = d[i];
								if (nbd[i] != 0)
								{
									if (a1 < 0 & nbd[i] <= 2)
									{
										a2 = l[i] - x[i];
										if (a2 >= 0)
										{
											stpmx = 0;
										}
										else
										{
											if (a1 * stpmx < a2)
											{
												stpmx = a2 / a1;
											}
										}
									}
									else
									{
										if (a1 > 0 & nbd[i] >= 2)
										{
											a2 = u[i] - x[i];
											if (a2 <= 0)
											{
												stpmx = 0;
											}
											else
											{
												if (a1 * stpmx > a2)
												{
													stpmx = a2 / a1;
												}
											}
										}
									}
								}
							}
						}
					}
					if (iter == 0 & !boxed)
					{
						stp = Math.Min(1 / dnrm, stpmx);
					}
					else
					{
						stp = 1;
					}
					for (i_ = 1; i_ <= n; i_++)
					{
						t[i_] = x[i_];
					}
					for (i_ = 1; i_ <= n; i_++)
					{
						r[i_] = g[i_];
					}
					fold = f;
					ifun = 0;
					iback = 0;
					csave = 0;
				}
				v = 0.0;
				for (i_ = 1; i_ <= n; i_++)
				{
					v += g[i_] * d[i_];
				}
				gd = v;
				if (ifun == 0)
				{
					gdold = gd;
					if (gd >= 0)
					{
						info = -4;
						return;
					}
				}
				lbfgsbdcsrch(f, gd, ref stp, ftol, gtol, xtol, 0, stpmx, ref csave, ref isave, ref dsave, ref addinfo);
				xstep = stp * dnrm;
				if (csave != 4 & csave != 3)
				{
					task = 1;
					ifun = ifun + 1;
					nfgv = nfgv + 1;
					iback = ifun - 1;
					if (stp == 1)
					{
						for (i_ = 1; i_ <= n; i_++)
						{
							x[i_] = z[i_];
						}
					}
					else
					{
						for (i = 1; i <= n; i++)
						{
							x[i] = stp * d[i] + t[i];
						}
					}
				}
				else
				{
					task = 5;
				}
			}

			private void lbfgsbmatupd(int n,
				int m,
				ref double[,] ws,
				ref double[,] wy,
				ref double[,] sy,
				ref double[,] ss,
				ref double[] d,
				ref double[] r,
				ref int itail,
				int iupdat,
				ref int col,
				ref int head,
				ref double theta,
				double rr,
				double dr,
				double stp,
				double dtd)
			{
				int j = 0;
				int pointr = 0;
				double v = 0;
				int i_ = 0;
				int i1_ = 0;

				if (iupdat <= m)
				{
					col = iupdat;
					itail = (head + iupdat - 2) % m + 1;
				}
				else
				{
					itail = itail % m + 1;
					head = head % m + 1;
				}
				for (i_ = 1; i_ <= n; i_++)
				{
					ws[i_, itail] = d[i_];
				}
				for (i_ = 1; i_ <= n; i_++)
				{
					wy[i_, itail] = r[i_];
				}
				theta = rr / dr;
				if (iupdat > m)
				{
					for (j = 1; j <= col - 1; j++)
					{
						i1_ = (2) - (1);
						for (i_ = 1; i_ <= j; i_++)
						{
							ss[i_, j] = ss[i_ + i1_, j + 1];
						}
						i1_ = (j + 1) - (j);
						for (i_ = j; i_ <= col - 1; i_++)
						{
							sy[i_, j] = sy[i_ + i1_, j + 1];
						}
					}
				}
				pointr = head;
				for (j = 1; j <= col - 1; j++)
				{
					v = 0.0;
					for (i_ = 1; i_ <= n; i_++)
					{
						v += d[i_] * wy[i_, pointr];
					}
					sy[col, j] = v;
					v = 0.0;
					for (i_ = 1; i_ <= n; i_++)
					{
						v += ws[i_, pointr] * d[i_];
					}
					ss[j, col] = v;
					pointr = pointr % m + 1;
				}
				if (stp == 1)
				{
					ss[col, col] = dtd;
				}
				else
				{
					ss[col, col] = stp * stp * dtd;
				}
				sy[col, col] = dr;
			}

			private void lbfgsbprojgr(int n,
				ref double[] l,
				ref double[] u,
				ref int[] nbd,
				ref double[] x,
				ref double[] g,
				ref double sbgnrm)
			{
				int i = 0;
				double gi = 0;

				sbgnrm = 0;
				for (i = 1; i <= n; i++)
				{
					gi = g[i];
					if (nbd[i] != 0)
					{
						if (gi < 0)
						{
							if (nbd[i] >= 2)
							{
								gi = Math.Max(x[i] - u[i], gi);
							}
						}
						else
						{
							if (nbd[i] <= 2)
							{
								gi = Math.Min(x[i] - l[i], gi);
							}
						}
					}
					sbgnrm = Math.Max(sbgnrm, Math.Abs(gi));
				}
			}

			private void lbfgsbsubsm(int n,
				int m,
				int nsub,
				ref int[] ind,
				ref double[] l,
				ref double[] u,
				ref int[] nbd,
				ref double[] x,
				ref double[] d,
				ref double[,] ws,
				ref double[,] wy,
				double theta,
				int col,
				int head,
				ref int iword,
				ref double[] wv,
				ref double[,] wn,
				ref int info)
			{
				int pointr = 0;
				int m2 = 0;
				int col2 = 0;
				int ibd = 0;
				int jy = 0;
				int js = 0;
				int i = 0;
				int j = 0;
				int k = 0;
				double alpha = 0;
				double dk = 0;
				double temp1 = 0;
				double temp2 = 0;

				if (nsub <= 0)
				{
					return;
				}
				pointr = head;
				for (i = 1; i <= col; i++)
				{
					temp1 = 0;
					temp2 = 0;
					for (j = 1; j <= nsub; j++)
					{
						k = ind[j];
						temp1 = temp1 + wy[k, pointr] * d[j];
						temp2 = temp2 + ws[k, pointr] * d[j];
					}
					wv[i] = temp1;
					wv[col + i] = theta * temp2;
					pointr = pointr % m + 1;
				}
				m2 = 2 * m;
				col2 = 2 * col;
				lbfgsbdtrsl(ref wn, col2, ref wv, 11, ref info);
				if (info != 0)
				{
					return;
				}
				for (i = 1; i <= col; i++)
				{
					wv[i] = -wv[i];
				}
				lbfgsbdtrsl(ref wn, col2, ref wv, 1, ref info);
				if (info != 0)
				{
					return;
				}
				pointr = head;
				for (jy = 1; jy <= col; jy++)
				{
					js = col + jy;
					for (i = 1; i <= nsub; i++)
					{
						k = ind[i];
						d[i] = d[i] + wy[k, pointr] * wv[jy] / theta + ws[k, pointr] * wv[js];
					}
					pointr = pointr % m + 1;
				}
				for (i = 1; i <= nsub; i++)
				{
					d[i] = d[i] / theta;
				}
				alpha = 1;
				temp1 = alpha;
				for (i = 1; i <= nsub; i++)
				{
					k = ind[i];
					dk = d[i];
					if (nbd[k] != 0)
					{
						if (dk < 0 & nbd[k] <= 2)
						{
							temp2 = l[k] - x[k];
							if (temp2 >= 0)
							{
								temp1 = 0;
							}
							else
							{
								if (dk * alpha < temp2)
								{
									temp1 = temp2 / dk;
								}
							}
						}
						else
						{
							if (dk > 0 & nbd[k] >= 2)
							{
								temp2 = u[k] - x[k];
								if (temp2 <= 0)
								{
									temp1 = 0;
								}
								else
								{
									if (dk * alpha > temp2)
									{
										temp1 = temp2 / dk;
									}
								}
							}
						}
						if (temp1 < alpha)
						{
							alpha = temp1;
							ibd = i;
						}
					}
				}
				if (alpha < 1)
				{
					dk = d[ibd];
					k = ind[ibd];
					if (dk > 0)
					{
						x[k] = u[k];
						d[ibd] = 0;
					}
					else
					{
						if (dk < 0)
						{
							x[k] = l[k];
							d[ibd] = 0;
						}
					}
				}
				for (i = 1; i <= nsub; i++)
				{
					k = ind[i];
					x[k] = x[k] + alpha * d[i];
				}
				if (alpha < 1)
				{
					iword = 1;
				}
				else
				{
					iword = 0;
				}
			}

			private void lbfgsbdcsrch(double f,
				double g,
				ref double stp,
				double ftol,
				double gtol,
				double xtol,
				double stpmin,
				double stpmax,
				ref int task,
				ref int[] isave,
				ref double[] dsave,
				ref int addinfo)
			{
				bool brackt = new bool();
				int stage = 0;
				double finit = 0;
				double ftest = 0;
				double fm = 0;
				double fx = 0;
				double fxm = 0;
				double fy = 0;
				double fym = 0;
				double ginit = 0;
				double gtest = 0;
				double gm = 0;
				double gx = 0;
				double gxm = 0;
				double gy = 0;
				double gym = 0;
				double stx = 0;
				double sty = 0;
				double stmin = 0;
				double stmax = 0;
				double width = 0;
				double width1 = 0;
				double xtrapl = 0;
				double xtrapu = 0;

				xtrapl = 1.1E0;
				xtrapu = 4.0E0;
				while (true)
				{
					if (task == 0)
					{
						if (stp < stpmin)
						{
							task = 2;
							addinfo = 0;
						}
						if (stp > stpmax)
						{
							task = 2;
							addinfo = 0;
						}
						if (g >= 0)
						{
							task = 2;
							addinfo = 0;
						}
						if (ftol < 0)
						{
							task = 2;
							addinfo = 0;
						}
						if (gtol < 0)
						{
							task = 2;
							addinfo = 0;
						}
						if (xtol < 0)
						{
							task = 2;
							addinfo = 0;
						}
						if (stpmin < 0)
						{
							task = 2;
							addinfo = 0;
						}
						if (stpmax < stpmin)
						{
							task = 2;
							addinfo = 0;
						}
						if (task == 2)
						{
							return;
						}
						brackt = false;
						stage = 1;
						finit = f;
						ginit = g;
						gtest = ftol * ginit;
						width = stpmax - stpmin;
						width1 = width / 0.5;
						stx = 0;
						fx = finit;
						gx = ginit;
						sty = 0;
						fy = finit;
						gy = ginit;
						stmin = 0;
						stmax = stp + xtrapu * stp;
						task = 1;
						break;
					}
					else
					{
						if (isave[1] == 1)
						{
							brackt = true;
						}
						else
						{
							brackt = false;
						}
						stage = isave[2];
						ginit = dsave[1];
						gtest = dsave[2];
						gx = dsave[3];
						gy = dsave[4];
						finit = dsave[5];
						fx = dsave[6];
						fy = dsave[7];
						stx = dsave[8];
						sty = dsave[9];
						stmin = dsave[10];
						stmax = dsave[11];
						width = dsave[12];
						width1 = dsave[13];
					}
					ftest = finit + stp * gtest;
					if (stage == 1 & f <= ftest & g >= 0)
					{
						stage = 2;
					}
					if (brackt & (stp <= stmin | stp >= stmax))
					{
						task = 3;
						addinfo = 1;
					}
					if (brackt & stmax - stmin <= xtol * stmax)
					{
						task = 3;
						addinfo = 2;
					}
					if (stp == stpmax & f <= ftest & g <= gtest)
					{
						task = 3;
						addinfo = 3;
					}
					if (stp == stpmin & (f > ftest | g >= gtest))
					{
						task = 3;
						addinfo = 4;
					}
					if (f <= ftest & Math.Abs(g) <= gtol * -ginit)
					{
						task = 4;
						addinfo = -1;
					}
					if (task == 3 | task == 4)
					{
						break;
					}
					if (stage == 1 & f <= fx & f > ftest)
					{
						fm = f - stp * gtest;
						fxm = fx - stx * gtest;
						fym = fy - sty * gtest;
						gm = g - gtest;
						gxm = gx - gtest;
						gym = gy - gtest;
						lbfgsbdcstep(ref stx, ref fxm, ref gxm, ref sty, ref fym, ref gym, ref stp, fm, gm, ref brackt, stmin, stmax);
						fx = fxm + stx * gtest;
						fy = fym + sty * gtest;
						gx = gxm + gtest;
						gy = gym + gtest;
					}
					else
					{
						lbfgsbdcstep(ref stx, ref fx, ref gx, ref sty, ref fy, ref gy, ref stp, f, g, ref brackt, stmin, stmax);
					}
					if (brackt)
					{
						if (Math.Abs(sty - stx) >= 0.66 * width1)
						{
							stp = stx + 0.5 * (sty - stx);
						}
						width1 = width;
						width = Math.Abs(sty - stx);
					}
					if (brackt)
					{
						stmin = Math.Min(stx, sty);
						stmax = Math.Max(stx, sty);
					}
					else
					{
						stmin = stp + xtrapl * (stp - stx);
						stmax = stp + xtrapu * (stp - stx);
					}
					stp = Math.Max(stp, stpmin);
					stp = Math.Min(stp, stpmax);
					if (brackt & (stp <= stmin | stp >= stmax) | brackt & stmax - stmin <= xtol * stmax)
					{
						stp = stx;
					}
					task = 1;
					break;
				}
				if (brackt)
				{
					isave[1] = 1;
				}
				else
				{
					isave[1] = 0;
				}
				isave[2] = stage;
				dsave[1] = ginit;
				dsave[2] = gtest;
				dsave[3] = gx;
				dsave[4] = gy;
				dsave[5] = finit;
				dsave[6] = fx;
				dsave[7] = fy;
				dsave[8] = stx;
				dsave[9] = sty;
				dsave[10] = stmin;
				dsave[11] = stmax;
				dsave[12] = width;
				dsave[13] = width1;
			}

			private void lbfgsbdcstep(ref double stx,
				ref double fx,
				ref double dx,
				ref double sty,
				ref double fy,
				ref double dy,
				ref double stp,
				double fp,
				double dp,
				ref bool brackt,
				double stpmin,
				double stpmax)
			{
				double gamma = 0;
				double p = 0;
				double q = 0;
				double r = 0;
				double s = 0;
				double sgnd = 0;
				double stpc = 0;
				double stpf = 0;
				double stpq = 0;
				double theta = 0;

				sgnd = dp * (dx / Math.Abs(dx));
				if (fp > fx)
				{
					theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
					s = Math.Max(Math.Abs(theta), Math.Max(Math.Abs(dx), Math.Abs(dp)));
					gamma = s * Math.Sqrt(Sqr(theta / s) - dx / s * (dp / s));
					if (stp < stx)
					{
						gamma = -gamma;
					}
					p = gamma - dx + theta;
					q = gamma - dx + gamma + dp;
					r = p / q;
					stpc = stx + r * (stp - stx);
					stpq = stx + dx / ((fx - fp) / (stp - stx) + dx) / 2 * (stp - stx);
					if (Math.Abs(stpc - stx) < Math.Abs(stpq - stx))
					{
						stpf = stpc;
					}
					else
					{
						stpf = stpc + (stpq - stpc) / 2;
					}
					brackt = true;
				}
				else
				{
					if (sgnd < 0)
					{
						theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
						s = Math.Max(Math.Abs(theta), Math.Max(Math.Abs(dx), Math.Abs(dp)));
						gamma = s * Math.Sqrt(Sqr(theta / s) - dx / s * (dp / s));
						if (stp > stx)
						{
							gamma = -gamma;
						}
						p = gamma - dp + theta;
						q = gamma - dp + gamma + dx;
						r = p / q;
						stpc = stp + r * (stx - stp);
						stpq = stp + dp / (dp - dx) * (stx - stp);
						if (Math.Abs(stpc - stp) > Math.Abs(stpq - stp))
						{
							stpf = stpc;
						}
						else
						{
							stpf = stpq;
						}
						brackt = true;
					}
					else
					{
						if (Math.Abs(dp) < Math.Abs(dx))
						{
							theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
							s = Math.Max(Math.Abs(theta), Math.Max(Math.Abs(dx), Math.Abs(dp)));
							gamma = s * Math.Sqrt(Math.Max(0, Sqr(theta / s) - dx / s * (dp / s)));
							if (stp > stx)
							{
								gamma = -gamma;
							}
							p = gamma - dp + theta;
							q = gamma + (dx - dp) + gamma;
							r = p / q;
							if (r < 0 & gamma != 0)
							{
								stpc = stp + r * (stx - stp);
							}
							else
							{
								if (stp > stx)
								{
									stpc = stpmax;
								}
								else
								{
									stpc = stpmin;
								}
							}
							stpq = stp + dp / (dp - dx) * (stx - stp);
							if (brackt)
							{
								if (Math.Abs(stpc - stp) < Math.Abs(stpq - stp))
								{
									stpf = stpc;
								}
								else
								{
									stpf = stpq;
								}
								if (stp > stx)
								{
									stpf = Math.Min(stp + 0.66 * (sty - stp), stpf);
								}
								else
								{
									stpf = Math.Max(stp + 0.66 * (sty - stp), stpf);
								}
							}
							else
							{
								if (Math.Abs(stpc - stp) > Math.Abs(stpq - stp))
								{
									stpf = stpc;
								}
								else
								{
									stpf = stpq;
								}
								stpf = Math.Min(stpmax, stpf);
								stpf = Math.Max(stpmin, stpf);
							}
						}
						else
						{
							if (brackt)
							{
								theta = 3 * (fp - fy) / (sty - stp) + dy + dp;
								s = Math.Max(Math.Abs(theta), Math.Max(Math.Abs(dy), Math.Abs(dp)));
								gamma = s * Math.Sqrt(Sqr(theta / s) - dy / s * (dp / s));
								if (stp > sty)
								{
									gamma = -gamma;
								}
								p = gamma - dp + theta;
								q = gamma - dp + gamma + dy;
								r = p / q;
								stpc = stp + r * (sty - stp);
								stpf = stpc;
							}
							else
							{
								if (stp > stx)
								{
									stpf = stpmax;
								}
								else
								{
									stpf = stpmin;
								}
							}
						}
					}
				}
				if (fp > fx)
				{
					sty = stp;
					fy = fp;
					dy = dp;
				}
				else
				{
					if (sgnd < 0)
					{
						sty = stx;
						fy = fx;
						dy = dx;
					}
					stx = stp;
					fx = fp;
					dx = dp;
				}
				stp = stpf;
			}

			private bool additionallbfgsbstoppingcriterion(int iter,
				ref double[] x,
				double f,
				ref double[] g)
			{
				bool result = new bool();

				result = false;
				return result;
			}

			private bool lbfgsbdpofa(ref double[,] a,
				int n)
			{
				bool result = new bool();
				double t = 0;
				double s = 0;
				double v = 0;
				int j = 0;
				int jm1 = 0;
				int k = 0;
				int i_ = 0;

				for (j = 1; j <= n; j++)
				{
					s = 0.0;
					jm1 = j - 1;
					if (jm1 >= 1)
					{
						for (k = 1; k <= jm1; k++)
						{
							v = 0.0;
							for (i_ = 1; i_ <= k - 1; i_++)
							{
								v += a[i_, k] * a[i_, j];
							}
							t = a[k, j] - v;
							t = t / a[k, k];
							a[k, j] = t;
							s = s + t * t;
						}
					}
					s = a[j, j] - s;
					if (s <= 0.0)
					{
						result = false;
						return result;
					}
					a[j, j] = Math.Sqrt(s);
				}
				result = true;
				return result;
			}

			private void lbfgsbdtrsl(ref double[,] t,
				int n,
				ref double[] b,
				int job,
				ref int info)
			{
				double temp = 0;
				double v = 0;
				int cse = 0;
				int j = 0;
				int jj = 0;
				int i_ = 0;

				for (j = 1; j <= n; j++)
				{
					if (t[j, j] == 0.0)
					{
						info = j;
						return;
					}
				}
				info = 0;
				cse = 1;
				if (job % 10 != 0)
				{
					cse = 2;
				}
				if (job % 100 / 10 != 0)
				{
					cse = cse + 2;
				}
				if (cse == 1)
				{
					b[1] = b[1] / t[1, 1];
					if (n < 2)
					{
						return;
					}
					for (j = 2; j <= n; j++)
					{
						temp = -b[j - 1];
						for (i_ = j; i_ <= n; i_++)
						{
							b[i_] = b[i_] + temp * t[i_, j - 1];
						}
						b[j] = b[j] / t[j, j];
					}
					return;
				}
				if (cse == 2)
				{
					b[n] = b[n] / t[n, n];
					if (n < 2)
					{
						return;
					}
					for (jj = 2; jj <= n; jj++)
					{
						j = n - jj + 1;
						temp = -b[j + 1];
						for (i_ = 1; i_ <= j; i_++)
						{
							b[i_] = b[i_] + temp * t[i_, j + 1];
						}
						b[j] = b[j] / t[j, j];
					}
					return;
				}
				if (cse == 3)
				{
					b[n] = b[n] / t[n, n];
					if (n < 2)
					{
						return;
					}
					for (jj = 2; jj <= n; jj++)
					{
						j = n - jj + 1;
						v = 0.0;
						for (i_ = j + 1; i_ <= j + 1 + jj - 1 - 1; i_++)
						{
							v += t[i_, j] * b[i_];
						}
						b[j] = b[j] - v;
						b[j] = b[j] / t[j, j];
					}
					return;
				}
				if (cse == 4)
				{
					b[1] = b[1] / t[1, 1];
					if (n < 2)
					{
						return;
					}
					for (j = 2; j <= n; j++)
					{
						v = 0.0;
						for (i_ = 1; i_ <= j - 1; i_++)
						{
							v += t[i_, j] * b[i_];
						}
						b[j] = b[j] - v;
						b[j] = b[j] / t[j, j];
					}
					return;
				}
			}

			private void lbfgsbnewiteration(ref double[] x,
				double f,
				ref double[] g)
			{
			}

			/*************************************************************************
				LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
								  JORGE NOCEDAL

			The subroutine minimizes function F(x) of N arguments by  using  a  quasi-
			Newton method (LBFGS scheme) which is optimized to use  a  minimum  amount
			of memory.

			The subroutine generates the approximation of an inverse Hessian matrix by
			using information about the last M steps of the algorithm  (instead of N).
			It lessens a required amount of memory from a value  of  order  N^2  to  a
			value of order 2*N*M.

			This subroutine uses the FuncGrad subroutine which calculates the value of
			the function F and gradient G in point X. The programmer should define the
			FuncGrad subroutine by himself.  It  should  be  noted that the subroutine
			doesn't need to waste time for memory allocation of array G,  because  the
			memory is allocated in calling the subroutine. Setting a dimension of array
			G  each  time  when  calling  a  subroutine  will excessively slow down an
			algorithm.

			The programmer could also redefine the LBFGSNewIteration subroutine  which
			is called on each new step. The current point X, the function value F  and
			the  gradient  G  are  passed  into  this  subroutine. It is reasonable to
			redefine the subroutine for better debugging, for  example,  to  visualize
			the solution process.

			Input parameters:
				N   -   problem dimension. N>0
				M   -   number of corrections in the BFGS scheme of Hessian
						approximation update. Recommended value:  3<=M<=7. The smaller
						value causes worse convergence, the bigger will  not  cause  a
						considerably better convergence, but will cause a fall in  the
						performance. M<=N.
				X   -   initial solution approximation.
						Array whose index ranges from 1 to N.
				EpsG -  positive number which  defines  a  precision  of  search.  The
						subroutine finishes its work if the condition ||G|| < EpsG  is
						satisfied, where ||.|| means Euclidian norm, G - gradient, X -
						current approximation.
				EpsF -  positive number which  defines  a  precision  of  search.  The
						subroutine finishes its work if on iteration  number  k+1  the
						condition |F(k+1)-F(k)| <= EpsF*max{|F(k)|, |F(k+1)|, 1}    is
						satisfied.
				EpsX -  positive number which  defines  a  precision  of  search.  The
						subroutine finishes its work if on iteration number k+1    the
						condition |X(k+1)-X(k)| <= EpsX is fulfilled.
				MaxIts- maximum number of iterations. If MaxIts=0, the number of
						iterations is unlimited.

			Output parameters:
				X   -   solution approximation. Array whose index ranges from 1 to N.
				Info-   a return code:
								* -1 wrong parameters were specified,
								* 0 interrupted by user,
								* 1 relative function decreasing is less or equal to EpsF,
								* 2 step is less or equal EpsX,
								* 4 gradient norm is less or equal to EpsG,
								* 5 number of iterations exceeds MaxIts.

			FuncGrad routine description. User-defined.
			Input parameters:
				X   -   array whose index ranges from 1 to N.
			Output parameters:
				F   -   function value at X.
				G   -   function gradient.
						Array whose index ranges from 1 to N.
			The memory for array G has already been allocated in the calling subroutine,
			and it isn't necessary to allocate it in the FuncGrad subroutine.
			*************************************************************************/

			private void lbfgsminimize(int n,
				int m,
				ref double[] x,
				double epsg,
				double epsf,
				double epsx,
				int maxits,
			ref int info)
			{
				double[] w = new double[0];
				double f = 0;
				double fold = 0;
				double tf = 0;
				double txnorm = 0;
				double v = 0;
				double[] xold = new double[0];
				double[] tx = new double[0];
				double[] g = new double[0];
				double[] diag = new double[0];
				double[] ta = new double[0];
				//bool finish = new bool();
				double gnorm = 0;
				double stp1 = 0;
				double ftol = 0;
				double stp = 0;
				double ys = 0;
				double yy = 0;
				double sq = 0;
				double yr = 0;
				double beta = 0;
				double xnorm = 0;
				int iter = 0;
				int nfun = 0;
				int point = 0;
				int ispt = 0;
				int iypt = 0;
				int maxfev = 0;
				int bound = 0;
				int npt = 0;
				int cp = 0;
				int i = 0;
				int nfev = 0;
				int inmc = 0;
				int iycn = 0;
				int iscn = 0;
				double xtol = 0;
				double gtol = 0;
				double stpmin = 0;
				double stpmax = 0;
				int i_ = 0;

				w = new double[n * (2 * m + 1) + 2 * m + 1];
				g = new double[n + 1];
				xold = new double[n + 1];
				tx = new double[n + 1];
				diag = new double[n + 1];
				ta = new double[n + 1];
				funcgrad(ref x, ref f, ref g);
				fold = f;
				iter = 0;
				info = 0;
				if (n <= 0 | m <= 0 | m > n | epsg < 0 | epsf < 0 | epsx < 0 | maxits < 0)
				{
					info = -1;
					return;
				}
				nfun = 1;
				point = 0;
				//finish = false;
				for (i = 1; i <= n; i++)
				{
					diag[i] = 1;
				}
				xtol = 100 * MachineEpsilon;
				gtol = 0.9;
				stpmin = Math.Pow(10, -20);
				stpmax = Math.Pow(10, 20);
				ispt = n + 2 * m;
				iypt = ispt + n * m;
				for (i = 1; i <= n; i++)
				{
					w[ispt + i] = -(g[i] * diag[i]);
				}
				gnorm = Math.Sqrt(lbfgsdotproduct(n, ref g, 1, ref g, 1));
				stp1 = 1 / gnorm;
				ftol = 0.0001;
				maxfev = 20;
				while (true)
				{
					for (i_ = 1; i_ <= n; i_++)
					{
						xold[i_] = x[i_];
					}
					iter = iter + 1;
					info = 0;
					bound = iter - 1;
					if (iter != 1)
					{
						if (iter > m)
						{
							bound = m;
						}
						ys = lbfgsdotproduct(n, ref w, iypt + npt + 1, ref w, ispt + npt + 1);
						yy = lbfgsdotproduct(n, ref w, iypt + npt + 1, ref w, iypt + npt + 1);
						for (i = 1; i <= n; i++)
						{
							diag[i] = ys / yy;
						}
						cp = point;
						if (point == 0)
						{
							cp = m;
						}
						w[n + cp] = 1 / ys;
						for (i = 1; i <= n; i++)
						{
							w[i] = -g[i];
						}
						cp = point;
						for (i = 1; i <= bound; i++)
						{
							cp = cp - 1;
							if (cp == -1)
							{
								cp = m - 1;
							}
							sq = lbfgsdotproduct(n, ref w, ispt + cp * n + 1, ref w, 1);
							inmc = n + m + cp + 1;
							iycn = iypt + cp * n;
							w[inmc] = w[n + cp + 1] * sq;
							lbfgslincomb(n, -w[inmc], ref w, iycn + 1, ref w, 1);
						}
						for (i = 1; i <= n; i++)
						{
							w[i] = diag[i] * w[i];
						}
						for (i = 1; i <= bound; i++)
						{
							yr = lbfgsdotproduct(n, ref w, iypt + cp * n + 1, ref w, 1);
							beta = w[n + cp + 1] * yr;
							inmc = n + m + cp + 1;
							beta = w[inmc] - beta;
							iscn = ispt + cp * n;
							lbfgslincomb(n, beta, ref w, iscn + 1, ref w, 1);
							cp = cp + 1;
							if (cp == m)
							{
								cp = 0;
							}
						}
						for (i = 1; i <= n; i++)
						{
							w[ispt + point * n + i] = w[i];
						}
					}
					nfev = 0;
					stp = 1;
					if (iter == 1)
					{
						stp = stp1;
					}
					for (i = 1; i <= n; i++)
					{
						w[i] = g[i];
					}
					lbfgsmcsrch(n, ref x, ref f, ref g, ref w, ispt + point * n + 1, ref stp, ftol, xtol, maxfev, ref info, ref nfev, ref diag, gtol, stpmin, stpmax);
					if (info != 1)
					{
						if (info == 0)
						{
							info = -1;
							return;
						}
					}
					nfun = nfun + nfev;
					npt = point * n;
					for (i = 1; i <= n; i++)
					{
						w[ispt + npt + i] = stp * w[ispt + npt + i];
						w[iypt + npt + i] = g[i] - w[i];
					}
					point = point + 1;
					if (point == m)
					{
						point = 0;
					}
					if (iter > maxits & maxits > 0)
					{
						info = 5;
						return;
					}
					lbfgsnewiteration(ref x, f, ref g);
					gnorm = Math.Sqrt(lbfgsdotproduct(n, ref g, 1, ref g, 1));
					if (gnorm <= epsg)
					{
						info = 4;
						return;
					}
					tf = Math.Max(Math.Abs(fold), Math.Max(Math.Abs(f), 1.0));
					if (fold - f <= epsf * tf)
					{
						info = 1;
						return;
					}
					for (i_ = 1; i_ <= n; i_++)
					{
						tx[i_] = xold[i_];
					}
					for (i_ = 1; i_ <= n; i_++)
					{
						tx[i_] = tx[i_] - x[i_];
					}
					xnorm = Math.Sqrt(lbfgsdotproduct(n, ref x, 1, ref x, 1));
					txnorm = Math.Max(xnorm, Math.Sqrt(lbfgsdotproduct(n, ref xold, 1, ref xold, 1)));
					txnorm = Math.Max(txnorm, 1.0);
					v = Math.Sqrt(lbfgsdotproduct(n, ref tx, 1, ref tx, 1));
					if (v <= epsx)
					{
						info = 2;
						return;
					}
					fold = f;
					for (i_ = 1; i_ <= n; i_++)
					{
						xold[i_] = x[i_];
					}
				}
			}

			private void lbfgslincomb(int n,
				double da,
				ref double[] dx,
				int sx,
				ref double[] dy,
				int sy)
			{
				int fx = 0;
				int fy = 0;
				int i_ = 0;
				int i1_ = 0;

				fx = sx + n - 1;
				fy = sy + n - 1;
				i1_ = (sx) - (sy);
				for (i_ = sy; i_ <= fy; i_++)
				{
					dy[i_] = dy[i_] + da * dx[i_ + i1_];
				}
			}

			private double lbfgsdotproduct(int n,
				ref double[] dx,
				int sx,
				ref double[] dy,
				int sy)
			{
				double result = 0;
				double v = 0;
				int fx = 0;
				int fy = 0;
				int i_ = 0;
				int i1_ = 0;

				fx = sx + n - 1;
				fy = sy + n - 1;
				i1_ = (sy) - (sx);
				v = 0.0;
				for (i_ = sx; i_ <= fx; i_++)
				{
					v += dx[i_] * dy[i_ + i1_];
				}
				result = v;
				return result;
			}

			private void lbfgsmcsrch(int n,
				ref double[] x,
				ref double f,
				ref double[] g,
				ref double[] s,
				int sstart,
				ref double stp,
				double ftol,
				double xtol,
				int maxfev,
				ref int info,
				ref int nfev,
				ref double[] wa,
				double gtol,
				double stpmin,
				double stpmax)
			{
				int infoc = 0;
				int j = 0;
				bool brackt = new bool();
				bool stage1 = new bool();
				double dg = 0;
				double dgm = 0;
				double dginit = 0;
				double dgtest = 0;
				double dgx = 0;
				double dgxm = 0;
				double dgy = 0;
				double dgym = 0;
				double finit = 0;
				double ftest1 = 0;
				double fm = 0;
				double fx = 0;
				double fxm = 0;
				double fy = 0;
				double fym = 0;
				double p5 = 0;
				double p66 = 0;
				double stx = 0;
				double sty = 0;
				double stmin = 0;
				double stmax = 0;
				double width = 0;
				double width1 = 0;
				double xtrapf = 0;
				double zero = 0;
				double mytemp = 0;

				sstart = sstart - 1;
				p5 = 0.5;
				p66 = 0.66;
				xtrapf = 4.0;
				zero = 0;
				funcgrad(ref x, ref f, ref g);
				infoc = 1;
				info = 0;
				if (n <= 0 | stp <= 0 | ftol < 0 | gtol < zero | xtol < zero | stpmin < zero | stpmax < stpmin | maxfev <= 0)
				{
					return;
				}
				dginit = 0;
				for (j = 1; j <= n; j++)
				{
					dginit = dginit + g[j] * s[j + sstart];
				}
				if (dginit >= 0)
				{
					return;
				}
				brackt = false;
				stage1 = true;
				nfev = 0;
				finit = f;
				dgtest = ftol * dginit;
				width = stpmax - stpmin;
				width1 = width / p5;
				for (j = 1; j <= n; j++)
				{
					wa[j] = x[j];
				}
				stx = 0;
				fx = finit;
				dgx = dginit;
				sty = 0;
				fy = finit;
				dgy = dginit;
				while (true)
				{
					if (brackt)
					{
						if (stx < sty)
						{
							stmin = stx;
							stmax = sty;
						}
						else
						{
							stmin = sty;
							stmax = stx;
						}
					}
					else
					{
						stmin = stx;
						stmax = stp + xtrapf * (stp - stx);
					}
					if (stp > stpmax)
					{
						stp = stpmax;
					}
					if (stp < stpmin)
					{
						stp = stpmin;
					}
					if (brackt & (stp <= stmin | stp >= stmax) | nfev >= maxfev - 1 | infoc == 0 | brackt & stmax - stmin <= xtol * stmax)
					{
						stp = stx;
					}
					for (j = 1; j <= n; j++)
					{
						x[j] = wa[j] + stp * s[j + sstart];
					}
					funcgrad(ref x, ref f, ref g);
					info = 0;
					nfev = nfev + 1;
					dg = 0;
					for (j = 1; j <= n; j++)
					{
						dg = dg + g[j] * s[j + sstart];
					}
					ftest1 = finit + stp * dgtest;
					if (brackt & (stp <= stmin | stp >= stmax) | infoc == 0)
					{
						info = 6;
					}
					if (stp == stpmax & f <= ftest1 & dg <= dgtest)
					{
						info = 5;
					}
					if (stp == stpmin & (f > ftest1 | dg >= dgtest))
					{
						info = 4;
					}
					if (nfev >= maxfev)
					{
						info = 3;
					}
					if (brackt & stmax - stmin <= xtol * stmax)
					{
						info = 2;
					}
					if (f <= ftest1 & Math.Abs(dg) <= -(gtol * dginit))
					{
						info = 1;
					}
					if (info != 0)
					{
						return;
					}
					mytemp = ftol;
					if (gtol < ftol)
					{
						mytemp = gtol;
					}
					if (stage1 & f <= ftest1 & dg >= mytemp * dginit)
					{
						stage1 = false;
					}
					if (stage1 & f <= fx & f > ftest1)
					{
						fm = f - stp * dgtest;
						fxm = fx - stx * dgtest;
						fym = fy - sty * dgtest;
						dgm = dg - dgtest;
						dgxm = dgx - dgtest;
						dgym = dgy - dgtest;
						lbfgsmcstep(ref stx, ref fxm, ref dgxm, ref sty, ref fym, ref dgym, ref stp, fm, dgm, ref brackt, stmin, stmax, ref infoc);
						fx = fxm + stx * dgtest;
						fy = fym + sty * dgtest;
						dgx = dgxm + dgtest;
						dgy = dgym + dgtest;
					}
					else
					{
						lbfgsmcstep(ref stx, ref fx, ref dgx, ref sty, ref fy, ref dgy, ref stp, f, dg, ref brackt, stmin, stmax, ref infoc);
					}
					if (brackt)
					{
						if (Math.Abs(sty - stx) >= p66 * width1)
						{
							stp = stx + p5 * (sty - stx);
						}
						width1 = width;
						width = Math.Abs(sty - stx);
					}
				}
			}

			private void lbfgsmcstep(ref double stx,
				ref double fx,
				ref double dx,
				ref double sty,
				ref double fy,
				ref double dy,
				ref double stp,
				double fp,
				double dp,
				ref bool brackt,
				double stmin,
				double stmax,
				ref int info)
			{
				bool bound = new bool();
				double gamma = 0;
				double p = 0;
				double q = 0;
				double r = 0;
				double s = 0;
				double sgnd = 0;
				double stpc = 0;
				double stpf = 0;
				double stpq = 0;
				double theta = 0;

				info = 0;
				if (brackt & (stp <= Math.Min(stx, sty) | stp >= Math.Max(stx, sty)) | dx * (stp - stx) >= 0 | stmax < stmin)
				{
					return;
				}
				sgnd = dp * (dx / Math.Abs(dx));
				if (fp > fx)
				{
					info = 1;
					bound = true;
					theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
					s = Math.Max(Math.Abs(theta), Math.Max(Math.Abs(dx), Math.Abs(dp)));
					gamma = s * Math.Sqrt(Sqr(theta / s) - dx / s * (dp / s));
					if (stp < stx)
					{
						gamma = -gamma;
					}
					p = gamma - dx + theta;
					q = gamma - dx + gamma + dp;
					r = p / q;
					stpc = stx + r * (stp - stx);
					stpq = stx + dx / ((fx - fp) / (stp - stx) + dx) / 2 * (stp - stx);
					if (Math.Abs(stpc - stx) < Math.Abs(stpq - stx))
					{
						stpf = stpc;
					}
					else
					{
						stpf = stpc + (stpq - stpc) / 2;
					}
					brackt = true;
				}
				else
				{
					if (sgnd < 0)
					{
						info = 2;
						bound = false;
						theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
						s = Math.Max(Math.Abs(theta), Math.Max(Math.Abs(dx), Math.Abs(dp)));
						gamma = s * Math.Sqrt(Sqr(theta / s) - dx / s * (dp / s));
						if (stp > stx)
						{
							gamma = -gamma;
						}
						p = gamma - dp + theta;
						q = gamma - dp + gamma + dx;
						r = p / q;
						stpc = stp + r * (stx - stp);
						stpq = stp + dp / (dp - dx) * (stx - stp);
						if (Math.Abs(stpc - stp) > Math.Abs(stpq - stp))
						{
							stpf = stpc;
						}
						else
						{
							stpf = stpq;
						}
						brackt = true;
					}
					else
					{
						if (Math.Abs(dp) < Math.Abs(dx))
						{
							info = 3;
							bound = true;
							theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
							s = Math.Max(Math.Abs(theta), Math.Max(Math.Abs(dx), Math.Abs(dp)));
							gamma = s * Math.Sqrt(Math.Max(0, Sqr(theta / s) - dx / s * (dp / s)));
							if (stp > stx)
							{
								gamma = -gamma;
							}
							p = gamma - dp + theta;
							q = gamma + (dx - dp) + gamma;
							r = p / q;
							if (r < 0 & gamma != 0)
							{
								stpc = stp + r * (stx - stp);
							}
							else
							{
								if (stp > stx)
								{
									stpc = stmax;
								}
								else
								{
									stpc = stmin;
								}
							}
							stpq = stp + dp / (dp - dx) * (stx - stp);
							if (brackt)
							{
								if (Math.Abs(stp - stpc) < Math.Abs(stp - stpq))
								{
									stpf = stpc;
								}
								else
								{
									stpf = stpq;
								}
							}
							else
							{
								if (Math.Abs(stp - stpc) > Math.Abs(stp - stpq))
								{
									stpf = stpc;
								}
								else
								{
									stpf = stpq;
								}
							}
						}
						else
						{
							info = 4;
							bound = false;
							if (brackt)
							{
								theta = 3 * (fp - fy) / (sty - stp) + dy + dp;
								s = Math.Max(Math.Abs(theta), Math.Max(Math.Abs(dy), Math.Abs(dp)));
								gamma = s * Math.Sqrt(Sqr(theta / s) - dy / s * (dp / s));
								if (stp > sty)
								{
									gamma = -gamma;
								}
								p = gamma - dp + theta;
								q = gamma - dp + gamma + dy;
								r = p / q;
								stpc = stp + r * (sty - stp);
								stpf = stpc;
							}
							else
							{
								if (stp > stx)
								{
									stpf = stmax;
								}
								else
								{
									stpf = stmin;
								}
							}
						}
					}
				}
				if (fp > fx)
				{
					sty = stp;
					fy = fp;
					dy = dp;
				}
				else
				{
					if (sgnd < 0.0)
					{
						sty = stx;
						fy = fx;
						dy = dx;
					}
					stx = stp;
					fx = fp;
					dx = dp;
				}
				stpf = Math.Min(stmax, stpf);
				stpf = Math.Max(stmin, stpf);
				stp = stpf;
				if (brackt & bound)
				{
					if (sty > stx)
					{
						stp = Math.Min(stx + 0.66 * (sty - stx), stp);
					}
					else
					{
						stp = Math.Max(stx + 0.66 * (sty - stx), stp);
					}
				}
			}

			private void lbfgsnewiteration(ref double[] x,
				double f,
				ref double[] g)
			{
			}
		}
	}
}
