﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System;
using FuncLib.Optimization.Ipopt;
using FuncLib.Optimization.Ipopt.Cureos;

namespace FuncLibGH.BenchSimpleElastica
{
    // Copyright (C) 2010 Anders Gustafsson and others. All Rights Reserved.
    // This code is published under the Eclipse Public License.
    //
    // Author:  Anders Gustafsson, Cureos AB 2010-12-07

    public class Program
    {
        public static void Main(string[] args)
        {
            /* create the IpoptProblem */
            HS071 p = new HS071();

            /* allocate space for the initial point and set the values */
            double[] x = { 1.0, 5.0, 5.0, 1.0 };

            IpoptReturnCode status;

            using (var problem = new IpoptProblem(p._n, p._x_L, p._x_U, p._m, p._g_L, p._g_U, p._nele_jac, p._nele_hess,
                p.eval_f, p.eval_g, p.eval_grad_f, p.eval_jac_g, p.eval_h))
            {
                /* Set some options.  The following ones are only examples,
                    they might not be suitable for your problem. */
                problem.AddOption("derivative_test", "second-order");
                problem.AddOption("tol", 1e-7);
                problem.AddOption("mu_strategy", "adaptive");
                problem.AddOption("output_file", "hs071.txt");

#if INTERMEDIATE
            problem.SetIntermediateCallback(p.intermediate);
#endif
                /* solve the problem */
                double obj;
                status = problem.SolveProblem(x, out obj, null, null, null, null);
            }

            Console.WriteLine("{0}{0}Optimization return status: {1}{0}{0}", Environment.NewLine, status);

            for (int i = 0; i < 4; ++i) Console.WriteLine("x[{0}]={1}", i, x[i]);

            Console.WriteLine("{0}Press <RETURN> to exit...", Environment.NewLine);
            Console.ReadLine();
        }
    }

    public class HS071
    {
        public int _n;
        public int _m;
        public int _nele_jac;
        public int _nele_hess;
        public double[] _x_L;
        public double[] _x_U;
        public double[] _g_L;
        public double[] _g_U;

        public HS071()
        {
            /* set the number of variables and allocate space for the bounds */
            /* set the values for the variable bounds */
            _n = 4;
            _x_L = new double[] { 1.0, 1.0, 1.0, 1.0 };
            _x_U = new double[] { 5.0, 5.0, 5.0, 5.0 };

            /* set the number of constraints and allocate space for the bounds */
            _m = 2;

            /* set the values of the constraint bounds */
            _g_L = new double[] { 25.0, 40.0 };
            _g_U = new double[] { IpoptProblem.PositiveInfinity, 40.0 };

            /* Number of nonzeros in the Jacobian of the constraints */
            _nele_jac = 8;

            /* Number of nonzeros in the Hessian of the Lagrangian (lower or
                upper triangual part only) */
            _nele_hess = 10;
        }

        public bool eval_f(int n, double[] x, bool new_x, out double obj_value)
        {
            obj_value = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];

            return true;
        }

        public bool eval_grad_f(int n, double[] x, bool new_x, double[] grad_f)
        {
            grad_f[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
            grad_f[1] = x[0] * x[3];
            grad_f[2] = x[0] * x[3] + 1;
            grad_f[3] = x[0] * (x[0] + x[1] + x[2]);

            return true;
        }

        public bool eval_g(int n, double[] x, bool new_x, int m, double[] g)
        {
            g[0] = x[0] * x[1] * x[2] * x[3];
            g[1] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3];

            return true;
        }

        public bool eval_jac_g(int n, double[] x, bool new_x, int m, int nele_jac, int[] iRow, int[] jCol, double[] values)
        {
            if (values == null)
            {
                /* set the structure of the jacobian */
                /* this particular jacobian is dense */

                iRow[0] = 0;
                jCol[0] = 0;
                iRow[1] = 0;
                jCol[1] = 1;
                iRow[2] = 0;
                jCol[2] = 2;
                iRow[3] = 0;
                jCol[3] = 3;
                iRow[4] = 1;
                jCol[4] = 0;
                iRow[5] = 1;
                jCol[5] = 1;
                iRow[6] = 1;
                jCol[6] = 2;
                iRow[7] = 1;
                jCol[7] = 3;
            }
            else
            {
                /* return the values of the jacobian of the constraints */

                values[0] = x[1] * x[2] * x[3]; /* 0,0 */
                values[1] = x[0] * x[2] * x[3]; /* 0,1 */
                values[2] = x[0] * x[1] * x[3]; /* 0,2 */
                values[3] = x[0] * x[1] * x[2]; /* 0,3 */

                values[4] = 2 * x[0];         /* 1,0 */
                values[5] = 2 * x[1];         /* 1,1 */
                values[6] = 2 * x[2];         /* 1,2 */
                values[7] = 2 * x[3];         /* 1,3 */
            }

            return true;
        }

        public bool eval_h(int n, double[] x, bool new_x, double obj_factor,
                    int m, double[] lambda, bool new_lambda,
                    int nele_hess, int[] iRow, int[] jCol,
                    double[] values)
        {
            if (values == null)
            {
                /* set the Hessian structure. This is a symmetric matrix, fill the lower left
                    * triangle only. */

                /* the hessian for this problem is actually dense */
                int idx = 0; /* nonzero element counter */
                for (int row = 0; row < 4; row++)
                {
                    for (int col = 0; col <= row; col++)
                    {
                        iRow[idx] = row;
                        jCol[idx] = col;
                        idx++;
                    }
                }

            }
            else
            {
                /* return the values. This is a symmetric matrix, fill the lower left
                    * triangle only */

                /* fill the objective portion */
                values[0] = obj_factor * (2 * x[3]);               /* 0,0 */

                values[1] = obj_factor * (x[3]);                 /* 1,0 */
                values[2] = 0;                                   /* 1,1 */

                values[3] = obj_factor * (x[3]);                 /* 2,0 */
                values[4] = 0;                                   /* 2,1 */
                values[5] = 0;                                   /* 2,2 */

                values[6] = obj_factor * (2 * x[0] + x[1] + x[2]); /* 3,0 */
                values[7] = obj_factor * (x[0]);                 /* 3,1 */
                values[8] = obj_factor * (x[0]);                 /* 3,2 */
                values[9] = 0;                                   /* 3,3 */


                /* add the portion for the first constraint */
                values[1] += lambda[0] * (x[2] * x[3]);          /* 1,0 */

                values[3] += lambda[0] * (x[1] * x[3]);          /* 2,0 */
                values[4] += lambda[0] * (x[0] * x[3]);          /* 2,1 */

                values[6] += lambda[0] * (x[1] * x[2]);          /* 3,0 */
                values[7] += lambda[0] * (x[0] * x[2]);          /* 3,1 */
                values[8] += lambda[0] * (x[0] * x[1]);          /* 3,2 */

                /* add the portion for the second constraint */
                values[0] += lambda[1] * 2;                      /* 0,0 */

                values[2] += lambda[1] * 2;                      /* 1,1 */

                values[5] += lambda[1] * 2;                      /* 2,2 */

                values[9] += lambda[1] * 2;                      /* 3,3 */
            }

            return true;
        }

#if INTERMEDIATE
    public bool intermediate(IpoptAlgorithmMode alg_mod, int iter_count, double obj_value, double inf_pr, double inf_du,
        double mu, double d_norm, double regularization_size, double alpha_du, double alpha_pr, int ls_trials)
    {
        Console.WriteLine("Intermediate callback method at iteration {0} in {1} with d_norm {2}",
            iter_count, alg_mod, d_norm);
        return iter_count < 5;
    }
#endif

    }
}
