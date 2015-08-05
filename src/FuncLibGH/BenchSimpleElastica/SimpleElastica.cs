using FuncLib.Utilities;
using FuncLib.Functions;
using FuncLib.Optimization;
using FuncLib.Optimization.Ipopt;
using FuncLib.Optimization.Ipopt.Cureos;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FuncLibGH.BenchSimpleElastica
{
    public abstract class SimpleElastica
    {
        protected int n;                        // number of elements
        public Vector3d[] X, Xi, X_bar;     // xi in X
        protected Vector3d[] E, Ei, E_bar;     // ei in E
        protected Vector3d[] K, Ki, K_bar;     // ki in K : curvature at node i

        protected double Ey = 1;
        protected double S = 100;
        protected double I = 1;
        protected double ES;
        protected double EI;

        public string info;
        
        protected SimpleElastica(Polyline X_bar, Polyline Xi)
        {
            n = X_bar.Count - 1; // n elements => n+1 nodes
            ES = Ey * S;
            EI = Ey * I;

            // INITIALISATION OF X, X_0, X_bar
            this.X = new Vector3d[n + 1];
            this.Xi = new Vector3d[n + 1];
            this.X_bar = new Vector3d[n + 1];
            for (int i = 0; i < n + 1; i++)
            {
                this.X_bar[i] = new Vector3d(X_bar[i]);
                this.Xi[i] = new Vector3d(Xi[i]);
                this.X[i] = new Vector3d(Xi[i]); 
            }

            // INITIALISATION OF E, E_0, E_bar
            this.E = new Vector3d[n];
            this.Ei = new Vector3d[n];
            this.E_bar = new Vector3d[n];
            for (int i = 0; i < n; i++)
            {
                this.E_bar[i] = new Vector3d(X_bar[i + 1] - X_bar[i]);
                this.Ei[i] = new Vector3d(Xi[i + 1] - Xi[i]);
                this.E[i] = new Vector3d(Xi[i + 1] - Xi[i]);
            }

            // INITIALISATION OF K, K_0, K_bar
            this.K_bar = ComputeK(E_bar);
            this.Ki = ComputeK(Ei);
            this.K = ComputeK(E);           
        }
        public abstract void Run();
        
        protected Vector3d[] ComputeK(Vector3d[] E, int type = 1)
        {
            Vector3d[] K = new Vector3d[n + 1];
            K[0] = new Vector3d(0, 0, 0);
            K[n] = new Vector3d(0, 0, 0);
            
            switch (type)
            {
                case 1: // tangent-tangent curvature (Audoly) => pb de définition si ||ei-1|| <> ||ei||
                    for (int i = 1; i < n; i++)
                    {
                        K[i] = 2 * Vector3d.CrossProduct(E[i - 1], E[i]) / (E[i - 1].Length * E[i].Length + E[i - 1] * E[i]);
                    }
                    break;
                
                case 2: // 3-pts curvature (Douthe)
                    for (int i = 1; i < n; i++)
                    {
                        K[i] = 2 * Vector3d.CrossProduct(E[i - 1], E[i]) / ((E[i - 1] + E[i]).Length);
                    }
                    break;
                
                default:
                    goto case 1;
            }
            
            return K;
        }
        protected Function[] ComputeK2(Variable[] x, Variable[] y, Variable[] z, int type = 1)
        {
            Function[] K2 = new Function[n + 1];
            K2[0] = 0;
            K2[n] = 0;
            switch (type)
            {
                case 1: // tangent-tangent curvature^2 (Audoly) => pb de définition si ||ei-1|| <> ||ei||
                    for (int i = 1; i < n; i++)
                    {
                        K2[i] = 4 * (
                            Function.Pow((y[i] - y[i - 1]) * (z[i + 1] - z[i]) - (z[i] - z[i - 1]) * (y[i + 1] - y[i]), 2) +
                            Function.Pow((z[i] - z[i - 1]) * (x[i + 1] - x[i]) - (x[i] - x[i - 1]) * (z[i + 1] - z[i]), 2) +
                            Function.Pow((x[i] - x[i - 1]) * (y[i + 1] - y[i]) - (y[i] - y[i - 1]) * (x[i + 1] - x[i]), 2));
                        K2[i] = K2[i] / Function.Pow(
                            (E_bar[i - 1].Length * E_bar[i].Length +
                            (x[i] - x[i - 1]) * (x[i + 1] - x[i]) +
                            (y[i] - y[i - 1]) * (y[i + 1] - y[i]) +
                            (z[i] - z[i - 1]) * (z[i + 1] - z[i])), 2);
                    }
                    break;

                case 2: // 3-pts curvature^2 (Douthe)
                    for (int i = 1; i < n; i++)
                    {
                        K2[i] = 4 * (
                            Function.Pow((y[i] - y[i - 1]) * (z[i + 1] - z[i]) - (z[i] - z[i - 1]) * (y[i + 1] - y[i]), 2) +
                            Function.Pow((z[i] - z[i - 1]) * (x[i + 1] - x[i]) - (x[i] - x[i - 1]) * (z[i + 1] - z[i]), 2) +
                            Function.Pow((x[i] - x[i - 1]) * (y[i + 1] - y[i]) - (y[i] - y[i - 1]) * (x[i + 1] - x[i]), 2));
                        K2[i] = K2[i] / (
                            Function.Pow(x[i + 1] - x[i - 1], 2) +
                            Function.Pow(y[i + 1] - y[i - 1], 2) +
                            Function.Pow(z[i + 1] - z[i - 1], 2)
                            );
                    }
                    break;

                default:
                    goto case 1;
            }

            return K2;
        }
    }
    
    public class SimpleElastica_FuncLib_BFGS : SimpleElastica
    {
        // FuncLib IPOPT Variables
        Variable[] X_bfgs;
        Variable[] x, y, z;
        VariableAssignment[] Xi_ipopt;

        // FuncLib IPOPT Solver
        public BfgsOptimizer opt;
        PreparedOptimizer opt_prepared;
        IOptimizerResult opt_res;

        public SimpleElastica_FuncLib_BFGS(Polyline X_bar, Polyline Xi)
            : base(X_bar, Xi)
        {
            Rhino.RhinoApp.WriteLine("X_bar : " + X_bar.Length);
            Rhino.RhinoApp.WriteLine("X_i : " + Xi.Length);

            // FuncLib definition of variables and assignement       
            X_bfgs = new Variable[3 * (n + 1)];
            x = new Variable[n + 1]; y = new Variable[n + 1]; z = new Variable[n + 1];
            for (int i = 0; i < n + 1; i++)
            {
                x[i] = new Variable("x_" + i);
                y[i] = new Variable("y_" + i);
                z[i] = new Variable("z_" + i);
                X_bfgs[3 * i] = x[i];
                X_bfgs[3 * i + 1] = y[i];
                X_bfgs[3 * i + 2] = z[i];
            }

            // Axial Potential Energy
            Function Ea = 0;
            for (int i = 0; i < n; i++) // Ea[i] = 1/2 x (ES/||ei_bar||) x (||ei||-||ei_bar||)^2 , for ei tq i = 0 ... n-1
            {
                Ea += (ES / (2 * E_bar[i].Length)) * Function.Pow(Function.Sqrt((x[i + 1] - x[i]) * (x[i + 1] - x[i]) + (y[i + 1] - y[i]) * (y[i + 1] - y[i]) + (z[i + 1] - z[i]) * (z[i + 1] - z[i])) - E_bar[i].Length, 2);
            }

            // Bending Potential Energy
            Function[] K2 = ComputeK2(x, y, z, 1); // Eb[i] = 1.2 x EI x ki^2 , for xi tq i = 1 ... n-1

            Function Eb = 0;
            for (int i = 1; i < n; i++)
            {
                Eb += EI * K2[i] / ((E_bar[i - 1].Length + E_bar[i].Length) / 2);
            }

            // Constraints : supports
            VariableEqualityConstraint[] C_supports = new VariableEqualityConstraint[6];
            C_supports[0] = new VariableEqualityConstraint(x[0], Xi[0].X);
            C_supports[1] = new VariableEqualityConstraint(y[0], Xi[0].Y);
            C_supports[2] = new VariableEqualityConstraint(z[0], Xi[0].Z);
            C_supports[3] = new VariableEqualityConstraint(x[n], Xi[n].X);
            C_supports[4] = new VariableEqualityConstraint(y[n], Xi[n].Y);
            C_supports[5] = new VariableEqualityConstraint(z[n], Xi[n].Z);

            // setup a new FuncLib BFGS solver
            opt = new BfgsOptimizer();
            opt.MaxIterations = 1000;
            opt.Variables.Add(X_bfgs);

            opt.ObjectiveFunction = Ea + Eb;
            opt.VariableEqualityConstraints.AddRange(C_supports);

            // prepare the solver for multiple run
            opt_prepared = opt.Prepare();
        }

        public override void Run()
        {
            Xi_ipopt = new VariableAssignment[3 * (n + 1)];
            for (int i = 0; i < n + 1; i++)
            {
                Xi_ipopt[3 * i] = new VariableAssignment(x[i], base.Xi[i].X);
                Xi_ipopt[3 * i + 1] = new VariableAssignment(y[i], base.Xi[i].Y);
                Xi_ipopt[3 * i + 2] = new VariableAssignment(z[i], base.Xi[i].Z);
            }

            opt_res = opt_prepared.Run(Xi_ipopt);

            for (int i = 0; i < n + 1; i++)
            {
                base.X[i] = new Vector3d(opt_res.OptimalPoint[x[i]], opt_res.OptimalPoint[y[i]], opt_res.OptimalPoint[z[i]]);
            }

            info = "";
            info += "OPTIMAL VALUE = " + opt_res.OptimalValue + System.Environment.NewLine;
            info += "CVG STATUT = " + opt_res.Status + System.Environment.NewLine;
        }
        

    }

   
}
