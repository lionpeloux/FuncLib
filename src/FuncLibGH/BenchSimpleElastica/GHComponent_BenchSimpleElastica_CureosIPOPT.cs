using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using FuncLib.Functions;
using FuncLib.Optimization.Ipopt;
using FuncLib.Optimization;
using FuncLibGH;
using FuncLib.Utilities;
using FuncLibGH.BenchSimpleElastica;
using FuncLib.Optimization.Ipopt.Cureos;

namespace FuncLibGH
{
    public class GHComponent_BenchSimpleElastica_CureosIPOPT : GH_Component
    {
        // FIELDS
        private int iter_max;
        private bool reset = true;
        private bool reset_cache = true;
        private Curve Xi_crv;
        private Curve X_bar_crv;
        private Polyline Xi, X_bar;

        private SimpleElastica_Cureos_IPOPT simpleElastica_cureos_ipopt;


        public GHComponent_BenchSimpleElastica_CureosIPOPT()
            : base("Cureos IPOPT","Creos IPOPT","Simple Elastica using Cureos IPOPT solver", "OptiELA", "Simple Elastica Bench")
        {
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("{D13D932B-D126-4995-9FA8-06F72F2F384C}"); }
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Xi", "Xi", "initial position", GH_ParamAccess.item);
            pManager.AddCurveParameter("X_bar", "X_bar", "initial position", GH_ParamAccess.item);
            pManager.AddBooleanParameter("reset", "reset", "restart the solver", GH_ParamAccess.item);
            pManager.AddIntegerParameter("iter_max", "iter_max", "iter_max", GH_ParamAccess.item,0);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.Register_StringParam("info", "info", "solver info");
            pManager.Register_PointParam("X", "X", "results");
            pManager.Register_VectorParam("e", "e", "edges");
            pManager.Register_DoubleParam("e_bar_norm", "e_bar_norm", "e norm au repos");
            pManager.Register_DoubleParam("l_bar", "l_bar", "longueur associée au noeud i");
            pManager.Register_VectorParam("grad_f", "grad_f", "last gradient");
            pManager.Register_VectorParam("gradk_dot_k_prev", "gradk_dot_k_prev", "last gradient");
            pManager.Register_VectorParam("gradk_dot_k_next", "gradk_dot_k_next", "last gradient");
            pManager.Register_VectorParam("k", "k", "curvature");
            pManager.Register_DoubleParam("k2", "k2", "square curvature norm");
            pManager.Register_VectorParam("check_grad_f", "check_grad_f", "check gradient");
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            if (!DA.GetData(0, ref Xi_crv)) { return; }
            if (!DA.GetData(1, ref X_bar_crv)) { return; }
            if (!DA.GetData(2, ref reset)) { return; }
            if (!DA.GetData(3, ref iter_max)) { return; }

            if (!Xi_crv.TryGetPolyline(out Xi)) return;
            if (!X_bar_crv.TryGetPolyline(out X_bar)) return;


            simpleElastica_cureos_ipopt = new SimpleElastica_Cureos_IPOPT(X_bar, Xi);

            simpleElastica_cureos_ipopt.iter_max = iter_max;
            simpleElastica_cureos_ipopt.Run();

            int j;

            Vector3d[] _e = new Vector3d[Xi.Count-1];
            for (int i = 0; i < Xi.Count - 1; i++)
            {
                j = 3 * i;
                _e[i] = new Vector3d(
                    simpleElastica_cureos_ipopt.e[j],
                    simpleElastica_cureos_ipopt.e[j + 1],
                    simpleElastica_cureos_ipopt.e[j + 2]);
            }

            Vector3d[] grad_f = new Vector3d[Xi.Count];
            Vector3d[] chek_grad_f = new Vector3d[Xi.Count];
            Vector3d[] gradk_dot_k_prev = new Vector3d[Xi.Count];
            Vector3d[] gradk_dot_k_next = new Vector3d[Xi.Count];
            Vector3d[] k = new Vector3d[Xi.Count];
            for (int i = 0; i < Xi.Count; i++)
            {
                j = 3 * i;
                grad_f[i] = new Vector3d(
                    simpleElastica_cureos_ipopt._grad_f[j],
                    simpleElastica_cureos_ipopt._grad_f[j + 1],
                    simpleElastica_cureos_ipopt._grad_f[j + 2]);

                chek_grad_f[i] = new Vector3d(
                    simpleElastica_cureos_ipopt._check_grad_f[j],
                    simpleElastica_cureos_ipopt._check_grad_f[j + 1],
                    simpleElastica_cureos_ipopt._check_grad_f[j + 2]);

                gradk_dot_k_prev[i] = new Vector3d(
                    simpleElastica_cureos_ipopt.gradk_dot_k_prev[j],
                    simpleElastica_cureos_ipopt.gradk_dot_k_prev[j + 1],
                    simpleElastica_cureos_ipopt.gradk_dot_k_prev[j + 2]);

                gradk_dot_k_next[i] = new Vector3d(
                    simpleElastica_cureos_ipopt.gradk_dot_k_next[j],
                    simpleElastica_cureos_ipopt.gradk_dot_k_next[j + 1],
                    simpleElastica_cureos_ipopt.gradk_dot_k_next[j + 2]);

                k[i] = new Vector3d(
                    simpleElastica_cureos_ipopt.k[j],
                    simpleElastica_cureos_ipopt.k[j + 1],
                    simpleElastica_cureos_ipopt.k[j + 2]);
                
            }

            string res = simpleElastica_cureos_ipopt.intermediateEventsList.ToStringTable(
            new string[] { "it", "objective", "inf_pr", "inf_du", "lg(mu)", "d", "lg(rg)", "alpha_pr", "alpha_du", "ls" },
            e => e.Iteration,
            e => e.ObjectiveValue.ToString("E1"),
            e => e.PrimalInfeasibility.ToString("E1"),
            e => e.DualInfeasibility.ToString("E1"),
            e => e.BarrierParameter.ToString("E1"),
            e => e.PrimalStep.ToString("E1"),
            e => e.RegularizationTerm.ToString(""),
            e => e.PriamlStepSize.ToString(""),
            e => e.DualStepSize.ToString(""),
            e => e.LineSearchStepsCount);

            res += System.Environment.NewLine;
            res += "CVG STATUT = " + simpleElastica_cureos_ipopt.status + System.Environment.NewLine;

            DA.SetData(0, res);
            DA.SetDataList(1, simpleElastica_cureos_ipopt.X);
            DA.SetDataList(2, _e);
            DA.SetDataList(3, simpleElastica_cureos_ipopt.e_bar_norm);
            DA.SetDataList(4, simpleElastica_cureos_ipopt.l_bar);
            DA.SetDataList(5, grad_f);
            DA.SetDataList(6, gradk_dot_k_prev);
            DA.SetDataList(7, gradk_dot_k_next);
            DA.SetDataList(8, k);
            DA.SetDataList(9, simpleElastica_cureos_ipopt.k2);
            DA.SetDataList(10, chek_grad_f);
        }

        public class SimpleElastica_Cureos_IPOPT : SimpleElastica
        {
            public int iter_max;
            public int _n;
            public int _m;
            public int _nele_jac;
            public int _nele_hess;
            public double[] _x_L;
            public double[] _x_U;
            public double[] _g_L;
            public double[] _g_U;
            public double[] _grad_f;

            public IpoptProblem problem;
            public List<IpoptIntermediateEventArgs> intermediateEventsList;
            public IpoptReturnCode status;

            // INTERNAL DATA
            int j;
            public double[] e;              // ei = xi+1 - xi | i = 0 ... n-1 (tableau de dim 3n)
            public double[] e_bar_norm;     // |ei| | i = 0 ... n-1 (tableau de dim n)
            public double[] l_bar;          // |ei-1|/2 + |ei|/2 | i = 0 ... n-1 (tableau de dim n+1)

            public double[] k;              // ki | xi , i = 0 ... n (tableau de dim 3(n+1))
            public double[] k_e_prod;       // |ei-1||ei| | xi , i = 1 ... n-1 & 0 en i = 0 , n (tableau de dim n+1)
            public double[] k_denom;        // |ei-1||ei| + ei-1.ei | xi , i = 0 ... n (tableau de dim n+1)
            public double[] k_num;          // 2 ei-1 x ei | xi , i = 0 ... n (tableau de dim 3(n+1))

            public double[] k2;             // ki^2 = ki.ki | xi , i = 0 ... n (tableau de dim n+1)

            public double[] gradk_dot_k_prev;      // gradk_dot_k_prev = (EI/li) x [grad_i-1(ki)]^T.ki est un vecteur de R3 (tableau de dim 3(n+1))
            public double[] gradk_dot_k_next;      // gradk_dot_k_next = (EI/li) x [grad_i+1(ki)]^T.ki est un vecteur de R3 (tableau de dim 3(n+1))


            // CHEK EVAL F
            public Function Eb;
            public Function[] grad_Eb;
            public double[] _check_grad_f;
            VariableAssignment[] Xi_ipopt;
            Variable[] X_ipopt;
            Variable[] x;
            Variable[] y;
            Variable[] z;

            public SimpleElastica_Cureos_IPOPT(Polyline X_bar, Polyline Xi)
                : base(X_bar, Xi)
            {
                Rhino.RhinoApp.WriteLine("CUREOS IPOPT : " + "START INIT");
                intermediateEventsList = new List<IpoptIntermediateEventArgs>();
                /* set the number of variables and allocate space for the bounds */
                /* set the values for the variable bounds */
                _n = 3 * (n + 1);
                Rhino.RhinoApp.WriteLine("CUREOS IPOPT : " + "_n = " + _n);
                Rhino.RhinoApp.WriteLine("CUREOS IPOPT : " + "_X_L ...");
                Rhino.RhinoApp.WriteLine("CUREOS IPOPT : " + "_X_U ...");
                _x_L = new double[_n];
                _x_U = new double[_n];

                _x_L[0] = Xi[0].X;
                _x_L[1] = Xi[0].Y;
                _x_L[2] = Xi[0].Z;

                _x_U[0] = Xi[0].X;
                _x_U[1] = Xi[0].Y;
                _x_U[2] = Xi[0].Z;

                for (int i = 3; i < _n - 3; i++)
                {
                    _x_L[i] = IpoptProblem.NegativeInfinity;
                    _x_U[i] = IpoptProblem.PositiveInfinity;
                }
                _x_L[_n - 3] = Xi[n].X;
                _x_L[_n - 2] = Xi[n].Y;
                _x_L[_n - 1] = Xi[n].Z;

                _x_U[_n - 3] = Xi[n].X;
                _x_U[_n - 2] = Xi[n].Y;
                _x_U[_n - 1] = Xi[n].Z;


                /* set the number of constraints and allocate space for the bounds */
                _m = n;
                Rhino.RhinoApp.WriteLine("CUREOS IPOPT : " + "m = " + _m);

                /* set the values of the constraint bounds */
                Rhino.RhinoApp.WriteLine("CUREOS IPOPT : " + "_g_L ...");
                Rhino.RhinoApp.WriteLine("CUREOS IPOPT : " + "_g_U ...");
                _g_L = new double[n];
                _g_U = new double[n];
                for (int i = 0; i < n; i++)
                {
                    _g_L[i] = E_bar[i].Length;
                    _g_U[i] = E_bar[i].Length;
                }

                /* Number of nonzeros in the Jacobian of the constraints */
                _nele_jac = 6 * n;
                Rhino.RhinoApp.WriteLine("CUREOS IPOPT : " + "_nele_jac = " + _nele_jac);

                /* Number of nonzeros in the Hessian of the Lagrangian (lower or
                    upper triangual part only) */
                _nele_hess = 3 * (n + 1) * (1 + 3 * (n + 1)) / 2; // on suppose que tous les termes sont "non zero"
                Rhino.RhinoApp.WriteLine("CUREOS IPOPT : " + "_nele_hess = " + _nele_hess);

                // Instanciation des tableaux intermédiares don la taille est fixe :
                e = new double[3 * n];
                e_bar_norm = new double[n];
                l_bar = new double[n + 1];

                k = new double[3 * (n + 1)];
                k_e_prod = new double[n + 1];
                k_num = new double[3 * (n + 1)];
                k_denom = new double[n + 1];

                k2 = new double[n + 1];

                gradk_dot_k_prev = new double[3 * (n + 1)];
                gradk_dot_k_next = new double[3 * (n + 1)];

                Rhino.RhinoApp.WriteLine("CUREOS IPOPT : " + "e_bar_norm ...");
                // Initialisation des grandeurs imuables au cour de la simulation
                for (int i = 0; i < n; i++) { e_bar_norm[i] = E_bar[i].Length; }

                Rhino.RhinoApp.WriteLine("CUREOS IPOPT : " + "l_bar ...");
                for (int i = 1; i < n + 1; i++)
                {
                    l_bar[i] += e_bar_norm[i - 1] / 2;
                    l_bar[i - 1] += l_bar[i];
                }

                Rhino.RhinoApp.WriteLine("CUREOS IPOPT : " + "k_e_prod ...");
                k_e_prod[0] = 0; k_e_prod[n] = 0;
                for (int i = 1; i < n; i++) { k_e_prod[i] = e_bar_norm[i - 1] * e_bar_norm[i]; }

                // la coubure est toujours nulle en x_0 et x_n :
                k2[0] = 0; k2[n] = 0;
                k_num[0] = 0; k_num[1] = 0; k_num[2] = 0; k_num[3 * n + 2] = 0; k_num[3 * n + 1] = 0; k_num[3 * n] = 0;
                k_denom[0] = 1; k_denom[n] = 1;
                k[0] = 0; k[1] = 0; k[2] = 0; k[3 * n + 2] = 0; k[3 * n + 1] = 0; k[3 * n] = 0;

                check_eval_f();

                Rhino.RhinoApp.WriteLine("CUREOS IPOPT : " + "END INIT");
            }

            private void check_eval_f()
            {
                // FuncLib definition of variables and assignement       
                X_ipopt = new Variable[3 * (n + 1)];
                x = new Variable[n + 1];
                y = new Variable[n + 1];
                z = new Variable[n + 1];
                for (int i = 0; i < n + 1; i++)
                {
                    x[i] = new Variable("x_" + i);
                    y[i] = new Variable("y_" + i);
                    z[i] = new Variable("z_" + i);
                    X_ipopt[3 * i] = x[i];
                    X_ipopt[3 * i + 1] = y[i];
                    X_ipopt[3 * i + 2] = z[i];
                }

                // Bending Potential Energy
                Function[] K2 = ComputeK2(x, y, z, 1); // Eb[i] = 1.2 x EI x ki^2 , for xi tq i = 1 ... n-1

                Eb = 0;
                for (int i = 1; i < n; i++)
                {
                    Eb += 0.5 * EI * K2[i] / ((E_bar[i - 1].Length + E_bar[i].Length) / 2);
                }

                grad_Eb = new Function[3 * (n + 1)];
                for (int i = 0; i < 3 * (n + 1); i++)
                {
                    grad_Eb[i] = Eb.Derivative(X_ipopt[i]);
                }

            }
            private void eval_k(double[] x)
            {
                //Rhino.RhinoApp.WriteLine("CUREOS IPOPT : " + "eval_k ...");
                // precalcul des grandeurs necessaires aux différentes évaluations

                // ei
                for (int i = 0; i < 3 * n; i++) { e[i] = x[i + 3] - x[i]; }

                // k & k2 (Audoly) => pb de définition si ||ei-1|| <> ||ei||
                for (int i = 1; i < n; i++)
                {
                    j = 3 * i;

                    k_num[j] = 2 * (e[j - 2] * e[j + 2] - e[j - 1] * e[j + 1]);
                    k_num[j + 1] = 2 * (e[j - 1] * e[j] - e[j - 3] * e[j + 2]);
                    k_num[j + 2] = 2 * (e[j - 3] * e[j + 1] - e[j - 2] * e[j]);

                    k_denom[i] = k_e_prod[i] + e[j - 3] * e[j] + e[j - 2] * e[j + 1] + e[j - 1] * e[j + 2];

                    k[j] = k_num[j] / k_denom[i];
                    k[j + 1] = k_num[j + 1] / k_denom[i];
                    k[j + 2] = k_num[j + 2] / k_denom[i];

                    k2[i] = Math.Pow(k[j], 2) + Math.Pow(k[j + 1], 2) + Math.Pow(k[j + 2], 2);
                }

                // chek eval_f
                Xi_ipopt = new VariableAssignment[3 * (n + 1)];
                for (int i = 0; i < n + 1; i++)
                {
                    j = 3 * i;
                    Xi_ipopt[j] = new VariableAssignment(this.x[i], x[j]);
                    Xi_ipopt[j + 1] = new VariableAssignment(y[i], x[j + 1]);
                    Xi_ipopt[j + 2] = new VariableAssignment(z[i], x[j + 2]);
                }
            }

            private void eval_grad_k_dot_k()
            {
                //Rhino.RhinoApp.WriteLine("CUREOS IPOPT : " + "eval_grad_k_dot_k");
                // calculé pour pour xi , i = 1 ... n-1 et en prenant 0 pour i = 0 et i = n          

                double gamma;

                for (int i = 1; i < n; i++)
                {
                    j = 3 * i;

                    gamma = EI / (l_bar[i] * k_denom[i]); // (EI/li) / (|ei-1||ei| + ei-1.ei)

                    gradk_dot_k_prev[j + 0] = gamma * (-2 * (e[j + 1] * k[j + 2] - e[j + 2] * k[j + 1]) + k2[i] * e[j]);
                    gradk_dot_k_prev[j + 1] = gamma * (-2 * (e[j + 2] * k[j] - e[j] * k[j + 2]) + k2[i] * e[j + 1]);
                    gradk_dot_k_prev[j + 2] = gamma * (-2 * (e[j] * k[j + 1] - e[j + 1] * k[j]) + k2[i] * e[j + 2]);

                    gradk_dot_k_next[j + 0] = gamma * (-2 * (e[j - 2] * k[j + 2] - e[j - 1] * k[j + 1]) - k2[i] * e[j - 3]);
                    gradk_dot_k_next[j + 1] = gamma * (-2 * (e[j - 1] * k[j] - e[j - 3] * k[j + 2]) - k2[i] * e[j - 2]);
                    gradk_dot_k_next[j + 2] = gamma * (-2 * (e[j - 3] * k[j + 1] - e[j - 2] * k[j]) - k2[i] * e[j - 1]);
                }
            }
            public bool eval_f(int n, double[] x, bool new_x, out double obj_value)
            {
                //Rhino.RhinoApp.WriteLine("CUREOS IPOPT : " + "eval_f");
                obj_value = 0;
                for (int i = 1; i < this.n; i++)
                {
                    obj_value += k2[i] / l_bar[i];
                }
                obj_value = (EI / 2) * obj_value;
                return true;
            }
            public bool eval_g(int n, double[] x, bool new_x, int m, double[] g)
            {
                //Rhino.RhinoApp.WriteLine("CUREOS IPOPT : " + "eval_g");
                if (new_x) // le "new x" est toujours appellé par eval_g => precompute k et grad_k_dot_k
                {
                    //Rhino.RhinoApp.WriteLine("NEW X");
                    eval_k(x);
                    eval_grad_k_dot_k();
                }
                for (int i = 0; i < this.n; i++)
                {
                    j = 3 * i;
                    g[i] = Math.Pow(e[j], 2) + Math.Pow(e[j + 1], 2) + Math.Pow(e[j + 2], 2);
                }
                return true;
            }
            public bool eval_jac_g(int n, double[] x, bool new_x, int m, int nele_jac, int[] iRow, int[] jCol, double[] values)
            {
                int count = 0;
                if (values == null)
                {
                    //Rhino.RhinoApp.WriteLine("CUREOS IPOPT : " + "eval_jac_g => STRUCTURE");
                    /* set the structure of the jacobian */
                    /* this particular jacobian is dense */
                    for (int i = 0; i < this.n; i++)
                    {
                        j = 3 * i;

                        iRow[count + 0] = i; jCol[count] = j;
                        iRow[count + 1] = i; jCol[count + 1] = j + 1;
                        iRow[count + 2] = i; jCol[count + 2] = j + 2;

                        iRow[count + 3] = i; jCol[count + 3] = j + 3;
                        iRow[count + 4] = i; jCol[count + 4] = j + 4;
                        iRow[count + 5] = i; jCol[count + 5] = j + 5; ;

                        count += 6;
                    }
                }
                else
                {
                    //Rhino.RhinoApp.WriteLine("CUREOS IPOPT : " + "eval_jac_g => VALUES");
                    /* return the values of the jacobian of the constraints */
                    for (int i = 0; i < this.n; i++)
                    {
                        
                        j = 3 * i;

                        // dg_i/dx_i+1 = 2 x ei^T
                        values[count + 0] = -2 * e[j + 0];
                        //Rhino.RhinoApp.WriteLine("CUREOS IPOPT : values[" + i + "," + j + "] = -2e[" + j + "]");
                        values[count + 1] = -2 * e[j + 1];
                        //Rhino.RhinoApp.WriteLine("CUREOS IPOPT : values[" + i + "," + (j+1) + "] = -2e[" + (j+1) + "]");
                        values[count + 2] = -2 * e[j + 2];
                        //Rhino.RhinoApp.WriteLine("CUREOS IPOPT : values[" + i + "," + (j +2)+ "] = -2e[" + (j+2) + "]");

                        // dg_i/dx_i = -dg_i/dx_i+1 = - 2 x ei^T
                        values[count + 3] = 2 * e[j + 0];
                        //Rhino.RhinoApp.WriteLine("CUREOS IPOPT : values[" + i + "," + (j +3)+ "] = 2e[" + (j) + "]");
                        values[count + 4] = 2 * e[j + 1];
                        //Rhino.RhinoApp.WriteLine("CUREOS IPOPT : values[" + i + "," + (j +4)+ "] = 2e[" + (j+1) + "]");
                        values[count + 5] = 2 * e[j + 2];
                        //Rhino.RhinoApp.WriteLine("CUREOS IPOPT : values[" + i + "," + (j +5)+ "] = 2e[" + (j+2) + "]");

                        count += 6;
                    }
                }
                //Rhino.RhinoApp.WriteLine("count : " + nele_jac + " = " + count);
                // verifier que count = nele_jac à la fin de la boucle
                return true;
            }
            public bool eval_grad_f(int n, double[] x, bool new_x, double[] grad_f)
            {
                //Rhino.RhinoApp.WriteLine("CUREOS IPOPT : " + "eval_grad_f");

                j = 0;
                grad_f[j + 0] = -gradk_dot_k_next[j + 0] + gradk_dot_k_prev[j + 3];
                grad_f[j + 1] = -gradk_dot_k_next[j + 1] + gradk_dot_k_prev[j + 4];
                grad_f[j + 2] = -gradk_dot_k_next[j + 2] + gradk_dot_k_prev[j + 5];

                for (int i = 1; i < this.n; i++)
                {
                    j = 3 * i;

                    grad_f[j + 0] = gradk_dot_k_next[j - 3] - gradk_dot_k_prev[j + 0] - gradk_dot_k_next[j + 0] + gradk_dot_k_prev[j + 3];
                    grad_f[j + 1] = gradk_dot_k_next[j - 2] - gradk_dot_k_prev[j + 1] - gradk_dot_k_next[j + 1] + gradk_dot_k_prev[j + 4];
                    grad_f[j + 2] = gradk_dot_k_next[j - 1] - gradk_dot_k_prev[j + 2] - gradk_dot_k_next[j + 2] + gradk_dot_k_prev[j + 5];
                }

                j = 3 * this.n;
                grad_f[j + 0] = gradk_dot_k_next[j - 3] - gradk_dot_k_prev[j + 0];
                grad_f[j + 1] = gradk_dot_k_next[j - 2] - gradk_dot_k_prev[j + 1];
                grad_f[j + 2] = gradk_dot_k_next[j - 1] - gradk_dot_k_prev[j + 2];

                _grad_f = grad_f;
                return true;
            }
            public bool eval_h(int n, double[] x, bool new_x, double obj_factor,
                        int m, double[] lambda, bool new_lambda,
                        int nele_hess, int[] iRow, int[] jCol,
                        double[] values)
            {
                Rhino.RhinoApp.WriteLine("CUREOS IPOPT : " + "eval_h");
                // use hessian L-BFGS evaluation
                return false;
            }

            public bool intermediate(IpoptAlgorithmMode alg_mod, int iter_count, double obj_value, double inf_pr, double inf_du,
            double mu, double d_norm, double regularization_size, double alpha_du, double alpha_pr, int ls_trials)
            {

                //Rhino.RhinoApp.WriteLine("---------------------------------");
                //Rhino.RhinoApp.WriteLine("       ITERATION " + iter_count);
                //Rhino.RhinoApp.WriteLine("---------------------------------");
                //Rhino.RhinoApp.WriteLine("objective = " + obj_value);
                //Rhino.RhinoApp.WriteLine("check = " + Eb.Value(Xi_ipopt));
                //Rhino.RhinoApp.WriteLine("delta (obj) = " + Math.Abs(obj_value - Eb.Value(Xi_ipopt)));

                _check_grad_f = new double[3 * (n + 1)];
                //double norm_delta_grad = 0;
                //for (int i = 0; i < 3 * (n + 1); i++)
                //{
                //    _check_grad_f[i] = grad_Eb[i].Value(Xi_ipopt);
                //    norm_delta_grad += Math.Pow(_check_grad_f[i] - _grad_f[i], 2);
                //}
                //Rhino.RhinoApp.WriteLine("delta (grad) = " + norm_delta_grad);

                intermediateEventsList.Add(new IpoptIntermediateEventArgs(
                    alg_mod,
                    iter_count,
                    obj_value,
                    inf_pr,
                    inf_du,
                    mu,
                    d_norm,
                    regularization_size,
                    alpha_du,
                    alpha_pr,
                    ls_trials
                    ));

                //Rhino.RhinoApp.WriteLine("");
                //Rhino.RhinoApp.WriteLine("---------------------------------");
                //Rhino.RhinoApp.WriteLine("       ITERATION "+iter_count);
                //Rhino.RhinoApp.WriteLine("---------------------------------");
                info += string.Format("Intermediate callback method at iteration {0} in {1} with d_norm {2}",
                    iter_count, alg_mod, d_norm);
                info += Environment.NewLine;
                return iter_count < iter_max;
            }

            public override void Run()
            {
                /* create the IpoptProblem */
                intermediateEventsList.Clear();
                /* allocate space for the initial point and set the values */
                double[] x = new double[3 * (n + 1)];
                for (int i = 0; i < n + 1; i++)
                {
                    j = 3 * i;
                    x[j] = Xi[i].X;
                    x[j + 1] = Xi[i].Y;
                    x[j + 2] = Xi[i].Z;
                }
                eval_k(x);
                eval_grad_k_dot_k();

                using (problem = new IpoptProblem(this._n, this._x_L, this._x_U, this._m, this._g_L, this._g_U, this._nele_jac, this._nele_hess,
                    this.eval_f, this.eval_g, this.eval_grad_f, this.eval_jac_g, this.eval_h))
                {
                    /* Set some options.  The following ones are only examples,
                        they might not be suitable for your problem. */
                    problem.AddOption("derivative_test", "first-order");
                    //problem.AddOption("tol", 1e-7);
                    //problem.AddOption("mu_strategy", "adaptive");
                    problem.AddOption("output_file", "ipopt_elastica.txt");
                    problem.AddOption("hessian_approximation", "limited-memory");
                    problem.AddOption("derivative test print all", "yes");
                    problem.SetIntermediateCallback(this.intermediate);

                    /* solve the problem */
                    double obj;
                    status = problem.SolveProblem(x, out obj, null, null, null, null);
         
                }

                info += String.Format("{0}{0}Optimization return status: {1}{0}{0}", Environment.NewLine, status);

                for (int i = 0; i < n + 1; i++)
                {
                    j = 3 * i;
                    X[i].X = x[j];
                    X[i].Y = x[j + 1];
                    X[i].Z = x[j + 2];
                }
                //for (int i = 0; i < 4; ++i) Console.WriteLine("x[{0}]={1}", i, x[i]);

                //Console.WriteLine("{0}Press <RETURN> to exit...", Environment.NewLine);
                //Console.ReadLine();
            }
        }


    }
}
