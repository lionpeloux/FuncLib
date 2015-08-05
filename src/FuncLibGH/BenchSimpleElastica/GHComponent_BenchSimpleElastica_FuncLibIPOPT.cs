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

namespace FuncLibGH
{
    public class GHComponent_BenchSimpleElastica_FuncLibIPOPT : GH_Component
    {
        // FIELDS
        private bool reset = true;
        private bool reset_cache = true;
        private bool enable_quasi_newton;
        private bool enable_inextensibility_as_constraint;
        private Curve Xi_crv;
        private Curve X_bar_crv;
        private Polyline Xi, X_bar;

        private SimpleElastica_FuncLib_IPOPT simpleElastica_funclib_ipopt;


        public GHComponent_BenchSimpleElastica_FuncLibIPOPT()
            : base("FuncLib IPOPT","FuncLib IPOPT","Simple Elastica using FuncLib IPOPT solver", "OptiELA", "Simple Elastica Bench")
        {
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("{9074C9D4-0830-414E-80CF-B13CF1A07042}"); }
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Xi", "Xi", "initial position", GH_ParamAccess.item);
            pManager.AddCurveParameter("X_bar", "X_bar", "initial position", GH_ParamAccess.item);
            pManager.AddBooleanParameter("enableQuasiNewton", "enableQuasiNewton", "enable Quasi-Newton approximation of the hessian", GH_ParamAccess.item);
            pManager.AddBooleanParameter("enableInextensibilityAsConstraint", "enableInextensibilityAsConstraint", "enable inextensibility to be trated as constrained and not axial enery", GH_ParamAccess.item);
            pManager.AddBooleanParameter("reset", "reset", "restart the solver", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.Register_StringParam("info", "info", "solver info");
            pManager.Register_PointParam("X", "X", "results");
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            if (!DA.GetData(0, ref Xi_crv)) { return; }
            if (!DA.GetData(1, ref X_bar_crv)) { return; }
            if (!DA.GetData(2, ref enable_quasi_newton)) { return; }
            if (!DA.GetData(3, ref enable_inextensibility_as_constraint)) { return; }
            if (!DA.GetData(4, ref reset)) { return; }

            if (!Xi_crv.TryGetPolyline(out Xi)) return;
            if (!X_bar_crv.TryGetPolyline(out X_bar)) return;

            simpleElastica_funclib_ipopt = new SimpleElastica_FuncLib_IPOPT(X_bar, Xi, enable_inextensibility_as_constraint);
            simpleElastica_funclib_ipopt.opt.EnableQuasiNewton = enable_quasi_newton;
            simpleElastica_funclib_ipopt.Run();
            foreach (IpoptOption option in simpleElastica_funclib_ipopt.opt.Options)
	        {
                Rhino.RhinoApp.WriteLine(option.Name + " = " + option.Value);
	        }

            
            DA.SetData(0, simpleElastica_funclib_ipopt.info);
            DA.SetDataList(1, simpleElastica_funclib_ipopt.X);
        }

        public class SimpleElastica_FuncLib_IPOPT : SimpleElastica
        {
            // FuncLib IPOPT Variables
            Variable[] X_ipopt;
            Variable[] x, y, z;
            VariableAssignment[] Xi_ipopt;

            // FuncLib IPOPT Solver
            public IpoptOptimizer opt;
            PreparedOptimizer opt_prepared;
            IOptimizerResult opt_res;
            List<IpoptIntermediateEventArgs> intermediateEventsList;

            public SimpleElastica_FuncLib_IPOPT(Polyline X_bar, Polyline Xi, bool enable_inextensibility_as_constraint)
                : base(X_bar, Xi)
            {
                Rhino.RhinoApp.WriteLine("X_bar : " + X_bar.Length);
                Rhino.RhinoApp.WriteLine("X_i : " + Xi.Length);

                // FuncLib definition of variables and assignement       
                X_ipopt = new Variable[3 * (n + 1)];
                x = new Variable[n + 1]; y = new Variable[n + 1]; z = new Variable[n + 1];
                for (int i = 0; i < n + 1; i++)
                {
                    x[i] = new Variable("x_" + i);
                    y[i] = new Variable("y_" + i);
                    z[i] = new Variable("z_" + i);
                    X_ipopt[3 * i] = x[i];
                    X_ipopt[3 * i + 1] = y[i];
                    X_ipopt[3 * i + 2] = z[i];
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
                    Eb += 0.5 * EI * K2[i] / ((E_bar[i - 1].Length + E_bar[i].Length) / 2);
                }

                // Constraints : supports
                //FunctionConstraint[] C_supports = new FunctionConstraint[6];
                //C_supports[0] = new FunctionConstraint(x[0], Xi[0].X, Xi[0].X);
                //C_supports[1] = new FunctionConstraint(y[0], Xi[0].Y, Xi[0].Y);
                //C_supports[2] = new FunctionConstraint(z[0], Xi[0].Z, Xi[0].Z);
                //C_supports[3] = new FunctionConstraint(x[n], Xi[n].X, Xi[n].X);
                //C_supports[4] = new FunctionConstraint(y[n], Xi[n].Y, Xi[n].Y);
                //C_supports[5] = new FunctionConstraint(z[n], Xi[n].Z, Xi[n].Z);

                // Constraints : suppport as 
                VariableConstraint[] C_supports = new VariableConstraint[6];
                C_supports[0] = new VariableConstraint(x[0], Xi[0].X, Xi[0].X);
                C_supports[1] = new VariableConstraint(y[0], Xi[0].Y, Xi[0].Y);
                C_supports[2] = new VariableConstraint(z[0], Xi[0].Z, Xi[0].Z);
                C_supports[3] = new VariableConstraint(x[n], Xi[n].X, Xi[n].X);
                C_supports[4] = new VariableConstraint(y[n], Xi[n].Y, Xi[n].Y);
                C_supports[5] = new VariableConstraint(z[n], Xi[n].Z, Xi[n].Z);

                // Constraints : inextensibility
                FunctionConstraint[] C_inextensibility = new FunctionConstraint[n];
                for (int i = 0; i < n; i++)
                {
                    C_inextensibility[i] = new FunctionConstraint(
                        Function.Pow(x[i + 1] - x[i], 2) + Function.Pow(y[i + 1] - y[i], 2) + Function.Pow(z[i + 1] - z[i], 2),
                        Math.Pow(E_bar[i].Length, 2), Math.Pow(E_bar[i].Length, 2));
                }

                // setup a new FuncLib ipopt solver
                opt = new IpoptOptimizer();
                intermediateEventsList = new List<IpoptIntermediateEventArgs>();
                opt.Intermediate += new EventHandler<IpoptIntermediateEventArgs>(IntermediateOutput);
                opt.PrintLevel = 5;
                opt.MaxIterations = 1000;
                opt.EnableQuasiNewton = true;
                opt.Variables.Add(X_ipopt);

                if (enable_inextensibility_as_constraint)
                {
                    opt.ObjectiveFunction = Eb;
                    opt.Constraints.Add(C_supports);
                    opt.Constraints.Add(C_inextensibility);
                }
                else
                {
                    opt.ObjectiveFunction = Ea + Eb;
                    opt.Constraints.Add(C_supports);
                }

                // prepare the solver for multiple run
                opt_prepared = opt.Prepare();
            }

            public override void Run()
            {
                intermediateEventsList.Clear();
                Xi_ipopt = new VariableAssignment[3 * (n + 1)];
                for (int i = 0; i < n + 1; i++)
                {
                    Xi_ipopt[3 * i] = new VariableAssignment(x[i], base.Xi[i].X);
                    Xi_ipopt[3 * i + 1] = new VariableAssignment(y[i], base.Xi[i].Y);
                    Xi_ipopt[3 * i + 2] = new VariableAssignment(z[i], base.Xi[i].Z);
                }

                opt_prepared = opt.Prepare();
                opt_res = opt_prepared.Run(Xi_ipopt);

                for (int i = 0; i < n + 1; i++)
                {
                    base.X[i] = new Vector3d(opt_res.OptimalPoint[x[i]], opt_res.OptimalPoint[y[i]], opt_res.OptimalPoint[z[i]]);
                }

                info = intermediateEventsList.ToStringTable(
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
                        e => e.LineSearchStepsCount
                );

                info += System.Environment.NewLine;
                info += "OPTIMAL VALUE = " + opt_res.OptimalValue + System.Environment.NewLine;
                info += "CVG STATUT = " + opt_res.Status + System.Environment.NewLine;
            }
            public void IntermediateOutput(Object sender, IpoptIntermediateEventArgs e)
            {
                // store the event with current iteration informations
                intermediateEventsList.Add(e);
            }
        }

    }
}
