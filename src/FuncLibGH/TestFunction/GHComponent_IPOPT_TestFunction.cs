using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using FuncLib.Functions;
using FuncLib.Optimization.Ipopt;
using FuncLib.Optimization;
using FuncLibGH;
using FuncLib.Utilities;

namespace FuncLibGH
{
    public class GHComponent_IPOPT_TestFunction : GH_Component
    {
        // FIELDS
        private bool loop_reset = true;
        private bool loop_reset_cache = true;
        private int bfgs_iteration_max;
        private double bfgs_cvg;
        private int test_function_id, test_function_id_cache;
        private double x_start, x_start_cache;
        private double y_start, y_start_cache;
        private string info;

        Variable[] X;
        Variable x, y;
        VariableAssignment x0, y0;
        Function F;
        IpoptOptimizer opt;
        PreparedOptimizer opt_prepared;
        IOptimizerResult opt_res;
        List<IpoptIntermediateEventArgs> intermediateEventsList;

        private Vector3d X_res;

        private bool normalized;

        public GHComponent_IPOPT_TestFunction()
            : base("IPOPT over Test Function", "IPOPT",
                "Show the BFGS solver on classic test functions",
                "OptiELA", "Test Functions")
        {
            Rhino.RhinoApp.WriteLine("DEBUG A");

            x = new Variable("x");
            y = new Variable("y");
            X = new Variable[] { x, y };
            x0 = new VariableAssignment(x, x_start);
            y0 = new VariableAssignment(y, y_start);
            F = TestFunction.GetTestFunction(TestFunctionType.Sphere, X);

            opt = new IpoptOptimizer();
            intermediateEventsList = new List<IpoptIntermediateEventArgs>();
            opt.Intermediate += new EventHandler<IpoptIntermediateEventArgs>(IntermediateOutput);
            opt.PrintLevel = 5;
            opt.Variables.Add(x, y);
            opt.Constraints.Add(x >= y);
            opt.ObjectiveFunction = F;
            opt_prepared = opt.Prepare();

            intermediateEventsList.Clear();
            opt_res = opt_prepared.Run(new VariableAssignment[] { x0, y0 });

            X_res = new Vector3d(opt_res.OptimalPoint[x], opt_res.OptimalPoint[y], opt_res.OptimalValue);
            Rhino.RhinoApp.WriteLine("OPTIMAL POINT : " + X_res);
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("{818A481B-5E60-45D7-8E2C-19DE2B40E96A}"); }
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddIntegerParameter("test_funtion_id", "test_funtion_id", "id of test function", GH_ParamAccess.item, test_function_id);
            pManager.AddNumberParameter("x_start","x_start","initial x coordinate",GH_ParamAccess.item, 1.00);
            pManager.AddNumberParameter("y_start","y_start","initial y coordinate",GH_ParamAccess.item, 1.00);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.Register_StringParam("info", "info", "solver info");
            pManager.Register_PointParam("Xk", "Xk", "sequence of Xk");
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Rhino.RhinoApp.WriteLine("DEBUG C");
            if (!DA.GetData(0, ref test_function_id)) { return; }
            if (!DA.GetData(1, ref x_start)) { return; }
            if (!DA.GetData(2, ref y_start)) { return; }
            Rhino.RhinoApp.WriteLine("DEBUG D");

            if (test_function_id != test_function_id_cache)
            {
                Rhino.RhinoApp.WriteLine("IPOPT COMPONENT EVENT : test_function_id has changed");

                F = TestFunction.GetTestFunction(test_function_id, X);
                opt.ObjectiveFunction = F;

                intermediateEventsList.Clear();
                opt_prepared = opt.Prepare();
                opt_res = opt_prepared.Run(new VariableAssignment[] { x0, y0 });

                X_res = new Vector3d(opt_res.OptimalPoint[x], opt_res.OptimalPoint[y], opt_res.OptimalValue);
                Rhino.RhinoApp.WriteLine("OPTIMAL POINT : " + X_res);

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

                test_function_id_cache = test_function_id;
            }
            else if(x_start != x_start_cache || y_start != y_start_cache)
            {
                Rhino.RhinoApp.WriteLine("IPOPT COMPONENT EVENT : x_start or y_start has changed");

                x0 = new VariableAssignment(x, x_start);
                y0 = new VariableAssignment(y, y_start);

                opt_res = opt_prepared.Run(new VariableAssignment[] { x0, y0 });

                X_res = new Vector3d(opt_res.OptimalPoint[x], opt_res.OptimalPoint[y], opt_res.OptimalValue);
                Rhino.RhinoApp.WriteLine("OPTIMAL POINT : " + X_res);

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

                x_start_cache = x_start;
                y_start_cache = y_start;
            }

            DA.SetData(0, info);
            DA.SetData(1, X_res);
        }

        public void IntermediateOutput(Object sender, IpoptIntermediateEventArgs e)
        {
            // store the event with current iteration informations
            intermediateEventsList.Add(e);
        }

    }
}
