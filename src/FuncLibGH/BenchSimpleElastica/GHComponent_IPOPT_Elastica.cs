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
    public class GHComponent_IPOPT_Elastica : GH_Component
    {
        // FIELDS
        private bool reset = true;
        private bool reset_cache = true;
        private Curve crv;
        private Polyline Xi, Xres;
        private int n;
        private double Lr, lr;
        private string info;

        IpoptOptimizer opt;
        PreparedOptimizer opt_prepared;
        IOptimizerResult opt_res;
        List<IpoptIntermediateEventArgs> intermediateEventsList;


        private bool normalized;

        public GHComponent_IPOPT_Elastica()
            : base("IPOPT over Test Function", "ELASTICA",
                "Show the IPOP solver on simple elastica",
                "OptiELA", "Test Functions")
        {
   
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("{95D43E08-AA49-4D61-A863-0C3DD96AC4B5}"); }
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Xi", "Xi", "initial position", GH_ParamAccess.item);
            pManager.AddIntegerParameter("n", "n", "number of elements", GH_ParamAccess.item);
            pManager.AddNumberParameter("Lr", "Lr", "rest length", GH_ParamAccess.item);
            pManager.AddBooleanParameter("reset", "reset", "restart the solver", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.Register_StringParam("info", "info", "solver info");
            pManager.Register_PointParam("Xf", "Xf", "results");
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Rhino.RhinoApp.WriteLine("DEBUG A");
            if (!DA.GetData(0, ref crv)) { return; }
            if (!DA.GetData(1, ref n)) { return; }
            if (!DA.GetData(2, ref Lr)) { return; }
            if (!DA.GetData(3, ref reset)) { return; }

            Rhino.RhinoApp.WriteLine("DEBUG B");
            if (!crv.TryGetPolyline(out Xi)) return;
            lr = Lr / n;

            for (int i = 0; i < Xi.Length; i++)
            {
                Rhino.RhinoApp.WriteLine("Xi["+i+"] = " + Xi[i]);
                
            }

            Rhino.RhinoApp.WriteLine("DEBUG C");
            VariableAssignment[] X0 = new VariableAssignment[2 * (n + 1)];
            Variable[] X = new Variable[2 * (n + 1)];
            Variable[] x = new Variable[n + 1];
            Variable[] y = new Variable[n + 1];
            for (int i = 0; i < n + 1; i++)
            {
                x[i] = new Variable("x_" + i);
                y[i] = new Variable("y_" + i);
                X[2 * i] = x[i];
                X[2 * i + 1] = y[i];
                X0[2 * i] = new VariableAssignment(x[i], Xi[i].X);
                X0[2 * i + 1] = new VariableAssignment(y[i], Xi[i].Y);
            }

            Rhino.RhinoApp.WriteLine("DEBUG D");
            // courbure au noeud i (portée par z) élevée au carré
            Function[] k2 = new Function[n - 1];

            for (int i = 0; i < n - 1; i++)
            {
                int j = i + 1;
                k2[i] = 4 * Function.Pow((x[j] - x[j - 1]) * (y[j + 1] - y[j]) - (y[j] - y[j - 1]) * (x[j + 1] - x[j]), 2);
                k2[i] = k2[i] / Function.Pow((Math.Pow(lr, 2) + (x[j] - x[j - 1]) * (x[j + 1] - x[j]) + (y[j] - y[j - 1]) * (y[j + 1] - y[j])), 2);
            }

            Rhino.RhinoApp.WriteLine("DEBUG E");
            // energy objective function
            double EI = 3.00;
            double alpha = EI/(2*lr);
            Function E = 0;
            for (int i = 0; i < n - 1; i++)
            {
                E += alpha * k2[i];
            }

            Rhino.RhinoApp.WriteLine("DEBUG F");
            // constraints
            double c = Math.Pow(lr, 2);
            FunctionConstraint[] C = new FunctionConstraint[n + 4];
            for (int i = 0; i < n; i++)
            {
                C[i] = new FunctionConstraint(Function.Pow(x[i + 1] - x[i], 2) + Function.Pow(y[i + 1] - y[i], 2), c, c);
            }
            C[n] = new FunctionConstraint(x[0], Xi[0].X, Xi[0].X);
            C[n + 1] = new FunctionConstraint(y[0], Xi[0].Y, Xi[0].Y);
            C[n + 2] = new FunctionConstraint(x[n], Xi[n].X, Xi[n].X);
            C[n + 3] = new FunctionConstraint(y[n], Xi[n].Y, Xi[n].Y);

            Rhino.RhinoApp.WriteLine("DEBUG G");
            opt = new IpoptOptimizer();
            intermediateEventsList = new List<IpoptIntermediateEventArgs>();
            opt.Intermediate += new EventHandler<IpoptIntermediateEventArgs>(IntermediateOutput);
            opt.PrintLevel = 5;
            opt.MaxIterations = 100;
            opt.Variables.Add(X);
            opt.Constraints.Add(C);
            opt.ObjectiveFunction = E;
            opt_prepared = opt.Prepare();

            Rhino.RhinoApp.WriteLine("DEBUG H");
            Rhino.RhinoApp.WriteLine(opt_prepared.ToString());
            opt_res = opt_prepared.Run(X0);

            Rhino.RhinoApp.WriteLine("DEBUG I");
            Point3d[] pts = new Point3d[n + 1];
            for (int i = 0; i < n + 1; i++)
            {
                pts[i].X = opt_res.OptimalPoint[x[i]];
                pts[i].Y = opt_res.OptimalPoint[y[i]];
                pts[i].Z = 0;
            }

            Rhino.RhinoApp.WriteLine("DEBUG J");
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

            DA.SetData(0, info);
            DA.SetDataList(1, pts);
        }

        public void IntermediateOutput(Object sender, IpoptIntermediateEventArgs e)
        {
            // store the event with current iteration informations
            intermediateEventsList.Add(e);
        }

        public double Curvature2(int j)
        {
            double k;
            k = 4 * Math.Pow((Xi[j].X - Xi[j - 1].X) * (Xi[j + 1].Y - Xi[j].Y) - (Xi[j].Y - Xi[j - 1].Y) * (Xi[j + 1].X - Xi[j].X), 2);
            k = k / Math.Pow((Math.Pow(lr, 2) + (Xi[j].X - Xi[j - 1].X) * (Xi[j + 1].X - Xi[j].X) + (Xi[j].Y - Xi[j - 1].Y) * (Xi[j + 1].Y - Xi[j].Y)), 2);
            return k;
        }       

    }
}
