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
    public class GHComponent_BenchSimpleElastica_FuncLibBFGS : GH_Component
    {
        // FIELDS
        private bool reset = true;
        private bool reset_cache = true;
        private Curve Xi_crv;
        private Curve X_bar_crv;
        private Polyline Xi, X_bar;

        private SimpleElastica_FuncLib_BFGS simpleElastica_funclib_bfgs;


        public GHComponent_BenchSimpleElastica_FuncLibBFGS()
            : base("FuncLib BFGS","FuncLib BFGS","Simple Elastica using FuncLib BFGS solver with equality constraint for suppirts", "OptiELA", "Simple Elastica Bench")
        {
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("{991A228C-6890-4C18-9336-B00D1B3D2815}"); }
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Xi", "Xi", "initial position", GH_ParamAccess.item);
            pManager.AddCurveParameter("X_bar", "X_bar", "initial position", GH_ParamAccess.item);
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
            if (!DA.GetData(2, ref reset)) { return; }

            if (!Xi_crv.TryGetPolyline(out Xi)) return;
            if (!X_bar_crv.TryGetPolyline(out X_bar)) return;

            if (reset != reset_cache)
            {
                Rhino.RhinoApp.WriteLine("DEBUG B");
                simpleElastica_funclib_bfgs.Run();
                reset_cache = reset;
            }
            else
            {
                Rhino.RhinoApp.WriteLine("DEBUG A");

                simpleElastica_funclib_bfgs = new SimpleElastica_FuncLib_BFGS(X_bar, Xi);
                simpleElastica_funclib_bfgs.opt.EpsG = 1e-10;
                simpleElastica_funclib_bfgs.Run();
            }
            
            DA.SetData(0, simpleElastica_funclib_bfgs.info);
            DA.SetDataList(1, simpleElastica_funclib_bfgs.X);
        }


    }
}
