using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using FuncLib.Functions;
using FuncLib.Optimization.Ipopt;
using FuncLib.Optimization;

namespace FuncLibGH
{
    public class GHComponent_Plot_TestFunction : GH_Component
    {
        // FIELDS
        private int test_function_id, test_function_id_cache;
        Interval x_domain, x_domain_cache, y_domain, y_domain_cache;
        private int nx, nx_cache;
        private int ny, ny_cache;
        private List<Point3d> pts_grid;
        private double z_min, z_max;
        double xi, yj, zij;
        private bool normalized;

        Variable[] X;
        Variable x, y;
        Function F;

        public GHComponent_Plot_TestFunction()
            : base("Plot Test Function", "Plot",
                "Plot a given test function",
                "OptiELA", "Test Functions")
        {
            x = new Variable("x");
            y = new Variable("y");
            X = new Variable[] { x, y };
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("{58039541-D1DE-489F-91D7-5FFAA1CFD000}"); }
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddIntegerParameter("test_funtion_id", "test_funtion_id", "id of test function", GH_ParamAccess.item, test_function_id);
            pManager.AddIntervalParameter("x_domain", "x_domain", "display F(x,y) over the x domain", GH_ParamAccess.item, new Interval(-1.00, 1.00));
            pManager.AddIntervalParameter("y_domain", "y_domain", "display F(x,y) over the y domain", GH_ParamAccess.item, new Interval(-1.00, 1.00));
            pManager.AddIntegerParameter("nx", "nx", "number of x_domain subdivisions", GH_ParamAccess.item, 20);
            pManager.AddIntegerParameter("ny", "ny", "number of y_domain subdivisions", GH_ParamAccess.item, 20);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.Register_StringParam("info", "info", "solver info");
            pManager.Register_PointParam("grid_pts", "grid_pts", "gridpoints plot representing z = f(x,y) surface");
            pManager.Register_DoubleParam("z_min", "z_min", "minimum of zij = f(xi,yj) over the gridpoints plot");
            pManager.Register_DoubleParam("z_max", "z_max", "maximum of zij = f(xi,yj) over the gridpoints plot");
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            if (!DA.GetData(0, ref test_function_id)) { return; }
            if (!DA.GetData(1, ref x_domain)) { return; }
            if (!DA.GetData(2, ref y_domain)) { return; }
            if (!DA.GetData(3, ref nx)) { return; }
            if (!DA.GetData(4, ref ny)) { return; }

            if (test_function_id != test_function_id_cache)
            {
                F = TestFunction.GetTestFunction(test_function_id, X);
                test_function_id_cache = test_function_id;
                RefreshPointsGrid();
            }
            else if (x_domain != x_domain_cache || y_domain != y_domain_cache)
            {
                x_domain_cache = x_domain;
                y_domain_cache = y_domain;
                RefreshPointsGrid();
            }
            else if (nx != nx_cache || ny != ny_cache)
            {
                nx_cache = nx;
                ny_cache = ny;
                RefreshPointsGrid();
            }

            DA.SetDataList(1, pts_grid);
            DA.SetData(2, z_min);
            DA.SetData(3, z_max);
        }

        private void RefreshPointsGrid()
        {
            pts_grid = new List<Point3d>(nx * ny);
            z_min = F.Value(new VariableAssignment[] { new VariableAssignment(x, x_domain.Min), new VariableAssignment(y, y_domain.Min) });
            z_max = F.Value(new VariableAssignment[] { new VariableAssignment(x, x_domain.Min), new VariableAssignment(y, y_domain.Min) });
            for (int i = 0; i < nx; i++)
            {
                xi = x_domain.Min + ((double)i / (nx - 1)) * x_domain.Length;
                for (int j = 0; j < ny; j++)
                {
                    yj = y_domain.Min + ((double)j / (ny - 1)) * y_domain.Length;
                    zij = F.Value(new VariableAssignment[] { new VariableAssignment(x, xi), new VariableAssignment(y, yj) });
                    if (zij < z_min) z_min = zij;
                    if (zij > z_max) z_max = zij;
                    pts_grid.Add(new Point3d(xi, yj, zij));
                }
            }
            if(normalized)
            {
                double zij_bar;
                double[] I = new double[]{z_min, z_max};
                double[] I_bar = new double[]{0, 1};
                for (int i = 0; i < nx; i++)
                {
                    for (int j = 0; j < ny; j++)
                    {
                        zij_bar = RemapNumber(pts_grid[nx * i + j].Z, I, I_bar);
                        pts_grid[nx * i + j] = new Point3d(pts_grid[nx * i + j].X, pts_grid[nx * i + j].Y, zij_bar);
                    }
                }
            }
        }
        private double RemapNumber(double z, double[] I, double[] I_bar)
        {
            double l = I[1] - I[0];
            double l_bar = I_bar[1] - I_bar[0];
            double z_bar = z * l_bar / l;
            return z_bar;
        }
    }
}
