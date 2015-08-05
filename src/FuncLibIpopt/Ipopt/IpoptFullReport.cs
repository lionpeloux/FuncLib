using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FuncLib.Optimization.Ipopt
{
    public class IpoptFullReport
    {
        public int _n;
        public int _m;
        public int _nele_jac;
        public int _nele_hess;
        public double[] _x_L;
        public double[] _x_U;
        public double[] _g_L;
        public double[] _g_U;

        public List<IpoptIntermediateEventArgs> intermediate_results;
        public IpoptOptimizerResult final_result;

        public IpoptFullReport(int _n, int _m, int _nele_jac, int _nele_hess, double[] _x_L, double[] _x_U,double[] _g_L,double[] _g_U)
        {
            this._n = _n;
            this._m = _m;
            this._nele_jac = _nele_jac;
            this._nele_hess = _nele_hess;
            this._x_L = _x_L;
            this._x_U = _x_U;
            this._g_L = _g_L;
            this._g_U = _g_U;

            this.intermediate_results = new List<IpoptIntermediateEventArgs>();
        }
    }
}
