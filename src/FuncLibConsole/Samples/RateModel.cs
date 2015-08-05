using System;

using FuncLib.Functions;
using FuncLib.Functions.Compilation;
using FuncLib.Optimization.Ipopt;

namespace FuncLibConsole.Samples
{
	public class RateModel
	{
		private Variable kappa, theta, sigma, r0;

		public void Run()
		{
			kappa = new Variable();
			theta = new Variable();
			sigma = new Variable();
			r0 = new Variable();

			// Use the objective function and a single non-linear constraint defined in the following methods.
			Function f = ObjectiveFunction();
			Function g = Constraint();
			
			// Replace the long and slow objective function with it's compiled equivalent (including up to second order partial derivatives). This sample
			// objective function is not really slow.
			f = Compiler.Compile(f, new Variable[] { kappa, theta, sigma, r0 }, 2);

			IpoptOptimizer o = new IpoptOptimizer();
			o.MaxIterations = 10;
			o.Variables.Add(kappa, theta, sigma, r0);
			o.ObjectiveFunction = f;
			//o.VariableEqualityConstraints
			o.VariableConstraints.Add(kappa >= 0.0, theta >= 0.0, sigma == 0.02, r0 >= 0.0);
			o.Constraints.Add(g >= 0.15, g <= 0.2);
			o.Run(kappa | 0.02);
		}

		public Function Rate(double maturity)
		{
			// Formula from Björk (xx.xx).
			return 0;
		}

		public Function Volatility(double maturity)
		{
			return 0;
		}

		public Function ObjectiveFunction()
		{
			//int n = 2;
			//double[] t = new double[] { 0.5, 5.0 };
			//double[] r = new double[] { 0.04, 0.05 };

			return Function.Sqrt(Function.Sqr(Rate(0.5) - 0.04) + Function.Sqr(Rate(5.0) - 0.04));
		}

		public Function Constraint()
		{
			return Volatility(10.0);
		}

		public void TestTimeWithDerivatives(Function f)
		{
			Function f_kappa = f.Derivative(kappa);
			Function f_theta = f.Derivative(theta);

			// Run a few times before starting the watch.
		}

		private class DaiSingleton
		{
			private Variable kappa11, kappa12, kappa21, kappa22, kappa33;
			private Variable theta1, theta2;
			private Variable sigma31, sigma32;
			private Variable alpha3;
			private Variable beta11, beta22, beta32;
			private Variable delta0, delta1;
			private Variable y1, y2, y3;

			public void Run()
			{
				IpoptOptimizer o = new IpoptOptimizer();
				o.Variables.Add(kappa11, kappa12, kappa21, kappa22, kappa33);
				o.Variables.Add(theta1, theta2);
				o.Variables.Add(sigma31, sigma32);
				o.Variables.Add(alpha3);
				o.Variables.Add(beta11, beta22, beta32);
				o.Variables.Add(delta0, delta1);
				o.Variables.Add(y1, y2, y3);
				o.ObjectiveFunction = Compiler.Compile(ObjectiveFunction(), o.Variables.ToArray(), 2);
				o.VariableConstraints.Add(kappa12 <= 0.0, kappa21 <= 0.0, alpha3 >= 0.0, beta11 >= 0.0, beta22 >= 0.0, beta32 >= 0.0, y1 >= 0.0, y2 >= 0.0);
				o.Constraints.Add(kappa11 * theta1 + kappa12 * theta2 >= 0.0);
				o.Constraints.Add(kappa21 * theta1 + kappa22 * theta2 >= 0.0);
				o.Constraints.Add(Rate(30.0) >= 0.03);
				o.Run(kappa11 | 0.1, kappa12 | 0.0);
			}

			private Function ObjectiveFunction()
			{
				return null;
			}

			private Function Rate(double maturity)
			{
				return null;
			}
		}
	}
}
