// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Collections.Generic;
using FuncLib.Functions;
using FuncLib.Functions.Compilation;
using FuncLib.Optimization;
using FuncLib.Optimization.Ipopt;
using FuncLib.Utilities;

namespace FuncLibConsole.Samples
{
	public class Optimizers
	{
        public static List<IpoptIntermediateEventArgs> intermediateEventsList;
		
        public static void Sample1()
		{
			// Variables. The string names are added for convenience when printing out the optimal point with Evaluation.ToString().
			Variable x = new Variable("x");
			Variable y = new Variable("y");

			// Rosenbrock function http://en.wikipedia.org/wiki/Rosenbrock_function.
			Function f = Function.Sqr(1.0 - x) + 100.0 * Function.Sqr(y - Function.Sqr(x));

			// Use the BFGS optimizer.
			BfgsOptimizer o = new BfgsOptimizer();

			// Specify variables and objective functions and add constraints. Derivatives are computed automatically.
			o.Variables.Add(x, y);
			o.VariableConstraints.Add(x >= 0.0, y >= 0.0);
			o.ObjectiveFunction = f;

			// Start optimizer from a random point.
			Random r = new Random(1);
			IOptimizerResult or = o.Run(x | r.NextDouble(), y | r.NextDouble());

			// Show convergence status and optimal point and value.
			Console.WriteLine(or.Status);
			Console.WriteLine(or.OptimalPoint);
			Console.WriteLine(or.OptimalValue);
		}

		public static void Sample2()
		{
			// See IpoptHowTo.txt how to get Ipopt to work. It's quite easy and Ipopt is very powerful.

			// Variables. The string names are added for convenience when printing out the optimal point with Evaluation.ToString().
			Variable x = new Variable("x");
			Variable y = new Variable("y");

			// Rosenbrock function http://en.wikipedia.org/wiki/Rosenbrock_function.
			Function f = Function.Sqr(1.0 - x) + 100.0 * Function.Sqr(y - Function.Sqr(x));

			// Use the BFGS optimizer.
			IpoptOptimizer o = new IpoptOptimizer();

			// Specify variables and objective functions and add constraints. Derivatives are computed automatically.
			o.Variables.Add(x, y);
			o.Constraints.Add(x >= 0.0, y >= 0.0);
			o.ObjectiveFunction = f;

			// Prepare the optimizer.
			PreparedOptimizer po = o.Prepare();

			// Run the prepared optimizer from many different points.
			Random r = new Random(1);
			for (int i = 0; i < 100; i++)
			{
				IOptimizerResult or = po.Run(x | r.NextDouble(), y | r.NextDouble());
				Console.WriteLine(or.OptimalPoint);
			}
		}

		public static void Sample3()
		{
			// See IpoptHowTo.txt how to get Ipopt to work. It's quite easy and Ipopt is very powerful.

			// Variables. The string names are added for convenience when printing out the optimal point with Evaluation.ToString().
			Variable x1 = new Variable("x1");
			Variable x2 = new Variable("x2");
			Variable x3 = new Variable("x3");
			Variable x4 = new Variable("x4");

			// Objective function and non-linear constraints.
			Function f = x1 * x4 * (x1 + x2 + x3) + x3;
			Function g1 = x1 * x2 * x3 * x4;
			Function g2 = Function.Sqr(x1) + Function.Sqr(x2) + Function.Sqr(x3) + Function.Sqr(x4);

			//CompiledFunction[] c = Compiler.Compile(new Function[] { f, g1, g2 }, new Variable[] { x1, x2, x3, x4 }, 2);
			//f = c[0];
			//g1 = c[1];
			//g2 = c[2];

			// Prepare the optimizer.
			IpoptOptimizer o = new IpoptOptimizer();
			o.Variables.Add(x1, x2, x3, x4);
			o.ObjectiveFunction = f;
			o.Constraints.Add(g1 >= 25.0);
			o.Constraints.Add(g2 == 40.0);
			o.Constraints.Add(1.0 <= x1, x1 <= 5.0);
			o.Constraints.Add(1.0 <= x2, x2 <= 5.0);
			o.Constraints.Add(1.0 <= x3, x3 <= 5.0);
			o.Constraints.Add(1.0 <= x4, x4 <= 5.0);

			// Verbose mode. Show Ipopt convergence.
			o.PrintLevel = 5;

			// Run optimization starting from (x1, x2, x3, x4) = (1, 5, 5, 1).
			IOptimizerResult or = o.Run(x1 | 1.0, x2 | 5.0, x3 | 5.0, x4 | 1.0);

			Console.WriteLine("x1 = " + or.OptimalPoint[x1]);
			Console.WriteLine("x2 = " + or.OptimalPoint[x2]);
			Console.WriteLine("x3 = " + or.OptimalPoint[x3]);
			Console.WriteLine("x4 = " + or.OptimalPoint[x4]);
			Console.WriteLine("f = " + f.Value(or.OptimalPoint));
			Console.WriteLine("g1 = " + g1.Value(or.OptimalPoint));
			Console.WriteLine("g2 = " + g2.Value(or.OptimalPoint));
		}

        public static void Sample_LDP_1()
        {
            
            Console.WindowWidth = Console.LargestWindowWidth-20;
            Console.WindowHeight = Console.LargestWindowHeight-20;
            Console.SetWindowPosition(0, 0);

            // EXEMPLE : optimization sans contraintes sur la sphère de R2 

            // Variables. The string names are added for convenience when printing out the optimal point with Evaluation.ToString().
            int n = 2;
            Variable[] vec = new Variable[n];
            for (int i = 0; i < n; i++)
			{
                vec[i] = new Variable("x_" + i);
			}

            Function f = Function.Pow(vec[0], 2);
            for (int i = 1; i < n; i++)
            {
                f += Function.Pow(vec[i], 2);
            }

            // Use the BFGS optimizer.
            IpoptOptimizer o = new IpoptOptimizer();

            // Specify variables and objective functions and add constraints. Derivatives are computed automatically.
            for (int i = 0; i < n; i++)
            {
                o.Variables.Add(vec);
            }

            intermediateEventsList = new List<IpoptIntermediateEventArgs>();
            o.Constraints.Add(vec[0] >=  vec[1]);
            o.ObjectiveFunction = f;
            o.Intermediate += new EventHandler<IpoptIntermediateEventArgs>(IntermediateOutput);

            // Prepare the optimizer.

            o.PrintLevel = 6;
            o.MaxIterations = 2;
            PreparedOptimizer po = o.Prepare();
            
           
            // Run the prepared optimizer from many different points.
            VariableAssignment[] startpoint = new VariableAssignment[n];
            for (int j = 0; j < n; j++)
            {
                startpoint[j] = new VariableAssignment(vec[j], 1);
            }
            
            IOptimizerResult or = po.Run(startpoint);

            var table = intermediateEventsList.ToStringTable(
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

            
            Console.WriteLine(table);

            Console.WriteLine(or.OptimalPoint);
            Console.WriteLine(or.Status);
            
        }

        public static void IntermediateOutput(Object sender, IpoptIntermediateEventArgs e)
        {
            // store the event with current iteration informations
            intermediateEventsList.Add(e);
        }
	}
}
