// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;

using FuncLibConsole.Samples;

namespace FuncLibConsole
{
	public class Program
	{
		public static void Main(string[] args)
		{
			//HelloWorld.Sample1();
			//HelloWorld.Sample2();
			//MatrixAlgebra.Sample2();
			Optimizers.Sample1();
			//Optimizers.Sample2();
			//Optimizers.Sample3();
			//Dual.Sample1();
			//Dual.Sample2();
			//Dual.Sample3();
            //Optimizers.Sample_LDP_1();
			Console.WriteLine("Done");
			Console.ReadKey();
		}
	}
}
