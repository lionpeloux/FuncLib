using FuncLib.Functions;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FuncLibGH
{
    public enum TestFunctionType : int
    {
        Sphere = 1,
        Matyas = 2,
        McCormick = 3,
        Beale = 4,
        Goldestein = 5,
        Booth = 6
    }
    public static class TestFunction
    {
        public static Function GetTestFunction(int type, Variable[] X)
        {
            return GetTestFunction((TestFunctionType)type, X);
        }
        public static Function GetTestFunction(TestFunctionType type, Variable[] X)
        {
            switch (type)
            {
                case TestFunctionType.Sphere:
                    return F_Sphere(X);

                case TestFunctionType.Matyas:
                    return F_Matyas(X);

                case TestFunctionType.McCormick:
                    return F_McCormick(X);

                default:
                    return F_Sphere(X);
            }
        }

        // Sphere = 1
        public static Function F_Sphere(Variable[] X)
        {
            Function F = Function.Pow(X[0], 2); ;
            for (int i = 1; i < X.Length; i++)
            {
                F += Function.Pow(X[i], 2);
            }
            return F;
        }

        // Matyas = 2
        public static Function F_Matyas(Variable[] X)
        {
            Function F = 0.26 * (Function.Pow(X[0], 2) + Function.Pow(X[1], 2)) - 0.48 * X[0] * X[1];
            return F;
        }
        // McCormick = 3
        public static Function F_McCormick(Variable[] X)
        {
            Function F = Function.Sin(X[0] + X[1]) + Function.Pow(X[0] - X[1], 2) - 1.5 * X[0] + 2.5 * X[1] + 1;
            return F;
        }
    }

}
