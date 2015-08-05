// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.CodeDom.Compiler;

namespace FuncLib.Functions.Compilation
{
	[Serializable]
	public class CompilerException : Exception
	{
		private CompilerErrorCollection errors;

		public CompilerException(string message, CompilerErrorCollection errors)
			: base(message)
		{
			this.errors = errors;
		}

		public CompilerErrorCollection Errors
		{
			get
			{
				return errors;
			}
		}
	}
}
