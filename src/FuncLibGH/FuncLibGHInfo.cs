using System;
using System.Drawing;
using Grasshopper.Kernel;

namespace FuncLibGH
{
    public class FuncLibGHInfo : GH_AssemblyInfo
    {
        public override string Name
        {
            get
            {
                return "FuncLibGH";
            }
        }
        public override Bitmap Icon
        {
            get
            {
                //Return a 24x24 pixel bitmap to represent this GHA library.
                return null;
            }
        }
        public override string Description
        {
            get
            {
                //Return a short string describing the purpose of this GHA library.
                return "";
            }
        }
        public override Guid Id
        {
            get
            {
                return new Guid("205be62a-fadc-4bb2-bac7-56797a5406bb");
            }
        }

        public override string AuthorName
        {
            get
            {
                //Return a string identifying you or your company.
                return "";
            }
        }
        public override string AuthorContact
        {
            get
            {
                //Return a string representing your preferred contact details.
                return "";
            }
        }
    }
}
