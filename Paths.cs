/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace m_435_NCAMR
{
    static class Paths
    {
        private static string dist_repo = @"C:\Users\Will\Desktop\dirs\school\class-dirs\repo\distributions";
        private static string img_repo = @"C:\Users\Will\Desktop\dirs\school\class-dirs\repo\images";
        private static string vb_repo = @"C:\Users\Will\Desktop\dirs\school\class-dirs\repo\visual-basic-templates";
        private static string mat_repo = @"C:\Users\Will\Desktop\dirs\school\class-dirs\repo\matrices";
        private static string nav_stokes_solutions = @"C:\Users\Will\Desktop\dirs\school\class-dirs\repo\n-s-solutions";
        public static string NavSolutionsRepo
        {
            get { return nav_stokes_solutions; }
        }
        public static string DistributionRepo
        {
            get { return dist_repo; }
        }
        public static string ImageRepo
        {
            get { return img_repo; }
        }
        public static string MatrixRepo
        {
            get { return mat_repo; }
        }
    }
}
