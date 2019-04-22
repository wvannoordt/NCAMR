/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;
using System.IO;
using System.IO.Compression;
using _3DSimple;
using System.Diagnostics;

namespace m_435_NCAMR
{
    class Program
    {
        static void Main(string[] args)
        {
            NCGrid_Distribution funcp = new NCGrid_Distribution(new RBounds2D(0, 10, 0, 10), 400, 400);
            NCGrid_Distribution func = ExactFunctionGeneratorVB2D.GenerateFunctionToGrid("(-1*x) + (2*y)", funcp);
            DistributionSketchSettings settings = DistributionSketchSettings.Fancy();
            int n = 10;
            settings.XLabelCount = n;
            settings.YLabelCount = n;
            settings.XGridlineCount = n;
            settings.YGridlineCount = n;
            DistributionSketch2D sktch = new DistributionSketch2D(func, settings);
            sktch.CreateSketch(true);
            sktch.SaveImage(@"C:\Users\Will\Desktop\output.png", true);
        }
        static void run_presentation_transient(int numframes)
        {
            double t = 0;
            double dt = 0.1;
            int n = 25;
            double L = 10;
            double H = 10;
            RBounds2D domain = new RBounds2D(0, L, 0, H);
            NCGrid_Distribution soln = new NCGrid_Distribution(domain, n, n);
            double max1 = 4;
            double max2 = 3;
            double max3 = 5;
            double max4 = 2;
            double omega1 = 1;
            double omega2 = 0.6;
            double omega3 = 0.2;
            double omega4 = 1.2;
            double phi1 = -0.4;
            double phi2 = -0.9;
            double phi3 = 1.3;
            double phi4 = 0.2;
            double dx = L / (n - 1);
            double dy = H / (n - 1);
            double refinestep = 0.1;
            int refinecount = 3;

            for (int i = 0; i < numframes; i++)
            {
                BoundaryConditions cond = new BoundaryConditions(soln, BoundaryConditions.BoundaryConditionType.Dirichlet);
                t += dt;
                for (int q = 0; q < n; q++)
                {
                    double x = q * dx;
                    double y = q * dy;
                    cond[q, BoundaryConditions.Direction.Negative_X] = max1 * Math.Sin(4*Math.PI * x / L) * Math.Sin(omega1 * t - phi1);
                    cond[q, BoundaryConditions.Direction.Negative_Y] = max2 * Math.Sin(2*Math.PI * y / L) * Math.Sin(omega2 * t - phi2);
                    cond[q, BoundaryConditions.Direction.Positive_X] = max3 * Math.Sin(3*Math.PI * x / L) * Math.Sin(omega3 * t - phi3);
                    cond[q, BoundaryConditions.Direction.Positive_Y] = max4 * Math.Sin(5*Math.PI * y / L) * Math.Sin(omega4 * t - phi4);
                }
                BVPLinear2D problem = new BVPLinear2D(cond, LinearOperatorOrder2.Laplace, soln);
                NCGrid_Distribution sol = problem.SolveKaczMarzExt(100, 10, 20);
                sol.ApplyMeshMorphGA(refinestep);
                for (int z = 0; z < refinecount-1; z++)
                {
                    problem = new BVPLinear2D(cond, LinearOperatorOrder2.Laplace, soln);
                    sol = problem.Solve(Matrix.SystemSolvingScheme.Kaczmarz);
                    sol.ApplyMeshMorphGA(refinestep);
                }
                sol.WriteToFile("longtransient_" + i.ToString());
                sol.QuickSketch("sol_" + bufferint(i, 4));
            }
        }
        static void testfunction()
        {
            Console.WriteLine("Press enter to begin, or enter \"c\" to clear repos.");
            if (Console.ReadLine() == "c")
            {
                RepoManagement.ClearRepo(Paths.DistributionRepo, Paths.ImageRepo);
            }
            DateTime then = DateTime.Now;
            Console.WriteLine("Solving...");
            double L = 10;
            double H = 10;
            int n = 38;
            RBounds2D bounds = new RBounds2D(0, L, 0, H);
            NCGrid_Distribution dist = new NCGrid_Distribution(bounds, n, n);
            BoundaryConditions conditions = new BoundaryConditions(bounds, n, n, BoundaryConditions.BoundaryConditionType.Dirichlet);
            double dx = L / (n - 1);
            double max = 5;
            string fxn = string.Format("{0}*Sin({1}*x/{2})*(Exp({1}*y/{2}) - Exp(-1*{1}*y/{2}))/(Exp({1}*{3}/{2}) - Exp(-1*{1}*{3}/{2}))", max, Math.PI, L, H);
            conditions.SetConstant(0, BoundaryConditions.Direction.Negative_X);
            conditions.SetConstant(0, BoundaryConditions.Direction.Negative_Y);
            conditions.SetConstant(0, BoundaryConditions.Direction.Positive_Y);
            for (int i = 0; i < n; i++)
            {
                double x = i * dx;
                double z = max * Math.Sin(Math.PI * x / L);
                conditions[i, BoundaryConditions.Direction.Positive_Y] = z;
            }
            int solcount = 50;
            double[] errors = new double[solcount];
            double[] iteration = new double[solcount];
            DistributionSketchSettings S = DistributionSketchSettings.Fancy();
            S.SetFigureTitle("Temperature Distribution");
            for (int i = 0; i < solcount; i++)
            {
                Console.WriteLine(i.ToString() + " of " + solcount.ToString() + " iterations processed.");
                iteration[i] = i;
                BVPLinear2D problem = new BVPLinear2D(conditions, LinearOperatorOrder2.Laplace, dist);
                problem.EnableConsoleOutput();
                NCGrid_Distribution soln = problem.Solve(Matrix.SystemSolvingScheme.Kaczmarz);
                NCGrid_Distribution analytic = ExactFunctionGeneratorVB2D.GenerateFunctionToGrid(fxn, soln);
                soln.ApplyMeshMorphGA(15);
                errors[i] = ErrorEstimation.NormDifference(soln, analytic);
                string title = "iterative-" + i.ToString();
                soln.WriteToFile(title);
                DistributionSketch2D sketch = new DistributionSketch2D(soln, S);
                dist = soln.Clone();
                sketch.CreateSketch(true);
                sketch.SaveImage(title + "-plot", false);
            }
            List<string> filestuff = new List<string>();
            for (int i = 0; i < iteration.Length; i++)
            {
                filestuff.Add(iteration[i].ToString() + "," + errors[i].ToString());
            }
            File.WriteAllLines(@"C:\Users\Will\Desktop\Folders\MATH435\repo\curves-1d\errors-temp.csv", filestuff.ToArray());
            Console.WriteLine("Done.");
            Console.ReadLine();
        }
        static void testfunction2()
        {
            Console.WriteLine("Press enter to begin, or enter \"c\" to clear repos.");
            if (Console.ReadLine() == "c")
            {
                RepoManagement.ClearRepo(Paths.DistributionRepo, Paths.ImageRepo);
            }
            DateTime then = DateTime.Now;
            Console.WriteLine("Solving...");
            double L = 10;
            double H = 10;
            int n = 38;
            RBounds2D bounds = new RBounds2D(0, L, 0, H);
            NCGrid_Distribution dist = new NCGrid_Distribution(bounds, n, n);
            BoundaryConditions conditions = new BoundaryConditions(bounds, n, n, BoundaryConditions.BoundaryConditionType.Dirichlet);
            double dx = L / (n - 1);
            double omega = 4;
            double max = 8;
            string fxn = string.Format("{4}*Sin({0}*{1}*x/{2})*Exp(y/{3})", Math.PI,omega,L,H, max);
            ExactFunctionGeneratorVB2D.quickPlot(fxn, "analytic");
            conditions.SetConstant(0, BoundaryConditions.Direction.Negative_X);
            conditions.SetConstant(0, BoundaryConditions.Direction.Positive_X);
            for (int i = 0; i < n; i++)
            {
                double x = i * dx;
                double z_neg = max * Math.Sin(Math.PI * omega*x / L);
                double z_pos = max * Math.E*Math.Sin(Math.PI * omega*x / L);
                conditions[i, BoundaryConditions.Direction.Positive_Y] = z_pos;
                conditions[i, BoundaryConditions.Direction.Negative_Y] = z_neg;
            }
            int solcount = 15;
            LinearOperatorOrder2 op = new LinearOperatorOrder2(0, 0, 0, 1 / (H * H), Math.PI * Math.PI * omega * omega / (L * L), 0);
            double[] errors = new double[solcount];
            double[] iteration = new double[solcount];
            DistributionSketchSettings S = DistributionSketchSettings.Fancy();
            S.SetFigureTitle("Temperature Distribution");
            for (int i = 0; i < solcount; i++)
            {
                Console.WriteLine(i.ToString() + " of " + solcount.ToString() + " iterations processed.");
                iteration[i] = i;
                BVPLinear2D problem = new BVPLinear2D(conditions, op, dist);
                problem.EnableConsoleOutput();
                NCGrid_Distribution soln = problem.SolveSRDD();
                NCGrid_Distribution analytic = ExactFunctionGeneratorVB2D.GenerateFunctionToGrid(fxn, soln);
                soln.ApplyMeshMorphGA(15, 0.0000015);
                errors[i] = ErrorEstimation.NormDifference(soln, analytic);
                string title = "iterative-" + i.ToString();
                soln.WriteToFile(title);
                DistributionSketch2D sketch = new DistributionSketch2D(soln, S);
                dist = soln.Clone();
                sketch.CreateSketch(true);
                sketch.SaveImage(title + "-plot", false);
            }
            List<string> filestuff = new List<string>();
            for (int i = 0; i < iteration.Length; i++)
            {
                filestuff.Add(iteration[i].ToString() + "," + errors[i].ToString());
            }
            File.WriteAllLines(@"C:\Users\Will\Desktop\Folders\MATH435\repo\curves-1d\errors-temp.csv", filestuff.ToArray());
            Console.WriteLine("Done.");
            Console.ReadLine();
        }
        static void testfunction3()
        {
            Console.WriteLine("Press enter to begin, or enter \"c\" to clear repos.");
            if (Console.ReadLine() == "c")
            {
                RepoManagement.ClearRepo(Paths.DistributionRepo, Paths.ImageRepo);
            }
            DateTime then = DateTime.Now;
            Console.WriteLine("Solving...");
            double L = 10;
            double H = 10;
            int n = 26;
            RBounds2D bounds = new RBounds2D(0, L, 0, H);
            NCGrid_Distribution dist = new NCGrid_Distribution(bounds, n, n);
            BoundaryConditions conditions = new BoundaryConditions(bounds, n, n, BoundaryConditions.BoundaryConditionType.Dirichlet);
            double dx = L / (n - 1);
            double omega = 4;
            double max = 8;
            conditions.SetConstant(0, BoundaryConditions.Direction.Negative_X);
            conditions.SetConstant(0, BoundaryConditions.Direction.Negative_Y);
            for (int i = 0; i < n; i++)
            {
                double x = i * dx;
                double y = i * dx;
                double ybound = max / (1 + Math.Exp(-1 * (omega * x / L)));
                double xbound = max / (1 + Math.Exp(-1 * (omega * y / H)));
                conditions[i, BoundaryConditions.Direction.Positive_Y] = ybound;
                conditions[i, BoundaryConditions.Direction.Positive_X] = xbound;
            }
            int solcount = 15;
            LinearOperatorOrder2 op = LinearOperatorOrder2.Laplace;
            DistributionSketchSettings S = DistributionSketchSettings.Fancy();
            S.SetFigureTitle("Double Logistic Boundary");
            for (int i = 0; i < solcount; i++)
            {
                Console.WriteLine(i.ToString() + " of " + solcount.ToString() + " iterations processed.");
                BVPLinear2D problem = new BVPLinear2D(conditions, op, dist);
                problem.EnableConsoleOutput();
                NCGrid_Distribution soln = problem.SolveSRDD();
                string title = "iterative-" + i.ToString();
                soln.WriteToFile(title);
                DistributionSketch2D sketch = new DistributionSketch2D(soln, S);
                dist = soln.Clone();
                sketch.CreateSketch(true);
                sketch.SaveImage(title + "-plot", false);
                //dist.ApplyMeshMorphGA(2, 0.0019);
                Random R = new Random();
                for (int h = 1; h < dist.Xcount-1; h++)
                {
                    for (int k = 1; k < dist.Ycount-1; k++)
                    {
                        double ddx = (0.5-R.NextDouble()) * dx*0.6;
                        double ddy = (0.5-R.NextDouble()) * dx*0.6;
                        Vector3 move = new Vector3(ddx, ddy, 0);
                        dist[h, k] = dist[h, k] + move;
                    }
                }
            }
            Console.WriteLine("Done.");
            Console.ReadLine();
        }
        public static void printMatrix(Matrix M)
        {
            int shelfsize = 8;
            int trunc = 6;
            for (int i = 0; i < M.Rows; i++)
            {
                string line = string.Empty;
                for (int j = 0; j < M.Columns; j++)
                {
                    string entry = M[i, j].ToString();
                    string tr_entry = entry.Substring(0, Math.Min(trunc, entry.Length));
                    line += tr_entry.PadRight(shelfsize);
                }
                Console.WriteLine(line);
            }
        }
        static string bufferint(int n, int size)
        {
            string l = string.Empty;
            for (int i = 0; i < size + 1 - n.ToString().Length; i++)
            {
                l += "0";
            }
            l += n.ToString();
            return l;
        }
    }
}
