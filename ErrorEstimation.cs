/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using _3DSimple;


namespace m_435_NCAMR
{
    static class ErrorEstimation
    {
        public static NCGrid_Distribution CompareNiceClean(NCGrid_Distribution test, string vbs_test_function)
        {
            string title = "temp";
            string path = Paths.DistributionRepo + "\\" + title + ".dist";
            ExactFunctionGeneratorVB2D.GenerateFunction(vbs_test_function, test, title);
            NCGrid_Distribution true_dist = NCGrid_Distribution.from_file(title, false);
            File.Delete(path);
            NCGrid_Distribution dif = test - true_dist;
            dif.force_extrema_update();
            return dif;
        }
        public static double NormDifference(NCGrid_Distribution A, NCGrid_Distribution B)
        {
            NCGrid_Distribution newgrid = new NCGrid_Distribution(A.Bounds, A.Xcount+1, A.Ycount+1);
            int n = newgrid.Ycount;
            int m = newgrid.Xcount;
            double xmin = newgrid.Xmin;
            double xmax = newgrid.Xmax;
            double ymin = newgrid.Ymin;
            double ymax = newgrid.Ymax;
            NCAMRNode ll = new NCAMRNode(xmin, ymin, 0);
            NCAMRNode lr = new NCAMRNode(xmax, ymin, 0);
            NCAMRNode ul = new NCAMRNode(xmin, ymax, 0);
            NCAMRNode ur = new NCAMRNode(xmax, ymax, 0);
            newgrid[0, 0] = ll;
            newgrid[0, n-1] = ul;
            newgrid[m-1, 0] = lr;
            newgrid[m-1, n-1] = ur;
            for (int i = 0; i < A.Xcount-1; i++)
            {
                newgrid[i + 1, 0] = xyavg(A[i, 0], A[i + 1, 0]);
                newgrid[i + 1, A.Ycount] = xyavg(A[i, A.Ycount - 1], A[i + 1, A.Ycount - 1]);
            }
            for (int j = 0; j < A.Ycount-1; j++)
            {
                newgrid[0, j + 1] = xyavg(A[0, j], A[0, j + 1]);
                newgrid[A.Xcount, j + 1] = xyavg(A[A.Xcount - 1, j], A[A.Xcount - 1, j + 1]);
            }
            for (int i = 0; i < A.Xcount-1; i++)
            {
                for (int j = 0; j < A.Ycount-1; j++)
                {
                    newgrid[i + 1, j + 1] = xyavg(A[i, j], A[i + 1, j], A[i, j + 1], A[i + 1, j + 1]);
                }
            }
            double accumulator = 0;
            for (int i = 0; i < m-1; i++)
            {
                for (int j = 0; j < n-1; j++)
                {
                    NCAMRNode[] nodes =
                    {
                        newgrid[i,j],
                        newgrid[i+1,j],
                        newgrid[i+1,j+1],
                        newgrid[i,j+1]
                    };
                    double[] x = new double[nodes.Length];
                    double[] y = new double[nodes.Length];
                    for (int z = 0; z < nodes.Length; z++)
                    {
                        x[z] = nodes[z].X;
                        y[z] = nodes[z].Y;
                    }
                    accumulator += quad_area(x, y)*(A[i,j].Value-B[i,j].Value)*(A[i,j].Value-B[i,j].Value);
                }
            }
            return Math.Sqrt(accumulator);
        }
        private static double area_of_influence_region(NCGrid_Distribution A, int i, int j)
        {
            int imax = A.Xcount - 1;
            int jmax = A.Ycount - 1;
            double[] x = new double[4];
            double[] y = new double[4];
            NCAMRNode[] grid =
            {
                A.GetFake(i-1,j+1),
                A.GetFake(i,j+1),
                A.GetFake(i+1,j+1),
                A.GetFake(i-1,j),
                A.GetFake(i,j),
                A.GetFake(i+1,j),
                A.GetFake(i-1,j-1),
                A.GetFake(i,j-1),
                A.GetFake(i+1,j-1)
            };
            NCAMRNode[] poly =
            {
                xyavg(grid[3], grid[4], grid[6], grid[7]),
                xyavg(grid[4], grid[5], grid[7], grid[8]),
                xyavg(grid[1], grid[2], grid[4], grid[5]),
                xyavg(grid[0], grid[1], grid[3], grid[4])
            };
            for (int z = 0; z < poly.Length; z++)
            {
                x[z] = poly[z].X;
                y[z] = poly[z].Y;
            }
            return quad_area(x, y);
        }
        private static NCAMRNode xyavg(params NCAMRNode[] nodes)
        {
            int count = nodes.Length;
            double xacc = 0;
            double yacc = 0;
            for (int i = 0; i < count; i++)
            {
                xacc += nodes[i].X;
                yacc += nodes[i].Y;
            }
            return new NCAMRNode(xacc / count, yacc / count, 0);
        }
        private static double quad_area(double[] x_clockwise, double[] y_clockwise)
        {
            Point3[] poly =
            {
                new Point3(x_clockwise[0], y_clockwise[0], 0),
                new Point3(x_clockwise[1], y_clockwise[1], 0),
                new Point3(x_clockwise[2], y_clockwise[2], 0),
                new Point3(x_clockwise[3], y_clockwise[3], 0)
            };
            Vector3 A1 = new Vector3(poly[0], poly[1]);
            Vector3 A2 = new Vector3(poly[0], poly[3]);
            Vector3 B1 = new Vector3(poly[2], poly[1]);
            Vector3 B2 = new Vector3(poly[2], poly[3]);
            return 0.5 * ((A1 % A2).Norm + (B1 % B2).Norm);
        }
    }
}
