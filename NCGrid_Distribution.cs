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
    class NCGrid_Distribution
    {
        private Random R = new Random();
        private RBounds2D bounds;
        private double dx, dy;
        private const int additional_file_lines = 2;
        private int x_node_count, y_node_count;
        private double minval, maxval;
        private NCAMRNode[,] nodes;
        public RBounds2D Bounds
        {
            get { return bounds; }
        }
        public enum DerivativeEstimationMode
        {
            WEIGHTFUNCTION,
            GFDM,
            TAYLOR_SURPLUS_UPPER_RIGHT,
            NAIVE
        }
        public enum SurplusNodeAccessingMode
        {
            UPPER_RIGHT,
            RANDOM,
            MINIMAL_DISTANCE
        }
        public double MaxValue
        {
            get { return maxval; }
        }
        public double MinValue
        {
            get { return minval; }
        }
        public int Xcount
        {
            get { return x_node_count; }
        }
        public int Ycount
        {
            get { return y_node_count; }
        }
        public double Xmin
        {
            get { return bounds.Xmin; }
            set { bounds.Xmin = value; }
        }
        public double Xmax
        {
            get { return bounds.Xmax; }
            set { bounds.Xmax = value; }
        }
        public double Ymin
        {
            get { return bounds.Ymin; }
            set { bounds.Ymin = value; }
        }
        public double Ymax
        {
            get { return bounds.Ymax; }
            set { bounds.Ymax = value; }
        }
        public static NCGrid_Distribution FirstAvailable()
        {
            string[] allnames = Directory.GetFiles(Paths.DistributionRepo);
            if (allnames.Length != 0)
            {
                return from_file(allnames[0], true);
            }
            else
            {
                throw new Exception("No files available.");
            }
        }
        public double interpolate_at(double x, double y)
        {
            if (x < bounds.Xmin || x > bounds.Xmax || y < bounds.Ymin || y > bounds.Ymax) { throw new Exception("Error: Extrapolation prevented."); }
            else
            {
                int i = get_x_ll_index(x, y);
                int j = get_y_ll_index(x, y);
                if (i >= Xcount - 1) { i--; }
                if (j >= Ycount - 1) { j--; }
                NCAMRNode n1 = nodes[i, j];
                NCAMRNode n2 = nodes[i + 1, j];
                NCAMRNode n3 = nodes[i + 1, j + 1];
                NCAMRNode n4 = nodes[i, j + 1];
                double[] weights =
                {
                    Math.Exp(-sq(x - n1.X) - sq(y - n1.Y)),
                    Math.Exp(-sq(x - n2.X) - sq(y - n2.Y)),
                    Math.Exp(-sq(x - n3.X) - sq(y - n3.Y)),
                    Math.Exp(-sq(x - n4.X) - sq(y - n4.Y))
                };
                double[] vals =
                {
                    n1.Value,
                    n2.Value,
                    n3.Value,
                    n4.Value
                };
                double acc_top = 0;
                double acc_bottom = 0;
                for (int q = 0; q < vals.Length; q++)
                {
                    acc_bottom += weights[q];
                    acc_top += weights[q] * vals[q];
                }
                return acc_top / acc_bottom;
            }
        }
        private int get_x_ll_index(double x, double y)
        {
            int upper = Xcount;
            int lower = 0;
            int med = (upper + lower) / 2;
            while (upper - lower > 1)
            {
                if (in_x_polygon_assume_right(med, x, y))
                {
                    upper = med;
                    med = (upper + lower) / 2;
                }
                else
                {
                    lower = med;
                    med = (upper + lower) / 2;
                }
            }
            return lower;
        }
        private int get_y_ll_index(double x, double y)
        {
            int upper = Ycount;
            int lower = 0;
            int med = (upper + lower) / 2;
            while (upper - lower > 1)
            {
                if (in_y_polygon_assume_right(med, x, y))
                {
                    upper = med;
                    med = (upper + lower) / 2;
                }
                else
                {
                    lower = med;
                    med = (upper + lower) / 2;
                }
            }
            return lower;
        }
        private double sq(double toSq)
        {
            return toSq * toSq;
        }
        private bool in_x_polygon_assume_right(int lowerbound, double x, double y)
        {
            bool odd = false;
            for (int i = 0; i < Ycount; i++)
            {
                NCAMRNode tracker = nodes[lowerbound, i];
                double track_x = tracker.X;
                double track_y = tracker.Y;
                if (((track_y > y) != odd) && track_x > x)
                {
                    odd = !odd;
                }
            }
            return odd;
        }
        private bool in_y_polygon_assume_right(int lowerbound, double x, double y)
        {
            bool odd = false;
            for (int i = 0; i < Xcount; i++)
            {
                NCAMRNode tracker = nodes[i, lowerbound];
                double track_x = tracker.X;
                double track_y = tracker.Y;
                if (((track_x > x) != odd) && track_y > y)
                {
                    odd = !odd;
                }
            }
            return odd;
        }
        public NCGrid_Distribution(RBounds2D _bounds, int xcount, int ycount)
        {
            bounds = _bounds;
            x_node_count = xcount;
            y_node_count = ycount;
            nodes = new NCAMRNode[x_node_count, y_node_count];
            generate_nodes(0);
            minval = 0;
            maxval = 0;
        }
        public NCGrid_Distribution(RBounds2D _bounds, int xcount, int ycount, double presetvalue)
        {
            bounds = _bounds;
            x_node_count = xcount;
            y_node_count = ycount;
            nodes = new NCAMRNode[x_node_count, y_node_count];
            generate_nodes(presetvalue);
            minval = presetvalue;
            maxval = presetvalue;
        }
        public void ApplyMeshMorphGA(double size)
        {
            int surplus_i = 1;
            int surplus_j = 1;
            for (int i = 1; i < x_node_count-1; i++)
            {
                for (int j = 1; j < y_node_count-1; j++)
                { 
                    NCAMRNode[] stencil =
                                {
                                nodes[i,j],
                                nodes[i+1,j],
                                nodes[i,j+1],
                                nodes[i-1,j],
                                nodes[i,j-1],
                                nodes[i+surplus_i,j+surplus_j]
                            };
                    TaylorSystem t_system = new TaylorSystem(stencil);
                    Matrix Astar = t_system.TaylorCoefficients;
                    double k = 200;
                    double[] kstuff = { k, k, k * k, k * k, k * k };
                    for (int c = 0; c < 5; c++)
                    {
                        for (int d = 0; d < 5; d++)
                        {
                            Astar[d, c] = Astar[d, c] * kstuff[c];
                        }
                    }
                    Matrix K = new Matrix(5);
                    K[0, 0] = k;
                    K[1, 1] = k;
                    K[2, 2] = k * k;
                    K[3, 3] = k * k;
                    K[4, 4] = k * k;
                    Matrix b =  K * Astar.Inverse;
                    Matrix xi = b * t_system.RHS;
                    double ux = xi[0];
                    double uxx = xi[2];
                    double uy = xi[1];
                    double uyy = xi[3];
                    double uxy = xi[4];
                    double _dx = ((uxx * ux) + (uxy * uy));
                    double _dy = ((uyy * uy) + (ux * uxy));
                    //normalize to max of gradient?
                    double step_size = size;
                    
                    Vector3 move = new Vector3(step_size * _dx, step_size * _dy, 0);
                    double xnew = (this[i, j] + move).X;
                    double ynew = (this[i, j] + move).Y;
                    while (xnew < bounds.Xmin || xnew > bounds.Xmax || ynew > bounds.Ymax || ynew < bounds.Ymin || move.Norm > 0.2*(dx+dy))
                    {
                        move = 0.5 * move;
                        xnew = (this[i, j] + move).X;
                        ynew = (this[i, j] + move).Y;
                    }
                    this[i, j] = this[i, j] + move;
                }
            }
        }
        public void ApplyMeshMorphGA(int iterations, double size)
        {
            for (int i = 0; i < iterations; i++)
            {
                ApplyMeshMorphGA(size);
            }
        }
        public NCGrid_Distribution Clone()
        {
            NCGrid_Distribution n = new NCGrid_Distribution(bounds, x_node_count, y_node_count);
            for (int i = 0; i < x_node_count; i++)
            {
                for (int j = 0; j < y_node_count; j++)
                {
                    n.assign_value_at(i, j, this[i, j].Value);
                    n[i, j].X = this[i, j].X;
                    n[i, j].Y = this[i, j].Y;
                }
            }
            return n;
        }
        private void generate_nodes(double preset_value)
        {
            dx = (bounds.Xmax - bounds.Xmin) / (x_node_count - 1);
            dy = (bounds.Ymax - bounds.Ymin) / (y_node_count - 1);
            minval = preset_value;

            maxval = preset_value;
            for (int i = 0; i < x_node_count; i++)
            {
                for (int j = 0; j < y_node_count; j++)
                {
                    NCAMRNode newnode = new NCAMRNode(bounds.Xmin + i * dx, bounds.Ymin + j * dy, preset_value);
                    nodes[i, j] = newnode;
                }
            }
        }
        public void WriteToFile(string title)
        {
            bool okaytowrite = !File.Exists(Paths.DistributionRepo + "\\" + title + ".dist");
            if (okaytowrite) { WriteToFile(title, false, false, false); }
            else
            {
                int n = 0;
                okaytowrite = !File.Exists(Paths.DistributionRepo + "\\" + title + n.ToString() + ".dist");
                while (!okaytowrite)
                {
                    n++;
                    okaytowrite = !File.Exists(Paths.DistributionRepo + "\\" + title + n.ToString() + ".dist");
                }
                WriteToFile(title + n.ToString(), false, false, false);
            }
            
        }
        public void WriteToFile(string title, bool enable_console_output, bool overwrite, bool using_full_path)
        {
            string path = Paths.DistributionRepo + "\\" + title + ".dist";
            if (using_full_path) { path = title; }
            if (File.Exists(path) && !overwrite) { throw new Exception("File exists, overwriting permission not explicitly granted"); }
            DateTime then = DateTime.Now;
            string[] contents = new string[x_node_count * y_node_count + additional_file_lines];
            contents[0] = bounds.ToString() + ":max_value:" + maxval.ToString() + ":min_value:" + minval.ToString();
            contents[1] = "Nx:" + x_node_count.ToString() + ":Ny:" + y_node_count.ToString();
            for (int i = 0; i < x_node_count; i++)
            {
                for (int j = 0; j < y_node_count; j++)
                {
                    contents[additional_file_lines + i + x_node_count * j] = nodes[i, j].ToString();
                }
            }
            File.WriteAllLines(path, contents);
            DateTime now = DateTime.Now;
            if (enable_console_output) { Console.WriteLine("File output in " + (now - then).Milliseconds.ToString() + " milliseconds."); }
        }
        public Vector3 estimate_num_graident(double x, double y, DerivativeEstimationMode mode)
        {
            //will need to add in edge cases. Ignore for now.
            switch (mode)
            {
                case DerivativeEstimationMode.GFDM:
                    {
                        return Vector3.Zero;
                    }
                case DerivativeEstimationMode.WEIGHTFUNCTION:
                    {
                        double delta_x = (bounds.Xmax - bounds.Xmin) / (Xcount - 1);
                        double delta_y = (bounds.Ymax - bounds.Ymin) / (Ycount - 1);
                        double dzdx = (this.interpolate_at(x + delta_x, y) - this.interpolate_at(x - delta_x, y)) / (2 * delta_x);
                        double dzdy = (this.interpolate_at(x, y + delta_y) - this.interpolate_at(x, y - delta_y)) / (2 * delta_x);
                        return new Vector3(dzdx, dzdy, 0);
                    }
                default:
                    {
                        return Vector3.Zero;
                    }
            }
        }
        public Vector3 estimate_num_graident(int i, int j, DerivativeEstimationMode mode)
        {
            //will need to add in edge cases. Ignore for now.
            switch (mode)
            {
                case DerivativeEstimationMode.GFDM:
                    {
                        //lacks implementation
                        return Vector3.Zero;
                    }
                case DerivativeEstimationMode.TAYLOR_SURPLUS_UPPER_RIGHT:
                    {
                        //clean this up???
                        NCAMRNode[] stencil =
                        {
                            nodes[i,j],
                            nodes[i+1,j],
                            nodes[i,j+1],
                            nodes[i-1,j],
                            nodes[i,j-1],
                            nodes[i+1,j+1]
                        };
                        double[] deltas_u =
                        {
                            stencil[1].Value-stencil[0].Value,
                            stencil[2].Value-stencil[0].Value,
                            stencil[3].Value-stencil[0].Value,
                            stencil[4].Value-stencil[0].Value,
                            stencil[5].Value-stencil[0].Value
                        };
                        double[] deltas_x =
                        {
                            stencil[1].X-stencil[0].X,
                            stencil[2].X-stencil[0].X,
                            stencil[3].X-stencil[0].X,
                            stencil[4].X-stencil[0].X,
                            stencil[5].X-stencil[0].X
                        };
                        double[] deltas_y =
                        {
                            stencil[1].Y-stencil[0].Y,
                            stencil[2].Y-stencil[0].Y,
                            stencil[3].Y-stencil[0].Y,
                            stencil[4].Y-stencil[0].Y,
                            stencil[5].Y-stencil[0].Y
                        };
                        double[] matrix_contents =
                        {
                            deltas_x[0], deltas_y[0], 0.5*sq(deltas_x[0]), 0.5*sq(deltas_x[0]), deltas_x[0]*deltas_y[0],
                            deltas_x[1], deltas_y[1], 0.5*sq(deltas_x[1]), 0.5*sq(deltas_x[1]), deltas_x[1]*deltas_y[1],
                            deltas_x[2], deltas_y[2], 0.5*sq(deltas_x[2]), 0.5*sq(deltas_x[2]), deltas_x[2]*deltas_y[2],
                            deltas_x[3], deltas_y[3], 0.5*sq(deltas_x[3]), 0.5*sq(deltas_x[3]), deltas_x[3]*deltas_y[3],
                            deltas_x[4], deltas_y[4], 0.5*sq(deltas_x[4]), 0.5*sq(deltas_x[4]), deltas_x[4]*deltas_y[4]
                        };
                        Matrix delta = new Matrix(5, 1);
                        return Vector3.Zero;
                    }
                case DerivativeEstimationMode.WEIGHTFUNCTION:
                    {
                        return Vector3.Zero;
                    }
                case DerivativeEstimationMode.NAIVE:
                    {
                        double delta_x = (bounds.Xmax - bounds.Xmin) / (Xcount - 1);
                        double delta_y = (bounds.Ymax - bounds.Ymin) / (Ycount - 1);
                        NCAMRNode node = nodes[i, j];
                        double dzdx = (this.interpolate_at(node.X + delta_x, node.Y) - this.interpolate_at(node.X - delta_x, node.Y)) / (2 * delta_x);
                        double dzdy = (this.interpolate_at(node.X, node.Y + delta_y) - this.interpolate_at(node.X, node.Y - delta_y)) / (2 * delta_x);
                        return new Vector3(dzdx, dzdy, 0);
                    }
                default:
                    {
                        return Vector3.Zero;
                    }
            }
        }
        public NCAMRNode this[int i, int j]
        {
            get { return nodes[i, j]; }
            set
            {
                nodes[i, j] = value;
                if (value.Value > maxval)
                {
                    maxval = value.Value;
                }
                if (value.Value < minval)
                {
                    minval = value.Value;
                }
            }
        }
        public void QuickSketch(string title)
        {
            force_extrema_update();
            DistributionSketchSettings S = DistributionSketchSettings.Fancy();
            DistributionSketch2D f = new DistributionSketch2D(this, S);
            f.CreateSketch(false);
            bool oktowrite = !File.Exists(Paths.ImageRepo + "\\" + title + ".bmp");
            if (oktowrite)
            {
                f.SaveImage(title, false);
            }
            else
            {
                int n = 0;
                oktowrite = !File.Exists(Paths.ImageRepo + "\\" + title + n.ToString() + ".bmp");
                while (!oktowrite)
                {
                    n++;
                    oktowrite = !File.Exists(Paths.ImageRepo + "\\" + title + n.ToString() + ".bmp");
                }
                f.SaveImage(title + n.ToString(), false);
            }
        }
        public void QuickSketch(string title, bool using_full_path)
        {
            force_extrema_update();
            DistributionSketchSettings S = DistributionSketchSettings.Fancy();
            DistributionSketch2D f = new DistributionSketch2D(this, S);
            f.CreateSketch(false);
            f.SaveImage(title, using_full_path);
        }
        public void force_extrema_update()
        {
            maxval = double.NegativeInfinity;
            minval = double.PositiveInfinity;
            for (int i = 0; i < Xcount; i++)
            {
                for (int j = 0; j < Ycount; j++)
                {
                    if (this[i, j].Value > maxval) { maxval = this[i, j].Value; }
                    if (this[i, j].Value < minval) { minval = this[i, j].Value; }
                }
            }
        }
        public static NCGrid_Distribution MakeMagnitude(NCGrid_Distribution A, NCGrid_Distribution B)
        {
            NCGrid_Distribution output = A.Clone();
            for (int i = 0; i < A.Xcount; i++)
            {
                for (int j = 0; j < B.Ycount; j++)
                {
                    double x = A[i, j].Value;
                    double y = B[i, j].Value;
                    output.assign_value_at(i, j, Math.Sqrt((x * x) + (y * y)));
                }
            }
            output.force_extrema_update();
            return output;
        }
        public void MimicMorph(NCGrid_Distribution template)
        {
            for (int i = 0; i < x_node_count; i++)
            {
                for (int j = 0; j < y_node_count; j++)
                {
                    this[i, j].X = template[i, j].X;
                    this[i, j].Y = template[i, j].Y;
                }
            }
        }
        public static NCGrid_Distribution from_file(string title, bool using_full_path)
        {
            string path = Paths.DistributionRepo + "\\" + title + ".dist";
            if (using_full_path) { path = title; }
            string[] contents = File.ReadAllLines(path);
            int length = contents.Length;
            RBounds2D new_bounds = RBounds2D.fromstring(contents[0]);
            string[] numbers = contents[1].Split(':');
            int nx, ny;
            if (!(int.TryParse(numbers[1], out nx) && int.TryParse(numbers[3], out ny))) { throw new Exception("File format error"); }
            NCGrid_Distribution created = new NCGrid_Distribution(new_bounds, nx, ny);
            for (int q = additional_file_lines; q < contents.Length; q++)
            {
                int qprime = q - additional_file_lines;
                int i = qprime % nx;
                int j = (qprime - i) / nx;
                created[i, j] = NCAMRNode.fromstring(contents[q]);
            }
            return created;
        }
        public void SetConstant(double constant)
        {
            maxval = constant;
            minval = constant;
            for (int i = 0; i < x_node_count; i++)
            {
                for (int j = 0; j < y_node_count; j++)
                {
                    this[i, j].Value = constant;
                }
            }
        }
        public Matrix GetTaylorSystemCoeffs(int i, int j, int surplus_i, int surplus_j)
        {
            NCAMRNode[] stencil =
                        {
                            nodes[i,j],
                            nodes[i+1,j],
                            nodes[i,j+1],
                            nodes[i-1,j],
                            nodes[i,j-1],
                            nodes[i+surplus_i,j+surplus_j]
                        };
            TaylorSystem t_system = new TaylorSystem(stencil);
            Matrix Astar = t_system.TaylorCoefficients;
            double k = 20;
            double[] kstuff = { k, k, k * k, k * k, k * k };
            for (int c = 0; c < 5; c++)
            {
                for (int d = 0; d < 5; d++)
                {
                    Astar[d, c] = Astar[d, c] * kstuff[c];
                }
            }
            Matrix K = new Matrix(5);
            K[0, 0] = k;
            K[1, 1] = k;
            K[2, 2] = k * k;
            K[3, 3] = k * k;
            K[4, 4] = k * k;
            return K * Astar.Inverse;
        }
        public Matrix EstimateSecondDerivs(int i, int j, SurplusNodeAccessingMode mode)
        {
            int surplus_i = 1;
            int surplus_j = 1;
            switch (mode)
            {
                case SurplusNodeAccessingMode.RANDOM:
                    {
                        surplus_i = 1 - 2 * R.Next(0, 1);
                        surplus_j = 1 - 2 * R.Next(0, 1);
                        break;
                    }
            }
            NCAMRNode[] stencil =
                        {
                            nodes[i,j],
                            nodes[i+1,j],
                            nodes[i,j+1],
                            nodes[i-1,j],
                            nodes[i,j-1],
                            nodes[i+surplus_i,j+surplus_j]
                        };
            TaylorSystem t_system = new TaylorSystem(stencil);
            Matrix Astar = t_system.TaylorCoefficients;
            double k = 9999;
            double[] kstuff = { k, k, k * k, k * k, k * k };
            for (int c = 0; c < 5; c++)
            {
                for (int d = 0; d < 5; d++)
                {
                    Astar[d, c] = Astar[d, c] * kstuff[c];
                }
            }
            Matrix K = new Matrix(5);
            K[0, 0] = k;
            K[1, 1] = k;
            K[2, 2] = k * k;
            K[3, 3] = k * k;
            K[4, 4] = k * k;
            return K * Astar.Inverse * t_system.RHS;
        }
        public void WriteAllDerivs(string title, bool enable_console_output, bool using_full_path, bool allow_overwrite)
        {
            DateTime Then = DateTime.Now;
            NCGrid_Distribution this_x = new NCGrid_Distribution(bounds, this.Xcount, this.Ycount);
            NCGrid_Distribution this_y = new NCGrid_Distribution(bounds, this.Xcount, this.Ycount);
            NCGrid_Distribution this_xx = new NCGrid_Distribution(bounds, this.Xcount, this.Ycount);
            NCGrid_Distribution this_yy = new NCGrid_Distribution(bounds, this.Xcount, this.Ycount);
            NCGrid_Distribution this_xy = new NCGrid_Distribution(bounds, this.Xcount, this.Ycount);
            for (int i = 1; i < this.Xcount - 1; i++)
            {
                for (int j = 1; j < this.Ycount - 1; j++)
                {
                    Matrix xi = EstimateSecondDerivs(i, j, SurplusNodeAccessingMode.UPPER_RIGHT);
                    this_x[i, j].Value = xi[0];
                    this_y[i, j].Value = xi[1];
                    this_xx[i, j].Value = xi[2];
                    this_yy[i, j].Value = xi[3];
                    this_xy[i, j].Value = xi[4];
                }
            }
            DateTime Now = DateTime.Now;
            string title_x = title + "_x";
            string title_y = title + "_y";
            string title_xx = title + "_xx";
            string title_yy = title + "_yy";
            string title_xy = title + "_xy";
            if (using_full_path)
            {
                string basestring = title.Split('.')[0];
                title_x = basestring + "_x.dist";
                title_y = basestring + "_y.dist";
                title_xx = basestring + "_xx.dist";
                title_yy = basestring + "_yy.dist";
                title_xy = basestring + "_xy.dist";
            }
            this_x.WriteToFile(title_x, enable_console_output, allow_overwrite, using_full_path);
            this_y.WriteToFile(title_y, enable_console_output, allow_overwrite, using_full_path);
            this_xx.WriteToFile(title_xx, enable_console_output, allow_overwrite, using_full_path);
            this_yy.WriteToFile(title_yy, enable_console_output, allow_overwrite, using_full_path);
            this_xy.WriteToFile(title_xy, enable_console_output, allow_overwrite, using_full_path);
            if (enable_console_output) { Console.WriteLine((Now - Then).TotalMilliseconds); }
        }
        public static NCGrid_Distribution operator +(NCGrid_Distribution A, NCGrid_Distribution B)
        {
            if (A.Xcount == B.Xcount && A.Ycount == B.Ycount)
            {
                NCGrid_Distribution C = new NCGrid_Distribution(A.Bounds, A.Xcount, B.Ycount);
                for (int i = 0; i < A.Xcount; i++)
                {
                    for (int j = 0; j < A.Ycount; j++)
                    {
                        C[i, j] = A[i, j].Clone();
                        C[i, j].Value = A[i, j].Value + B[i, j].Value;
                    }
                }
                return C;
            }
            else
            {
                throw new Exception("Dimensions do not agree.");
            }
        }
        public static NCGrid_Distribution operator -(NCGrid_Distribution A, NCGrid_Distribution B)
        {
            if (A.Xcount == B.Xcount && A.Ycount == B.Ycount)
            {
                NCGrid_Distribution C = new NCGrid_Distribution(A.Bounds, A.Xcount, B.Ycount);
                for (int i = 0; i < A.Xcount; i++)
                {
                    for (int j = 0; j < A.Ycount; j++)
                    {
                        C[i, j] = A[i, j].Clone();
                        C[i, j].Value = A[i, j].Value - B[i, j].Value;
                    }
                }
                return C;
            }
            else
            {
                throw new Exception("Dimensions do not agree.");
            }
        }
        public static NCGrid_Distribution operator *(NCGrid_Distribution A, NCGrid_Distribution B)
        {
            if (A.Xcount == B.Xcount && A.Ycount == B.Ycount)
            {
                NCGrid_Distribution C = new NCGrid_Distribution(A.Bounds, A.Xcount, B.Ycount);
                for (int i = 0; i < A.Xcount; i++)
                {
                    for (int j = 0; j < A.Ycount; j++)
                    {
                        C[i, j] = A[i, j].Clone();
                        C[i, j].Value = A[i, j].Value * B[i, j].Value;
                    }
                }
                return C;
            }
            else
            {
                throw new Exception("Dimensions do not agree.");
            }
        }
        public static NCGrid_Distribution operator /(NCGrid_Distribution A, NCGrid_Distribution B)
        {
            if (A.Xcount == B.Xcount && A.Ycount == B.Ycount)
            {
                NCGrid_Distribution C = new NCGrid_Distribution(A.Bounds, A.Xcount, B.Ycount);
                for (int i = 0; i < A.Xcount; i++)
                {
                    for (int j = 0; j < A.Ycount; j++)
                    {
                        C[i, j] = A[i, j].Clone();
                        C[i, j].Value = A[i, j].Value / B[i, j].Value;
                    }
                }
                return C;
            }
            else
            {
                throw new Exception("Dimensions do not agree.");
            }
        }
        public static NCGrid_Distribution operator *(double scale, NCGrid_Distribution A)
        {
            NCGrid_Distribution output = A.Clone();
            for (int i = 0; i < A.Xcount; i++)
            {
                for (int j = 0; j < A.Ycount; j++)
                {
                    output.assign_value_at(i, j, scale * A[i, j].Value);
                }
            }
            return output;
        }
        public static NCGrid_Distribution NiceFunction(string vb_syntax_eval, RBounds2D bounds, int xcount, int ycount, bool keepfile, string title)
        {
            ExactFunctionGeneratorVB2D.GenerateFunction(vb_syntax_eval, bounds, xcount, ycount, title);
            NCGrid_Distribution n = NCGrid_Distribution.from_file(title, false);
            return n;
        }
        public static NCGrid_Distribution NiceFunction(string vb_syntax_eval, RBounds2D bounds, int xcount, int ycount)
        {
            ExactFunctionGeneratorVB2D.GenerateFunction(vb_syntax_eval, bounds, xcount, ycount, "temp");
            NCGrid_Distribution n = NCGrid_Distribution.from_file("temp", false);
            string to_delete = Paths.DistributionRepo + "\\temp.dist";
            File.Delete(to_delete);
            return n;
        }
        public void ApplyBoundary(BoundaryConditions B)
        {
            if (!check_compatibility(B, this)) { throw new Exception("Incompatible conditions."); }
            for (int i = 0; i < B.Xcount; i++)
            {
                double byn = B[i, BoundaryConditions.Direction.Negative_Y];
                double byp = B[i, BoundaryConditions.Direction.Positive_Y];
                this[i, 0].Value = byn;
                this[i, B.Xcount - 1].Value = byp;
                if (byp > maxval) { maxval = byp; }
                if (byn > maxval) { maxval = byn; }
                if (byp < minval) { minval = byp; }
                if (byn < minval) { minval = byn; }
            }
            for (int j = 0; j < B.Ycount; j++)
            {
                double bxp = B[j, BoundaryConditions.Direction.Positive_X];
                double bxn = B[j, BoundaryConditions.Direction.Negative_X];
                this[0, j].Value = bxn;
                this[B.Ycount - 1, j].Value = bxp;
                if (bxp > maxval) { maxval = bxp; }
                if (bxn > maxval) { maxval = bxn; }
                if (bxp < minval) { minval = bxp; }
                if (bxn < minval) { minval = bxn; }
            }
        }
        public void assign_value_at(int i, int j, double valuein)
        {
            this[i, j].Value = valuein;
            if (valuein > maxval) { maxval = valuein; }
            if (valuein < minval) { minval = valuein; }
        }
        public static bool check_compatibility(BoundaryConditions B, NCGrid_Distribution C)
        {
            bool[] stuff =
            {
                B.Bounds.Xmax == C.Xmax,
                B.Bounds.Xmin == C.Xmin,
                B.Bounds.Ymax == C.Ymax,
                B.Bounds.Ymin == C.Ymin,
                B.Ycount == C.Ycount,
                B.Xcount == C.Xcount
            };
            bool valid = true;
            foreach (bool i in stuff)
            {
                valid = valid && i;
            }
            return valid;
        }
        public NCAMRNode GetFake(int i, int j)
        {
            if (i >= 0 && i <= Xcount-1 && j >= 0 && j <= Ycount-1)
            {
                return this[i, j];
            }
            else
            {
                int imax = Xcount - 1;
                int jmax = Ycount - 1;
                double x = double.NegativeInfinity;
                double y = double.NegativeInfinity;
                if (i < 0 && j < 0) //ll
                {
                    x = Xmin;
                    y = Ymin;
                }
                if (i > imax && j < 0) //lr
                {
                    x = Xmax;
                    y = Ymin;
                }
                if (i > imax && j > jmax) //ur
                {
                    x = Xmax;
                    y = Ymax;
                }
                if (i < 0 && j > jmax) //ul
                {
                    y = Ymax;
                    x = Xmin;
                }
                if (i > imax && j <= jmax && j >= 0) //px
                {
                    x = Xmax;
                    y = this[imax, j].Y;
                }
                if (i < 0 && j <= jmax && j >= 0)
                {
                    x = Xmin;
                    y = this[0, j].Y;
                }
                if (j < 0 && i <= imax && i >= 0)
                {
                    x = this[i, 0].X;
                    y = Ymax;
                }
                if (j > jmax && i >= 0 && i <= imax)
                {
                    x = this[i, jmax].X;
                    y = Ymax;
                }
                return new NCAMRNode(x, y, 0);
            }
        }
    }
}
