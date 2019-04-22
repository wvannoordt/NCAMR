/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Threading.Tasks;

namespace m_435_NCAMR
{
    class BoundaryConditions
    {
        private string bc_repo = @"C:\Users\Will\Desktop\Folders\MATH435\repo\boundary-conditions";
        public enum Direction
        {
            Positive_X,
            Positive_Y,
            Negative_X,
            Negative_Y
        }
        public enum BoundaryConditionType
        {
            Dirichlet,
            Neumann
        }
        private RBounds2D bounds;
        private BoundaryConditionType positive_x;
        private BoundaryConditionType positive_y;
        private BoundaryConditionType negative_x;
        private BoundaryConditionType negative_y;
        public RBounds2D Bounds
        {
            get { return bounds; }
        }
        private int xcount, ycount;
        public int Xcount
        {
            get { return xcount; }
        }
        public int Ycount
        {
            get { return ycount; }
        }
        private double[] positive_x_vals;
        private double[] positive_y_vals;
        private double[] negative_x_vals;
        private double[] negative_y_vals;
        public double this[int i, Direction indicator]
        {
            get
            {
                switch (indicator)
                {
                    case Direction.Positive_X:
                        {
                            return positive_x_vals[i];
                        }
                    case Direction.Negative_X:
                        {
                            return negative_x_vals[i];
                        }
                    case Direction.Positive_Y:
                        {
                            return positive_y_vals[i];
                        }
                    case Direction.Negative_Y:
                        {
                            return negative_y_vals[i];
                        }
                    default:
                        {
                            throw new Exception("Invalid direction.");
                        }
                }
            }
            set
            {
                switch (indicator)
                {
                    case Direction.Positive_X:
                        {
                            positive_x_vals[i] = value; break;
                        }
                    case Direction.Negative_X:
                        {
                            negative_x_vals[i] = value; break;
                        }
                    case Direction.Positive_Y:
                        {
                            positive_y_vals[i] = value; break;
                        }
                    case Direction.Negative_Y:
                        {
                            negative_y_vals[i] = value; break;
                        }
                    default:
                        {
                            throw new Exception("Invalid direction.");
                        }
                }
            }
        }
        public BoundaryConditions(RBounds2D _bounds, int _xcount, int _ycount, BoundaryConditionType _type)
        {
            xcount = _xcount;
            ycount = _ycount;
            positive_x_vals = new double[xcount];
            positive_y_vals = new double[xcount];
            negative_x_vals = new double[ycount];
            negative_y_vals = new double[ycount];
            bounds = _bounds;
            SetAllBoundaryConditionTypes(_type);
        }
        public void SetAllBoundaryConditionTypes(BoundaryConditionType type)
        {
            positive_x = type;
            positive_y = type;
            negative_x = type;
            negative_y = type;
        }
        public BoundaryConditions(NCGrid_Distribution boundary_of, BoundaryConditionType _type)
        {
            bounds = boundary_of.Bounds;
            SetAllBoundaryConditionTypes(_type);
            xcount = boundary_of.Xcount;
            ycount = boundary_of.Ycount;
            positive_x_vals = new double[xcount];
            positive_y_vals = new double[xcount];
            negative_x_vals = new double[ycount];
            negative_y_vals = new double[ycount];
            for (int i = 0; i < xcount; i++)
            {
                positive_x_vals[i] = boundary_of[i, ycount - 1].Value;
                negative_x_vals[i] = boundary_of[i, 0].Value;
            }
            for (int j = 0; j < ycount; j++)
            {
                positive_y_vals[j] = boundary_of[xcount-1,j].Value;
                negative_y_vals[j] = boundary_of[0,j].Value;
            }
        }
        public void SetBoundaryConditionType(Direction indicator, BoundaryConditionType type)
        {
            switch (indicator)
            {
                case Direction.Positive_X:
                    {
                        positive_x = type;
                        break;
                    }
                case Direction.Negative_X:
                    {
                        negative_x = type;
                        break;
                    }
                case Direction.Positive_Y:
                    {
                        positive_y = type;
                        break;
                    }
                case Direction.Negative_Y:
                    {
                        negative_y = type;
                        break;
                    }
            }
        }
        public BoundaryConditionType GetBoundaryConditionType(Direction indicator)
        {
            switch (indicator)
            {
                case Direction.Positive_X:
                    {
                        return positive_x;
                    }
                case Direction.Negative_X:
                    {
                        return negative_x;
                    }
                case Direction.Positive_Y:
                    {
                        return positive_y;
                    }
                case Direction.Negative_Y:
                    {
                        return negative_y;
                    }
                default:
                    {
                        throw new Exception("Error: Invalid direction.");
                    }
            }
        }
        public void WriteFile(string title, bool using_full_path, bool allow_overwrite)
        {
            string path = bc_repo + "\\" + title + ".bc";
            if (using_full_path) { path = title; }
            if (File.Exists(path) && !allow_overwrite)
            {
                throw new Exception("Error: Overwriting permissions not granted.");
            }
            List<string> stuff = new List<string>();
            stuff.Add(bounds.ToString());
            stuff.Add(string.Format("Nx:{0}:Ny:{1}", xcount, ycount));

            string px = "positive_x:" + ConditionTypeToString(positive_x);
            stuff.Add(px);
            int px_count = positive_x_vals.Length;
            for (int i = 0; i < px_count; i++)
            {
                stuff.Add(positive_x_vals[i].ToString());
            }

            string nx = "negative_x:" + ConditionTypeToString(negative_x);
            stuff.Add(nx);
            int nx_count = negative_x_vals.Length;
            for (int i = 0; i < nx_count; i++)
            {
                stuff.Add(negative_x_vals[i].ToString());
            }

            string py = "positive_y:" + ConditionTypeToString(positive_y);
            stuff.Add(py);
            int py_count = positive_y_vals.Length;
            for (int i = 0; i < py_count; i++)
            {
                stuff.Add(positive_y_vals[i].ToString());
            }

            string ny = "negative_y:" + ConditionTypeToString(negative_y);
            stuff.Add(ny);
            int ny_count = negative_y_vals.Length;
            for (int i = 0; i < ny_count; i++)
            {
                stuff.Add(negative_y_vals[i].ToString());
            }

            File.WriteAllLines(path, stuff.ToArray());
        }
        public BoundaryConditions(string title, bool using_full_path)
        {
            string path = bc_repo + "//" + title + ".bc";
            if (using_full_path) { path = title; }
            string[] contents = File.ReadAllLines(path);
            bounds = RBounds2D.fromstring(contents[0]);
            if (!(int.TryParse(contents[1].Split(':')[1], out xcount) && int.TryParse(contents[1].Split(':')[3], out ycount))) { throw new Exception("Error: Improper file format."); }
            int px_offset = 3;
            int nx_offset = px_offset+1+xcount;
            int py_offset = nx_offset+1+xcount;
            int ny_offset = py_offset+1+ycount;
            int terminal = ny_offset+1+ycount;
            positive_x_vals = new double[xcount];
            positive_y_vals = new double[xcount];
            negative_x_vals = new double[ycount];
            negative_y_vals = new double[ycount];
            positive_x = ConditionTypeFromString(contents[px_offset - 1].Split(':')[1]);
            negative_x = ConditionTypeFromString(contents[nx_offset - 1].Split(':')[1]);
            positive_y = ConditionTypeFromString(contents[py_offset - 1].Split(':')[1]);
            negative_y = ConditionTypeFromString(contents[ny_offset - 1].Split(':')[1]);
            for (int i = px_offset; i < nx_offset-1; i++)
            {
                double z;
                if (!double.TryParse(contents[i], out z)) { throw new Exception("Error: Improper file format"); }
                positive_x_vals[i - px_offset] = z;
            }
            for (int i = nx_offset; i < py_offset-1; i++)
            {
                double z;
                if (!double.TryParse(contents[i], out z)) { throw new Exception("Error: Improper file format"); }
                negative_x_vals[i - nx_offset] = z;
            }
            for (int i = py_offset; i < ny_offset-1; i++)
            {
                double z;
                if (!double.TryParse(contents[i], out z)) { throw new Exception("Error: Improper file format"); }
                positive_y_vals[i - py_offset] = z;
            }
            for (int i = ny_offset; i < terminal-1; i++)
            {
                double z;
                if (!double.TryParse(contents[i], out z)) { throw new Exception("Error: Improper file format"); }
                negative_y_vals[i - ny_offset] = z;
            }
        }
        public void SetConstant(double val, Direction indicator)
        {
            switch (indicator)
            {
                case Direction.Positive_X:
                    {
                        for (int i = 0; i < xcount; i++)
                        {
                            positive_x_vals[i] = val;
                        }
                        break;
                    }
                case Direction.Negative_X:
                    {
                        for (int i = 0; i < xcount; i++)
                        {
                            negative_x_vals[i] = val;
                        }
                        break;
                    }
                case Direction.Positive_Y:
                    {
                        for (int i = 0; i < ycount; i++)
                        {
                            positive_y_vals[i] = val;
                        }
                        break;
                    }
                case Direction.Negative_Y:
                    {
                        for (int i = 0; i < ycount; i++)
                        {
                            negative_y_vals[i] = val;
                        }
                        break;
                    }
            }
        }
        private string ConditionTypeToString(BoundaryConditionType type)
        {
            switch (type)
            {
                case BoundaryConditionType.Dirichlet:
                    {
                        return "dirichlet";
                    }
                case BoundaryConditionType.Neumann:
                    {
                        return "neumann";
                    }
                default:
                    {
                        throw new Exception("Error: Invalid type.");
                    }
            }
        }
        private BoundaryConditionType ConditionTypeFromString(string source)
        {
            switch (source)
            {
                case "dirichlet":
                    {
                        return BoundaryConditionType.Dirichlet;
                    }
                case "neumann":
                    {
                        return BoundaryConditionType.Neumann;
                    }
                default:
                    {
                        throw new Exception("Error: Invalid type.");
                    }
            }
        }
    }
}
