/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace m_435_NCAMR
{
    class RBounds2D
    {
        private double xmin, xmax, ymin, ymax;
        public double Xmin
        {
            get { return xmin; }
            set { xmin = value; }
        }
        public double Xmax
        {
            get { return xmax; }
            set { xmax = value; }
        }
        public double Ymin
        {
            get { return ymin; }
            set { ymin = value; }
        }
        public double Ymax
        {
            get { return ymax; }
            set { ymax = value; }
        }
        public RBounds2D(double _xmin, double _xmax, double _ymin, double _ymax)
        {
            xmin = _xmin;
            ymin = _ymin;
            xmax = _xmax;
            ymax = _ymax;
        }
        public override string ToString()
        {
            return "bounds:" + xmin.ToString() + ":" + xmax.ToString() + ":" + ymin.ToString() + ":" + ymax.ToString();
        }
        public static RBounds2D fromstring(string source)
        {
            string[] ar = source.Split(':');
            double xn, xx, yn, yx;
            if (!(double.TryParse(ar[1], out xn) && double.TryParse(ar[2], out xx) && double.TryParse(ar[3], out yn) && double.TryParse(ar[4], out yx))) { throw new ArgumentException("String improperly formatted."); }
            else { return new RBounds2D(xn, xx, yn, yx); }
        }
        public static RBounds2D Square(double radius)
        {
            return new RBounds2D(-radius, radius, -radius, radius);
        }

    }
}
