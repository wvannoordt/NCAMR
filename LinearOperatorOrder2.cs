/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace m_435_NCAMR
{
    class LinearOperatorOrder2
    {
        public enum OperatorType
        {
            Parabolic,
            Hyperbolic,
            Elliptic
        }
        OperatorType type;
        double a, a_x, a_y, a_xx, a_yy, a_xy;
        public OperatorType Classification
        {
            get { return type; }
        }
        public double A
        {
            get { return a; }
            set { a = value; determine_type(); }
        }
        public double Ax
        {
            get { return a_x; }
            set { a_x = value; determine_type(); }
        }
        public double Ay
        {
            get { return a_y; }
            set { a_y = value; determine_type(); }
        }
        public double Axx
        {
            get { return a_xx; }
            set { a_xx = value; determine_type(); }
        }
        public double Ayy
        {
            get { return a_yy; }
            set { a_yy = value; determine_type(); }
        }
        public double Axy
        {
            get { return a_xy; }
            set { a_xy = value; determine_type(); }
        }
        public LinearOperatorOrder2(double[] coeffs_6x)
        {
            if (coeffs_6x.Length != 6) { throw new Exception("Error: Please provide exactly 6 coefficients."); }
            a = coeffs_6x[0];
            a_x = coeffs_6x[1];
            a_y = coeffs_6x[2];
            a_xx = coeffs_6x[3];
            a_yy = coeffs_6x[4];
            a_xy = coeffs_6x[5];
            determine_type();
        }
        public LinearOperatorOrder2(double _a, double _a_x, double _a_y, double _a_xx, double _a_yy, double _a_xy)
        {
            a = _a;
            a_x = _a_x;
            a_y = _a_y;
            a_xx = _a_xx;
            a_yy = _a_yy;
            a_xy = _a_xy;
            determine_type();
        }
        private void determine_type()
        {
            double z = (a_xy * a_xy) - (4 * a_xx * a_yy);
            if (Math.Abs(z) < 1E-10) { type = OperatorType.Parabolic; }
            else if (z > 0) { type = OperatorType.Hyperbolic; }
            else { type = OperatorType.Elliptic; }
        }
        public static LinearOperatorOrder2 Laplace = new LinearOperatorOrder2(0, 0, 0, 1, 1, 0);
        public double this[int i]
        {
            get
            {
                switch (i)
                {
                    case -1:
                        {
                            return a;
                        }
                    case 0:
                        {
                            return a_x;
                        }
                    case 1:
                        {
                            return a_y;
                        }
                    case 2:
                        {
                            return a_xx;
                        }
                    case 3:
                        {
                            return a_yy;
                        }
                    case 4:
                        {
                            return a_xy;
                        }
                    default:
                        {
                            throw new Exception("Index out of range.");
                        }
                }
            }
            set
            {
                switch (i)
                {
                    case -1:
                        {
                            a = value;
                            break;
                        }
                    case 0:
                        {
                            a_x = value;
                            break;
                        }
                    case 1:
                        {
                            a_y = value;
                            break;
                        }
                    case 2:
                        {
                            a_xx = value;
                            break;
                        }
                    case 3:
                        {
                            a_yy = value;
                            break;
                        }
                    case 4:
                        {
                            a_xy = value;
                            break;
                        }
                    default:
                        {
                            throw new Exception("Index out of range.");
                        }
                }
            }
        }
    }
}
