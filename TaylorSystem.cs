/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Matrix_v1;

namespace m_435_NCAMR
{
    class TaylorSystem
    {
        private Matrix taylor_coefs, rhs;
        public Matrix TaylorCoefficients
        {
            get { return taylor_coefs; }
        }
        public Matrix RHS
        {
            get { return rhs; }
        }
        public TaylorSystem(NCAMRNode[] stencil)
        {
            BuildTaylorCoeffs(stencil);
        }
        private void BuildTaylorCoeffs(NCAMRNode[] stencil)
        {
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
                            deltas_x[0], deltas_y[0], 0.5*sq(deltas_x[0]), 0.5*sq(deltas_y[0]), deltas_x[0]*deltas_y[0],
                            deltas_x[1], deltas_y[1], 0.5*sq(deltas_x[1]), 0.5*sq(deltas_y[1]), deltas_x[1]*deltas_y[1],
                            deltas_x[2], deltas_y[2], 0.5*sq(deltas_x[2]), 0.5*sq(deltas_y[2]), deltas_x[2]*deltas_y[2],
                            deltas_x[3], deltas_y[3], 0.5*sq(deltas_x[3]), 0.5*sq(deltas_y[3]), deltas_x[3]*deltas_y[3],
                            deltas_x[4], deltas_y[4], 0.5*sq(deltas_x[4]), 0.5*sq(deltas_y[4]), deltas_x[4]*deltas_y[4]
                        };
            Matrix delta = new Matrix(5, 1, deltas_u);
            Matrix A = new Matrix(5, 5, matrix_contents);
            taylor_coefs = A;
            rhs = delta;
            string[] filestuff = { taylor_coefs.Latexoutput(true) };
        }
        private double sq(double to)
        {
            return to * to;
        }
    }
}
