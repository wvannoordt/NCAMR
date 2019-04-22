/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace m_435_NCAMR
{
    class FluidProperties
    {
        private double rho, mu, nu;
        public double Density
        {
            get { return rho; }
        }
        public double DynamicViscosity
        {
            get { return mu; }
        }
        public double KinematicViscosity
        {
            get { return nu; }
        }
        public FluidProperties(double _density, double dynamic_viscosity)
        {
            mu = dynamic_viscosity;
            rho = _density;
            nu = mu / rho;
        }
    }
}
