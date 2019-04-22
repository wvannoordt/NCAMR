/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace m_435_NCAMR
{
    class NavierStokesProblemSingleFrame
    {
        private int problemsize;
        private Matrix system_matrix;
        private NCGrid_Distribution prev_u, prev_v, u_star, v_star, discretization;
        private NCGrid_Distribution[] solution;
        public NCGrid_Distribution[] Solution { get { return solution; } }
        private Matrix current_solution_raw;
        private Matrix rhs;
        private FluidProperties f_props;
        int eq_nodes_x;
        int eq_nodes_y;
        double deltat;
        private BoundaryConditions u_boundary_conditions, v_boundary_conditions, p_boundary_conditions;
        public NavierStokesProblemSingleFrame(NCGrid_Distribution ustar, NCGrid_Distribution vstar, NCGrid_Distribution previous_solution_u, NCGrid_Distribution previous_solution_v, BoundaryConditions u, BoundaryConditions v, BoundaryConditions p, FluidProperties props, double dt)
        {
            deltat = dt;
            f_props = props;
            u_star = ustar;
            v_star = vstar;
            discretization = previous_solution_u.Clone();
            prev_u = previous_solution_u;
            prev_v = previous_solution_v;
            u_boundary_conditions = u;
            v_boundary_conditions = v;
            p_boundary_conditions = p;
            problemsize = 3*(discretization.Xcount - 2) * (discretization.Ycount - 2);
            system_matrix = new Matrix(problemsize);
            rhs = new Matrix(problemsize, 1);
            eq_nodes_x = discretization.Xcount - 2;
            eq_nodes_y = discretization.Ycount - 2;
        }
        public NCGrid_Distribution[] Solve()
        {
            build_system();
            for (int i = 0; i < 2; i++)//number of refinements = 3 for now
            {
                current_solution_raw = Matrix.solve_system(system_matrix, rhs,Matrix.SystemSolvingScheme.Kaczmarz);
                solution = make_from_matrix(current_solution_raw);
                NCGrid_Distribution magnitude = NCGrid_Distribution.MakeMagnitude(solution[1], solution[2]);
                magnitude.ApplyMeshMorphGA(0.2);
                solution[0].MimicMorph(magnitude);
                solution[1].MimicMorph(magnitude);
                solution[2].MimicMorph(magnitude);
                discretization.MimicMorph(magnitude);
                u_star = solution[1].Clone();
                v_star = solution[2].Clone();
                build_system();
            }
            current_solution_raw = Matrix.solve_system(system_matrix, rhs, Matrix.SystemSolvingScheme.Kaczmarz);
            solution = make_from_matrix(current_solution_raw);
            solution[0].MimicMorph(discretization);
            solution[1].MimicMorph(discretization);
            solution[2].MimicMorph(discretization);
            return solution;
        }
        private void build_system()
        {
            double rho = f_props.Density;
            double nu = f_props.KinematicViscosity;
            //interior
            for (int i_global = 1; i_global < discretization.Xcount - 1; i_global++)
            {
                for (int j_global = 1; j_global < discretization.Ycount - 1; j_global++)
                {
                    int i_local = i_global - 1;
                    int j_local = j_global - 1;
                    int current_row_cluster = j_local + (i_local * eq_nodes_y);
                    int force_x = 3 * current_row_cluster;
                    int force_y = force_x + 1;
                    int continuity = force_x + 2;
                    Matrix b = discretization.GetTaylorSystemCoeffs(i_global, j_global, 1, 1);
                    double[] sums = new double[5];
                    for (int q = 0; q < 5; q++)
                    {
                        sums[q] = 0;
                    }
                    for (int q = 0; q < 5; q++)
                    {
                        sums[0] += b[0, q];
                        sums[1] += b[1, q];
                        sums[2] += b[2, q];
                        sums[3] += b[3, q];
                        sums[4] += b[4, q];
                    }
                    double[] P = {
                        -1*sums[0]/f_props.Density,
                        b[0, 0] / rho,
                        b[0, 1] / rho,
                        b[0, 2] / rho,
                        b[0, 3] / rho,
                        b[0, 4] / rho
                    };
                    double[] Q = {
                        (-1*u_star[i_global,j_global].Value*sums[0]) - (v_star[i_global,j_global].Value*sums[1]) + (nu*sums[2])+(nu*sums[3]) - (1/deltat),
                        (u_star[i_global,j_global].Value*b[0,0]) + (v_star[i_global,j_global].Value*b[1,0])-(nu*b[2,0])-(nu*b[3,0]),
                        (u_star[i_global,j_global].Value*b[0,1]) + (v_star[i_global,j_global].Value*b[1,1])-(nu*b[2,1])-(nu*b[3,1]),
                        (u_star[i_global,j_global].Value*b[0,2]) + (v_star[i_global,j_global].Value*b[1,2])-(nu*b[2,2])-(nu*b[3,2]),
                        (u_star[i_global,j_global].Value*b[0,3]) + (v_star[i_global,j_global].Value*b[1,3])-(nu*b[2,3])-(nu*b[3,3]),
                        (u_star[i_global,j_global].Value*b[0,4]) + (v_star[i_global,j_global].Value*b[1,4])-(nu*b[2,4])-(nu*b[3,4])
                    };
                    double[] R = {
                        -1*sums[1]/f_props.Density,
                        b[1, 0] / rho,
                        b[1, 1] / rho,
                        b[1, 2] / rho,
                        b[1, 3] / rho,
                        b[1, 4] / rho
                    };
                    double[] S = Q;
                    double[] T = {
                        sums[0],
                        b[0,0],
                        b[0,1],
                        b[0,2],
                        b[0,3],
                        b[0,4]
                    };
                    double[] W = {
                        sums[1],
                        b[1,0],
                        b[1,1],
                        b[1,2],
                        b[1,3],
                        b[1,4]
                    };
                    //build pressure'
                    int rowbase = 3 * current_row_cluster;
                    if (i_global > 1 && i_global < discretization.Xcount - 2 && j_global > 1 && j_global < discretization.Ycount - 2)
                    {
                        map_interior_node_eqs(rowbase, i_local, j_local, P, Q, R, S, T, W);
                    }

                    if (i_global == 1 && j_global != 1 && j_global != discretization.Ycount - 2)
                    {
                        map_left_node_eqs(rowbase, i_local, j_local, P, Q, R, S, T, W);
                    }
                    if (i_global == discretization.Xcount - 2 && j_global != 1 && j_global != discretization.Ycount - 2)
                    {
                        map_right_node_eqs(rowbase, i_local, j_local, P, Q, R, S, T, W);
                    }

                    if (j_global == 1 && i_global != 1 && i_global != discretization.Xcount - 2)
                    {
                        map_lower_node_eqs(rowbase, i_local, j_local, P, Q, R, S, T, W);
                    }
                    if (j_global == discretization.Ycount - 2 && i_global != 1 && i_global != discretization.Xcount - 2)
                    {
                        map_upper_node_eqs(rowbase, i_local, j_local, P, Q, R, S, T, W);
                    }

                    if (j_global == 1 && i_global == 1)
                    {
                        map_llcorner_node_eqs(rowbase, i_local, j_local, P, Q, R, S, T, W);
                    }
                    if (j_global == 1 && i_global == discretization.Xcount - 2)
                    {
                        map_lrcorner_node_eqs(rowbase, i_local, j_local, P, Q, R, S, T, W);
                    }
                    if (j_global == discretization.Ycount - 2 && i_global == 1)
                    {
                        map_ulcorner_node_eqs(rowbase, i_local, j_local, P, Q, R, S, T, W);
                    }
                    if (j_global == discretization.Ycount - 2 && i_global == discretization.Xcount - 2)
                    {
                        map_urcorner_node_eqs(rowbase, i_local, j_local, P, Q, R, S, T, W);
                    }
                }
            }
        }
        private NCGrid_Distribution[] make_from_matrix(Matrix doink)
        {
            NCGrid_Distribution pout = prev_u.Clone();
            NCGrid_Distribution uout = prev_u.Clone();
            NCGrid_Distribution vout = prev_u.Clone();
            pout.ApplyBoundary(p_boundary_conditions);
            uout.ApplyBoundary(u_boundary_conditions);
            vout.ApplyBoundary(v_boundary_conditions);

            for (int i = 0; i < doink.Rows; i+= 3)
            {
                double pij = doink[i];
                double uij = doink[i+1];
                double vij = doink[i+2];
                int[] pcoords = inverse_map_var(i, 'p');
                int[] ucoords = inverse_map_var(i+1, 'u');
                int[] vcoords = inverse_map_var(i+2, 'v');
                pout.assign_value_at(pcoords[0], pcoords[1], pij);
                uout.assign_value_at(ucoords[0], ucoords[1], uij);
                vout.assign_value_at(vcoords[0], vcoords[1], vij);
            }
            NCGrid_Distribution[] output = new NCGrid_Distribution[3];
            output[0] = pout;
            output[1] = uout;
            output[2] = vout;
            return output;
        }
        private void map_right_node_eqs(int rowbase, int i_local, int j_local, double[] P, double[] Q, double[] R, double[] S, double[] T, double[] W)
        {
            int global_i = i_local + 1;
            int global_j = j_local + 1;
            rhs[rowbase] = (-1 / deltat) * prev_u[i_local + 1, j_local + 1].Value;
            rhs[rowbase + 1] = (-1 / deltat) * prev_v[i_local + 1, j_local + 1].Value;
            rhs[rowbase + 2] = 0;
            //x-momentum
            system_matrix[rowbase, map_var(i_local, j_local, 'p')] = P[0];
            //system_matrix[rowbase, map_var(i_local + 1, j_local, 'p')] = P[1];
            system_matrix[rowbase, map_var(i_local, j_local + 1, 'p')] = P[2];
            system_matrix[rowbase, map_var(i_local - 1, j_local, 'p')] = P[3];
            system_matrix[rowbase, map_var(i_local, j_local - 1, 'p')] = P[4];
            //system_matrix[rowbase, map_var(i_local + 1, j_local + 1, 'p')] = P[5];

            system_matrix[rowbase, map_var(i_local, j_local, 'u')] = Q[0];
            //system_matrix[rowbase, map_var(i_local + 1, j_local, 'u')] = Q[1];
            system_matrix[rowbase, map_var(i_local, j_local + 1, 'u')] = Q[2];
            system_matrix[rowbase, map_var(i_local - 1, j_local, 'u')] = Q[3];
            system_matrix[rowbase, map_var(i_local, j_local - 1, 'u')] = Q[4];
            //system_matrix[rowbase, map_var(i_local + 1, j_local + 1, 'u')] = Q[5];

            rhs[rowbase] -= p_boundary_conditions[global_j, BoundaryConditions.Direction.Positive_X] * P[1];
            rhs[rowbase] -= p_boundary_conditions[global_j + 1, BoundaryConditions.Direction.Positive_X] * P[5];

            rhs[rowbase] -= u_boundary_conditions[global_j, BoundaryConditions.Direction.Positive_X] * Q[1];
            rhs[rowbase] -= u_boundary_conditions[global_j + 1, BoundaryConditions.Direction.Positive_X] * Q[5];

            //y-momentum
            system_matrix[rowbase + 1, map_var(i_local, j_local, 'p')] = R[0];
            //system_matrix[rowbase + 1, map_var(i_local + 1, j_local, 'p')] = R[1];
            system_matrix[rowbase + 1, map_var(i_local, j_local + 1, 'p')] = R[2];
            system_matrix[rowbase + 1, map_var(i_local - 1, j_local, 'p')] = R[3];
            system_matrix[rowbase + 1, map_var(i_local, j_local - 1, 'p')] = R[4];
            //system_matrix[rowbase + 1, map_var(i_local + 1, j_local + 1, 'p')] = R[5];

            system_matrix[rowbase + 1, map_var(i_local, j_local, 'v')] = S[0];
            //system_matrix[rowbase + 1, map_var(i_local + 1, j_local, 'v')] = S[1];
            system_matrix[rowbase + 1, map_var(i_local, j_local + 1, 'v')] = S[2];
            system_matrix[rowbase + 1, map_var(i_local - 1, j_local, 'v')] = S[3];
            system_matrix[rowbase + 1, map_var(i_local, j_local - 1, 'v')] = S[4];
            //system_matrix[rowbase + 1, map_var(i_local + 1, j_local + 1, 'v')] = S[5];

            rhs[rowbase + 1] -= p_boundary_conditions[global_j, BoundaryConditions.Direction.Positive_X] * R[1];
            rhs[rowbase + 1] -= p_boundary_conditions[global_j + 1, BoundaryConditions.Direction.Positive_X] * R[5];

            rhs[rowbase + 1] -= v_boundary_conditions[global_j, BoundaryConditions.Direction.Positive_X] * R[1];
            rhs[rowbase + 1] -= v_boundary_conditions[global_j + 1, BoundaryConditions.Direction.Positive_X] * R[5];

            //continuity
            system_matrix[rowbase + 2, map_var(i_local, j_local, 'u')] = T[0];
            //system_matrix[rowbase + 2, map_var(i_local + 1, j_local, 'u')] = T[1];
            system_matrix[rowbase + 2, map_var(i_local, j_local + 1, 'u')] = T[2];
            system_matrix[rowbase + 2, map_var(i_local - 1, j_local, 'u')] = T[3];
            system_matrix[rowbase + 2, map_var(i_local, j_local - 1, 'u')] = T[4];
            //system_matrix[rowbase + 2, map_var(i_local + 1, j_local + 1, 'u')] = T[5];

            system_matrix[rowbase + 2, map_var(i_local, j_local, 'v')] = W[0];
            //system_matrix[rowbase + 2, map_var(i_local + 1, j_local, 'v')] = W[1];
            system_matrix[rowbase + 2, map_var(i_local, j_local + 1, 'v')] = W[2];
            system_matrix[rowbase + 2, map_var(i_local - 1, j_local, 'v')] = W[3];
            system_matrix[rowbase + 2, map_var(i_local, j_local - 1, 'v')] = W[4];
            //system_matrix[rowbase + 2, map_var(i_local + 1, j_local + 1, 'v')] = W[5];

            rhs[rowbase + 2] -= u_boundary_conditions[global_j, BoundaryConditions.Direction.Positive_X] * T[1];
            rhs[rowbase + 2] -= u_boundary_conditions[global_j + 1, BoundaryConditions.Direction.Positive_X] * T[5];

            rhs[rowbase + 2] -= v_boundary_conditions[global_j, BoundaryConditions.Direction.Positive_X] * W[1];
            rhs[rowbase + 2] -= v_boundary_conditions[global_j + 1, BoundaryConditions.Direction.Positive_X] * W[5];


        }
        private void map_upper_node_eqs(int rowbase, int i_local, int j_local, double[] P, double[] Q, double[] R, double[] S, double[] T, double[] W)
        {
            int global_i = i_local + 1;
            int global_j = j_local + 1;
            rhs[rowbase] = (-1 / deltat) * prev_u[i_local + 1, j_local + 1].Value;
            rhs[rowbase + 1] = (-1 / deltat) * prev_v[i_local + 1, j_local + 1].Value;
            rhs[rowbase + 2] = 0;
            //x-momentum
            system_matrix[rowbase, map_var(i_local, j_local, 'p')] = P[0];
            system_matrix[rowbase, map_var(i_local + 1, j_local, 'p')] = P[1];
            //system_matrix[rowbase, map_var(i_local, j_local + 1, 'p')] = P[2];//
            system_matrix[rowbase, map_var(i_local - 1, j_local, 'p')] = P[3];
            system_matrix[rowbase, map_var(i_local, j_local - 1, 'p')] = P[4];
            //system_matrix[rowbase, map_var(i_local + 1, j_local + 1, 'p')] = P[5];//

            system_matrix[rowbase, map_var(i_local, j_local, 'u')] = Q[0];
            system_matrix[rowbase, map_var(i_local + 1, j_local, 'u')] = Q[1];
            //system_matrix[rowbase, map_var(i_local, j_local + 1, 'u')] = Q[2];//
            system_matrix[rowbase, map_var(i_local - 1, j_local, 'u')] = Q[3];
            system_matrix[rowbase, map_var(i_local, j_local - 1, 'u')] = Q[4];
            //system_matrix[rowbase, map_var(i_local + 1, j_local + 1, 'u')] = Q[5];//

            rhs[rowbase] -= p_boundary_conditions[global_i, BoundaryConditions.Direction.Positive_Y] * P[2];
            rhs[rowbase] -= p_boundary_conditions[global_i + 1, BoundaryConditions.Direction.Positive_Y] * P[5];

            rhs[rowbase] -= u_boundary_conditions[global_i, BoundaryConditions.Direction.Positive_Y] * Q[2];
            rhs[rowbase] -= u_boundary_conditions[global_i + 1, BoundaryConditions.Direction.Positive_Y] * Q[5];

            //y-momentum
            system_matrix[rowbase + 1, map_var(i_local, j_local, 'p')] = R[0];
            system_matrix[rowbase + 1, map_var(i_local + 1, j_local, 'p')] = R[1];
            //system_matrix[rowbase + 1, map_var(i_local, j_local + 1, 'p')] = R[2];//
            system_matrix[rowbase + 1, map_var(i_local - 1, j_local, 'p')] = R[3];
            system_matrix[rowbase + 1, map_var(i_local, j_local - 1, 'p')] = R[4];
            //system_matrix[rowbase + 1, map_var(i_local + 1, j_local + 1, 'p')] = R[5];//

            system_matrix[rowbase + 1, map_var(i_local, j_local, 'v')] = S[0];
            system_matrix[rowbase + 1, map_var(i_local + 1, j_local, 'v')] = S[1];
            //system_matrix[rowbase + 1, map_var(i_local, j_local + 1, 'v')] = S[2];//
            system_matrix[rowbase + 1, map_var(i_local - 1, j_local, 'v')] = S[3];
            system_matrix[rowbase + 1, map_var(i_local, j_local - 1, 'v')] = S[4];
            //system_matrix[rowbase + 1, map_var(i_local + 1, j_local + 1, 'v')] = S[5];//

            rhs[rowbase + 1] -= p_boundary_conditions[global_i, BoundaryConditions.Direction.Positive_Y] * R[2];
            rhs[rowbase + 1] -= p_boundary_conditions[global_i + 1, BoundaryConditions.Direction.Positive_Y] * R[5];

            rhs[rowbase + 1] -= v_boundary_conditions[global_i, BoundaryConditions.Direction.Positive_Y] * R[2];
            rhs[rowbase + 1] -= v_boundary_conditions[global_i + 1, BoundaryConditions.Direction.Positive_Y] * R[5];

            //continuity
            system_matrix[rowbase + 2, map_var(i_local, j_local, 'u')] = T[0];
            system_matrix[rowbase + 2, map_var(i_local + 1, j_local, 'u')] = T[1];
            //system_matrix[rowbase + 2, map_var(i_local, j_local + 1, 'u')] = T[2];//
            system_matrix[rowbase + 2, map_var(i_local - 1, j_local, 'u')] = T[3];
            system_matrix[rowbase + 2, map_var(i_local, j_local - 1, 'u')] = T[4];
            //system_matrix[rowbase + 2, map_var(i_local + 1, j_local + 1, 'u')] = T[5];//

            system_matrix[rowbase + 2, map_var(i_local, j_local, 'v')] = W[0];
            system_matrix[rowbase + 2, map_var(i_local + 1, j_local, 'v')] = W[1];
            //system_matrix[rowbase + 2, map_var(i_local, j_local + 1, 'v')] = W[2];//
            system_matrix[rowbase + 2, map_var(i_local - 1, j_local, 'v')] = W[3];
            system_matrix[rowbase + 2, map_var(i_local, j_local - 1, 'v')] = W[4];
            //system_matrix[rowbase + 2, map_var(i_local + 1, j_local + 1, 'v')] = W[5];//

            rhs[rowbase + 2] -= u_boundary_conditions[global_i, BoundaryConditions.Direction.Positive_Y] * T[2];
            rhs[rowbase + 2] -= u_boundary_conditions[global_i + 1, BoundaryConditions.Direction.Positive_Y] * T[5];

            rhs[rowbase + 2] -= v_boundary_conditions[global_i, BoundaryConditions.Direction.Positive_Y] * W[2];
            rhs[rowbase + 2] -= v_boundary_conditions[global_i + 1, BoundaryConditions.Direction.Positive_Y] * W[5];


        }
        private void map_lower_node_eqs(int rowbase, int i_local, int j_local, double[] P, double[] Q, double[] R, double[] S, double[] T, double[] W)
        {
            int global_i = i_local + 1;
            int global_j = j_local + 1;
            rhs[rowbase] = (-1 / deltat) * prev_u[i_local + 1, j_local + 1].Value;
            rhs[rowbase + 1] = (-1 / deltat) * prev_v[i_local + 1, j_local + 1].Value;
            rhs[rowbase + 2] = 0;
            //x-momentum
            system_matrix[rowbase, map_var(i_local, j_local, 'p')] = P[0];
            system_matrix[rowbase, map_var(i_local + 1, j_local, 'p')] = P[1];
            system_matrix[rowbase, map_var(i_local, j_local + 1, 'p')] = P[2];
            system_matrix[rowbase, map_var(i_local - 1, j_local, 'p')] = P[3];
            //system_matrix[rowbase, map_var(i_local, j_local - 1, 'p')] = P[4];
            system_matrix[rowbase, map_var(i_local + 1, j_local + 1, 'p')] = P[5]; //problem here

            system_matrix[rowbase, map_var(i_local, j_local, 'u')] = Q[0];
            system_matrix[rowbase, map_var(i_local + 1, j_local, 'u')] = Q[1];
            system_matrix[rowbase, map_var(i_local, j_local + 1, 'u')] = Q[2];
            system_matrix[rowbase, map_var(i_local - 1, j_local, 'u')] = Q[3];
            //system_matrix[rowbase, map_var(i_local, j_local - 1, 'u')] = Q[4];
            system_matrix[rowbase, map_var(i_local + 1, j_local + 1, 'u')] = Q[5];

            rhs[rowbase] -= p_boundary_conditions[global_i, BoundaryConditions.Direction.Negative_Y] * P[4];

            rhs[rowbase] -= u_boundary_conditions[global_i, BoundaryConditions.Direction.Negative_Y] * Q[4];

            //y-momentum
            system_matrix[rowbase + 1, map_var(i_local, j_local, 'p')] = R[0];
            system_matrix[rowbase + 1, map_var(i_local + 1, j_local, 'p')] = R[1];
            system_matrix[rowbase + 1, map_var(i_local, j_local + 1, 'p')] = R[2];
            system_matrix[rowbase + 1, map_var(i_local - 1, j_local, 'p')] = R[3];
            //system_matrix[rowbase + 1, map_var(i_local, j_local - 1, 'p')] = R[4];
            system_matrix[rowbase + 1, map_var(i_local + 1, j_local + 1, 'p')] = R[5];

            system_matrix[rowbase + 1, map_var(i_local, j_local, 'v')] = S[0];
            system_matrix[rowbase + 1, map_var(i_local + 1, j_local, 'v')] = S[1];
            system_matrix[rowbase + 1, map_var(i_local, j_local + 1, 'v')] = S[2];
            system_matrix[rowbase + 1, map_var(i_local - 1, j_local, 'v')] = S[3];
            //system_matrix[rowbase + 1, map_var(i_local, j_local - 1, 'v')] = S[4];
            system_matrix[rowbase + 1, map_var(i_local + 1, j_local + 1, 'v')] = S[5];

            rhs[rowbase+1] -= p_boundary_conditions[global_i, BoundaryConditions.Direction.Negative_Y] * R[4];

            rhs[rowbase+1] -= v_boundary_conditions[global_i, BoundaryConditions.Direction.Negative_Y] * S[4];

            //continuity
            system_matrix[rowbase + 2, map_var(i_local, j_local, 'u')] = T[0];
            system_matrix[rowbase + 2, map_var(i_local + 1, j_local, 'u')] = T[1];
            system_matrix[rowbase + 2, map_var(i_local, j_local + 1, 'u')] = T[2];
            system_matrix[rowbase + 2, map_var(i_local - 1, j_local, 'u')] = T[3];
            //system_matrix[rowbase + 2, map_var(i_local, j_local - 1, 'u')] = T[4];
            system_matrix[rowbase + 2, map_var(i_local + 1, j_local + 1, 'u')] = T[5];

            system_matrix[rowbase + 2, map_var(i_local, j_local, 'v')] = W[0];
            system_matrix[rowbase + 2, map_var(i_local + 1, j_local, 'v')] = W[1];
            system_matrix[rowbase + 2, map_var(i_local, j_local + 1, 'v')] = W[2];
            system_matrix[rowbase + 2, map_var(i_local - 1, j_local, 'v')] = W[3];
            //system_matrix[rowbase + 2, map_var(i_local, j_local - 1, 'v')] = W[4];
            system_matrix[rowbase + 2, map_var(i_local + 1, j_local + 1, 'v')] = W[5];

            rhs[rowbase+2] -= u_boundary_conditions[global_i, BoundaryConditions.Direction.Negative_Y] * T[4];

            rhs[rowbase+2] -= v_boundary_conditions[global_i, BoundaryConditions.Direction.Negative_Y] * W[4];


        }
        private void map_left_node_eqs(int rowbase, int i_local, int j_local, double[] P, double[] Q, double[] R, double[] S, double[] T, double[] W)
        {
            int global_i = i_local + 1;
            int global_j = j_local + 1;
            rhs[rowbase] = (-1 / deltat) * prev_u[i_local + 1, j_local + 1].Value;
            rhs[rowbase + 1] = (-1 / deltat) * prev_v[i_local + 1, j_local + 1].Value;
            rhs[rowbase + 2] = 0;
            //x-momentum
            system_matrix[rowbase, map_var(i_local, j_local, 'p')] = P[0];
            system_matrix[rowbase, map_var(i_local + 1, j_local, 'p')] = P[1];
            system_matrix[rowbase, map_var(i_local, j_local + 1, 'p')] = P[2];
            //system_matrix[rowbase, map_var(i_local - 1, j_local, 'p')] = P[3];
            system_matrix[rowbase, map_var(i_local, j_local - 1, 'p')] = P[4];
            system_matrix[rowbase, map_var(i_local + 1, j_local + 1, 'p')] = P[5];

            system_matrix[rowbase, map_var(i_local, j_local, 'u')] = Q[0];
            system_matrix[rowbase, map_var(i_local + 1, j_local, 'u')] = Q[1];
            system_matrix[rowbase, map_var(i_local, j_local + 1, 'u')] = Q[2];
            //system_matrix[rowbase, map_var(i_local - 1, j_local, 'u')] = Q[3];
            system_matrix[rowbase, map_var(i_local, j_local - 1, 'u')] = Q[4];
            system_matrix[rowbase, map_var(i_local + 1, j_local + 1, 'u')] = Q[5];

            rhs[rowbase] -= p_boundary_conditions[global_j, BoundaryConditions.Direction.Negative_X] * P[3];

            rhs[rowbase] -= u_boundary_conditions[global_j, BoundaryConditions.Direction.Negative_X] * Q[3];

            //y-momentum
            system_matrix[rowbase + 1, map_var(i_local, j_local, 'p')] = R[0];
            system_matrix[rowbase + 1, map_var(i_local + 1, j_local, 'p')] = R[1];
            system_matrix[rowbase + 1, map_var(i_local, j_local + 1, 'p')] = R[2];
            //system_matrix[rowbase + 1, map_var(i_local - 1, j_local, 'p')] = R[3];
            system_matrix[rowbase + 1, map_var(i_local, j_local - 1, 'p')] = R[4];
            system_matrix[rowbase + 1, map_var(i_local + 1, j_local + 1, 'p')] = R[5];

            system_matrix[rowbase + 1, map_var(i_local, j_local, 'v')] = S[0];
            system_matrix[rowbase + 1, map_var(i_local + 1, j_local, 'v')] = S[1];
            system_matrix[rowbase + 1, map_var(i_local, j_local + 1, 'v')] = S[2];
            //system_matrix[rowbase + 1, map_var(i_local - 1, j_local, 'v')] = S[3];
            system_matrix[rowbase + 1, map_var(i_local, j_local - 1, 'v')] = S[4];
            system_matrix[rowbase + 1, map_var(i_local + 1, j_local + 1, 'v')] = S[5];

            rhs[rowbase + 1] -= p_boundary_conditions[global_j, BoundaryConditions.Direction.Negative_X] * R[3];

            rhs[rowbase + 1] -= v_boundary_conditions[global_j, BoundaryConditions.Direction.Negative_X] * S[3];

            //continuity
            system_matrix[rowbase + 2, map_var(i_local, j_local, 'u')] = T[0];
            system_matrix[rowbase + 2, map_var(i_local + 1, j_local, 'u')] = T[1];
            system_matrix[rowbase + 2, map_var(i_local, j_local + 1, 'u')] = T[2];
            //system_matrix[rowbase + 2, map_var(i_local - 1, j_local, 'u')] = T[3];
            system_matrix[rowbase + 2, map_var(i_local, j_local - 1, 'u')] = T[4];
            system_matrix[rowbase + 2, map_var(i_local + 1, j_local + 1, 'u')] = T[5];

            system_matrix[rowbase + 2, map_var(i_local, j_local, 'v')] = W[0];
            system_matrix[rowbase + 2, map_var(i_local + 1, j_local, 'v')] = W[1];
            system_matrix[rowbase + 2, map_var(i_local, j_local + 1, 'v')] = W[2];
            //system_matrix[rowbase + 2, map_var(i_local - 1, j_local, 'v')] = W[3];
            system_matrix[rowbase + 2, map_var(i_local, j_local - 1, 'v')] = W[4];
            system_matrix[rowbase + 2, map_var(i_local + 1, j_local + 1, 'v')] = W[5];

            rhs[rowbase + 2] -= u_boundary_conditions[global_j, BoundaryConditions.Direction.Negative_X] * T[3];

            rhs[rowbase + 2] -= v_boundary_conditions[global_j, BoundaryConditions.Direction.Negative_X] * W[3];
        }

        private void map_ulcorner_node_eqs(int rowbase, int i_local, int j_local, double[] P, double[] Q, double[] R, double[] S, double[] T, double[] W)
        {
            int global_i = i_local + 1;
            int global_j = j_local + 1;
            rhs[rowbase] = (-1 / deltat) * prev_u[i_local + 1, j_local + 1].Value;
            rhs[rowbase + 1] = (-1 / deltat) * prev_v[i_local + 1, j_local + 1].Value;
            rhs[rowbase + 2] = 0;
            //x-momentum
            system_matrix[rowbase, map_var(i_local, j_local, 'p')] = P[0];
            system_matrix[rowbase, map_var(i_local + 1, j_local, 'p')] = P[1];
            //system_matrix[rowbase, map_var(i_local, j_local + 1, 'p')] = P[2];
            //system_matrix[rowbase, map_var(i_local - 1, j_local, 'p')] = P[3];
            system_matrix[rowbase, map_var(i_local, j_local - 1, 'p')] = P[4];
            //system_matrix[rowbase, map_var(i_local + 1, j_local + 1, 'p')] = P[5];

            system_matrix[rowbase, map_var(i_local, j_local, 'u')] = Q[0];
            system_matrix[rowbase, map_var(i_local + 1, j_local, 'u')] = Q[1];
            //system_matrix[rowbase, map_var(i_local, j_local + 1, 'u')] = Q[2];
            //system_matrix[rowbase, map_var(i_local - 1, j_local, 'u')] = Q[3];
            system_matrix[rowbase, map_var(i_local, j_local - 1, 'u')] = Q[4];
            //system_matrix[rowbase, map_var(i_local + 1, j_local + 1, 'u')] = Q[5];

            rhs[rowbase] -= p_boundary_conditions[global_i, BoundaryConditions.Direction.Positive_Y] * P[1];
            rhs[rowbase] -= p_boundary_conditions[global_i + 1, BoundaryConditions.Direction.Positive_Y] * P[5];
            rhs[rowbase] -= p_boundary_conditions[global_j, BoundaryConditions.Direction.Negative_X] * P[3];

            rhs[rowbase] -= u_boundary_conditions[global_i, BoundaryConditions.Direction.Positive_Y] * Q[1];
            rhs[rowbase] -= u_boundary_conditions[global_i + 1, BoundaryConditions.Direction.Positive_Y] * Q[5];
            rhs[rowbase] -= u_boundary_conditions[global_j, BoundaryConditions.Direction.Negative_X] *Q[3];

            //y-momentum
            system_matrix[rowbase + 1, map_var(i_local, j_local, 'p')] = R[0];
            system_matrix[rowbase + 1, map_var(i_local + 1, j_local, 'p')] = R[1];
            //system_matrix[rowbase + 1, map_var(i_local, j_local + 1, 'p')] = R[2];
            //system_matrix[rowbase + 1, map_var(i_local - 1, j_local, 'p')] = R[3];
            system_matrix[rowbase + 1, map_var(i_local, j_local - 1, 'p')] = R[4];
            //system_matrix[rowbase + 1, map_var(i_local + 1, j_local + 1, 'p')] = R[5];

            system_matrix[rowbase + 1, map_var(i_local, j_local, 'v')] = S[0];
            system_matrix[rowbase + 1, map_var(i_local + 1, j_local, 'v')] = S[1];
            //system_matrix[rowbase + 1, map_var(i_local, j_local + 1, 'v')] = S[2];
            //system_matrix[rowbase + 1, map_var(i_local - 1, j_local, 'v')] = S[3];
            system_matrix[rowbase + 1, map_var(i_local, j_local - 1, 'v')] = S[4];
            //system_matrix[rowbase + 1, map_var(i_local + 1, j_local + 1, 'v')] = S[5];

            rhs[rowbase+1] -= p_boundary_conditions[global_i, BoundaryConditions.Direction.Positive_Y] * R[1];
            rhs[rowbase+1] -= p_boundary_conditions[global_i + 1, BoundaryConditions.Direction.Positive_Y] * R[5];
            rhs[rowbase+1] -= p_boundary_conditions[global_j, BoundaryConditions.Direction.Negative_X] * R[3];
                       
            rhs[rowbase+1] -= v_boundary_conditions[global_i, BoundaryConditions.Direction.Positive_Y] * S[1];
            rhs[rowbase+1] -= v_boundary_conditions[global_i + 1, BoundaryConditions.Direction.Positive_Y] * S[5];
            rhs[rowbase+1] -= v_boundary_conditions[global_j, BoundaryConditions.Direction.Negative_X] * S[3];

            //continuity
            system_matrix[rowbase + 2, map_var(i_local, j_local, 'u')] = T[0];
            system_matrix[rowbase + 2, map_var(i_local + 1, j_local, 'u')] = T[1];
            //system_matrix[rowbase + 2, map_var(i_local, j_local + 1, 'u')] = T[2];
            //system_matrix[rowbase + 2, map_var(i_local - 1, j_local, 'u')] = T[3];
            system_matrix[rowbase + 2, map_var(i_local, j_local - 1, 'u')] = T[4];
            //system_matrix[rowbase + 2, map_var(i_local + 1, j_local + 1, 'u')] = T[5];

            system_matrix[rowbase + 2, map_var(i_local, j_local, 'v')] = W[0];
            system_matrix[rowbase + 2, map_var(i_local + 1, j_local, 'v')] = W[1];
            //system_matrix[rowbase + 2, map_var(i_local, j_local + 1, 'v')] = W[2];
            //system_matrix[rowbase + 2, map_var(i_local - 1, j_local, 'v')] = W[3];
            system_matrix[rowbase + 2, map_var(i_local, j_local - 1, 'v')] = W[4];
            //system_matrix[rowbase + 2, map_var(i_local + 1, j_local + 1, 'v')] = W[5];

            rhs[rowbase+2] -= u_boundary_conditions[global_i, BoundaryConditions.Direction.Positive_Y] * T[1];
            rhs[rowbase+2] -= u_boundary_conditions[global_i + 1, BoundaryConditions.Direction.Positive_Y] * T[5];
            rhs[rowbase+2] -= u_boundary_conditions[global_j, BoundaryConditions.Direction.Negative_X] * T[3];

            rhs[rowbase+2] -= v_boundary_conditions[global_i, BoundaryConditions.Direction.Positive_Y] * W[1];
            rhs[rowbase+2] -= v_boundary_conditions[global_i + 1, BoundaryConditions.Direction.Positive_Y] * W[5];
            rhs[rowbase+2] -= v_boundary_conditions[global_j, BoundaryConditions.Direction.Negative_X] * W[3];


        }
        private void map_urcorner_node_eqs(int rowbase, int i_local, int j_local, double[] P, double[] Q, double[] R, double[] S, double[] T, double[] W)
        {
            int global_i = i_local + 1;
            int global_j = j_local + 1;
            rhs[rowbase] = (-1 / deltat) * prev_u[i_local + 1, j_local + 1].Value;
            rhs[rowbase + 1] = (-1 / deltat) * prev_v[i_local + 1, j_local + 1].Value;
            rhs[rowbase + 2] = 0;
            //x-momentum
            system_matrix[rowbase, map_var(i_local, j_local, 'p')] = P[0];
            //system_matrix[rowbase, map_var(i_local + 1, j_local, 'p')] = P[1];
            //system_matrix[rowbase, map_var(i_local, j_local + 1, 'p')] = P[2];
            system_matrix[rowbase, map_var(i_local - 1, j_local, 'p')] = P[3];
            system_matrix[rowbase, map_var(i_local, j_local - 1, 'p')] = P[4];
            //system_matrix[rowbase, map_var(i_local + 1, j_local + 1, 'p')] = P[5];

            system_matrix[rowbase, map_var(i_local, j_local, 'u')] = Q[0];
            //system_matrix[rowbase, map_var(i_local + 1, j_local, 'u')] = Q[1];
            //system_matrix[rowbase, map_var(i_local, j_local + 1, 'u')] = Q[2];
            system_matrix[rowbase, map_var(i_local - 1, j_local, 'u')] = Q[3];
            system_matrix[rowbase, map_var(i_local, j_local - 1, 'u')] = Q[4];
            //system_matrix[rowbase, map_var(i_local + 1, j_local + 1, 'u')] = Q[5];

            rhs[rowbase] -= p_boundary_conditions[global_i, BoundaryConditions.Direction.Positive_Y] * P[1];
            rhs[rowbase] -= p_boundary_conditions[global_i + 1, BoundaryConditions.Direction.Positive_Y] * P[5];
            rhs[rowbase] -= p_boundary_conditions[global_j, BoundaryConditions.Direction.Positive_X] * P[2];

            rhs[rowbase] -= u_boundary_conditions[global_i, BoundaryConditions.Direction.Positive_Y] * Q[1];
            rhs[rowbase] -= u_boundary_conditions[global_i + 1, BoundaryConditions.Direction.Positive_Y] * Q[5];
            rhs[rowbase] -= u_boundary_conditions[global_j, BoundaryConditions.Direction.Positive_X] * Q[2];

            //y-momentum
            system_matrix[rowbase + 1, map_var(i_local, j_local, 'p')] = R[0];
            //system_matrix[rowbase + 1, map_var(i_local + 1, j_local, 'p')] = R[1];
            //system_matrix[rowbase + 1, map_var(i_local, j_local + 1, 'p')] = R[2];
            system_matrix[rowbase + 1, map_var(i_local - 1, j_local, 'p')] = R[3];
            system_matrix[rowbase + 1, map_var(i_local, j_local - 1, 'p')] = R[4];
            //system_matrix[rowbase + 1, map_var(i_local + 1, j_local + 1, 'p')] = R[5];

            system_matrix[rowbase + 1, map_var(i_local, j_local, 'v')] = S[0];
            //system_matrix[rowbase + 1, map_var(i_local + 1, j_local, 'v')] = S[1];
            //system_matrix[rowbase + 1, map_var(i_local, j_local + 1, 'v')] = S[2];
            system_matrix[rowbase + 1, map_var(i_local - 1, j_local, 'v')] = S[3];
            system_matrix[rowbase + 1, map_var(i_local, j_local - 1, 'v')] = S[4];
            //system_matrix[rowbase + 1, map_var(i_local + 1, j_local + 1, 'v')] = S[5];

            rhs[rowbase + 1] -= p_boundary_conditions[global_i, BoundaryConditions.Direction.Positive_Y] * R[1];
            rhs[rowbase + 1] -= p_boundary_conditions[global_i + 1, BoundaryConditions.Direction.Positive_Y] * R[5];
            rhs[rowbase + 1] -= p_boundary_conditions[global_j, BoundaryConditions.Direction.Positive_X] * R[2];

            rhs[rowbase + 1] -= v_boundary_conditions[global_i, BoundaryConditions.Direction.Positive_Y] * S[1];
            rhs[rowbase + 1] -= v_boundary_conditions[global_i + 1, BoundaryConditions.Direction.Positive_Y] * S[5];
            rhs[rowbase + 1] -= v_boundary_conditions[global_j, BoundaryConditions.Direction.Positive_X] * S[2];

            //continuity
            system_matrix[rowbase + 2, map_var(i_local, j_local, 'u')] = T[0];
            //system_matrix[rowbase + 2, map_var(i_local + 1, j_local, 'u')] = T[1];
            //system_matrix[rowbase + 2, map_var(i_local, j_local + 1, 'u')] = T[2];
            system_matrix[rowbase + 2, map_var(i_local - 1, j_local, 'u')] = T[3];
            system_matrix[rowbase + 2, map_var(i_local, j_local - 1, 'u')] = T[4];
            //system_matrix[rowbase + 2, map_var(i_local + 1, j_local + 1, 'u')] = T[5];

            system_matrix[rowbase + 2, map_var(i_local, j_local, 'v')] = W[0];
            //system_matrix[rowbase + 2, map_var(i_local + 1, j_local, 'v')] = W[1];
            //system_matrix[rowbase + 2, map_var(i_local, j_local + 1, 'v')] = W[2];
            system_matrix[rowbase + 2, map_var(i_local - 1, j_local, 'v')] = W[3];
            system_matrix[rowbase + 2, map_var(i_local, j_local - 1, 'v')] = W[4];
            //system_matrix[rowbase + 2, map_var(i_local + 1, j_local + 1, 'v')] = W[5];

            rhs[rowbase + 2] -= u_boundary_conditions[global_i, BoundaryConditions.Direction.Positive_Y] * T[1];
            rhs[rowbase + 2] -= u_boundary_conditions[global_i + 1, BoundaryConditions.Direction.Positive_Y] * T[5];
            rhs[rowbase + 2] -= u_boundary_conditions[global_j, BoundaryConditions.Direction.Positive_X] * T[2];

            rhs[rowbase + 2] -= v_boundary_conditions[global_i, BoundaryConditions.Direction.Positive_Y] * W[1];
            rhs[rowbase + 2] -= v_boundary_conditions[global_i + 1, BoundaryConditions.Direction.Positive_Y] * W[5];
            rhs[rowbase + 2] -= v_boundary_conditions[global_j, BoundaryConditions.Direction.Positive_X] * W[2];


        }
        private void map_llcorner_node_eqs(int rowbase, int i_local, int j_local, double[] P, double[] Q, double[] R, double[] S, double[] T, double[] W)
        {
            int global_i = i_local + 1;
            int global_j = j_local + 1;
            rhs[rowbase] = (-1 / deltat) * prev_u[i_local + 1, j_local + 1].Value;
            rhs[rowbase + 1] = (-1 / deltat) * prev_v[i_local + 1, j_local + 1].Value;
            rhs[rowbase + 2] = 0;
            //x-momentum
            system_matrix[rowbase, map_var(i_local, j_local, 'p')] = P[0];
            system_matrix[rowbase, map_var(i_local + 1, j_local, 'p')] = P[1];
            system_matrix[rowbase, map_var(i_local, j_local + 1, 'p')] = P[2];
            //system_matrix[rowbase, map_var(i_local - 1, j_local, 'p')] = P[3];
            //system_matrix[rowbase, map_var(i_local, j_local - 1, 'p')] = P[4];
            system_matrix[rowbase, map_var(i_local + 1, j_local + 1, 'p')] = P[5];

            system_matrix[rowbase, map_var(i_local, j_local, 'u')] = Q[0];
            system_matrix[rowbase, map_var(i_local + 1, j_local, 'u')] = Q[1];
            system_matrix[rowbase, map_var(i_local, j_local + 1, 'u')] = Q[2];
            //system_matrix[rowbase, map_var(i_local - 1, j_local, 'u')] = Q[3];
            //system_matrix[rowbase, map_var(i_local, j_local - 1, 'u')] = Q[4];
            system_matrix[rowbase, map_var(i_local + 1, j_local + 1, 'u')] = Q[5];

            rhs[rowbase] -= p_boundary_conditions[global_i, BoundaryConditions.Direction.Negative_Y] * P[4];
            rhs[rowbase] -= p_boundary_conditions[global_j, BoundaryConditions.Direction.Negative_X] * P[3];

            rhs[rowbase] -= u_boundary_conditions[global_i, BoundaryConditions.Direction.Negative_Y] * Q[4];
            rhs[rowbase] -= u_boundary_conditions[global_j, BoundaryConditions.Direction.Negative_X] * Q[3];

            //y-momentum
            system_matrix[rowbase + 1, map_var(i_local, j_local, 'p')] = R[0];
            system_matrix[rowbase + 1, map_var(i_local + 1, j_local, 'p')] = R[1];
            system_matrix[rowbase + 1, map_var(i_local, j_local + 1, 'p')] = R[2];
            //system_matrix[rowbase + 1, map_var(i_local - 1, j_local, 'p')] = R[3];
            //system_matrix[rowbase + 1, map_var(i_local, j_local - 1, 'p')] = R[4];
            system_matrix[rowbase + 1, map_var(i_local + 1, j_local + 1, 'p')] = R[5];

            system_matrix[rowbase + 1, map_var(i_local, j_local, 'v')] = S[0];
            system_matrix[rowbase + 1, map_var(i_local + 1, j_local, 'v')] = S[1];
            system_matrix[rowbase + 1, map_var(i_local, j_local + 1, 'v')] = S[2];
            //system_matrix[rowbase + 1, map_var(i_local - 1, j_local, 'v')] = S[3];
            //system_matrix[rowbase + 1, map_var(i_local, j_local - 1, 'v')] = S[4];
            system_matrix[rowbase + 1, map_var(i_local + 1, j_local + 1, 'v')] = S[5];

            rhs[rowbase + 1] -= p_boundary_conditions[global_i, BoundaryConditions.Direction.Negative_Y] * R[4];
            rhs[rowbase + 1] -= p_boundary_conditions[global_j, BoundaryConditions.Direction.Negative_X] * R[3];

            rhs[rowbase + 1] -= v_boundary_conditions[global_i, BoundaryConditions.Direction.Negative_Y] * S[4];
            rhs[rowbase + 1] -= v_boundary_conditions[global_j, BoundaryConditions.Direction.Negative_X] * S[3];

            //continuity
            system_matrix[rowbase + 2, map_var(i_local, j_local, 'u')] = T[0];
            system_matrix[rowbase + 2, map_var(i_local + 1, j_local, 'u')] = T[1];
            system_matrix[rowbase + 2, map_var(i_local, j_local + 1, 'u')] = T[2];
            //system_matrix[rowbase + 2, map_var(i_local - 1, j_local, 'u')] = T[3];
            //system_matrix[rowbase + 2, map_var(i_local, j_local - 1, 'u')] = T[4];
            system_matrix[rowbase + 2, map_var(i_local + 1, j_local + 1, 'u')] = T[5];

            system_matrix[rowbase + 2, map_var(i_local, j_local, 'v')] = W[0];
            system_matrix[rowbase + 2, map_var(i_local + 1, j_local, 'v')] = W[1];
            system_matrix[rowbase + 2, map_var(i_local, j_local + 1, 'v')] = W[2];
            //system_matrix[rowbase + 2, map_var(i_local - 1, j_local, 'v')] = W[3];
            //system_matrix[rowbase + 2, map_var(i_local, j_local - 1, 'v')] = W[4];
            system_matrix[rowbase + 2, map_var(i_local + 1, j_local + 1, 'v')] = W[5];

            rhs[rowbase + 2] -= u_boundary_conditions[global_i, BoundaryConditions.Direction.Negative_Y] * T[4];
            rhs[rowbase + 2] -= u_boundary_conditions[global_j, BoundaryConditions.Direction.Negative_X] * T[3];

            rhs[rowbase + 2] -= v_boundary_conditions[global_i, BoundaryConditions.Direction.Negative_Y] * W[4];
            rhs[rowbase + 2] -= v_boundary_conditions[global_j, BoundaryConditions.Direction.Negative_X] * W[3];


        }
        private void map_lrcorner_node_eqs(int rowbase, int i_local, int j_local, double[] P, double[] Q, double[] R, double[] S, double[] T, double[] W)
        {
            int global_i = i_local + 1;
            int global_j = j_local + 1;
            rhs[rowbase] = (-1 / deltat) * prev_u[i_local + 1, j_local + 1].Value;
            rhs[rowbase + 1] = (-1 / deltat) * prev_v[i_local + 1, j_local + 1].Value;
            rhs[rowbase + 2] = 0;
            //x-momentum
            system_matrix[rowbase, map_var(i_local, j_local, 'p')] = P[0];
            //system_matrix[rowbase, map_var(i_local + 1, j_local, 'p')] = P[1];
            system_matrix[rowbase, map_var(i_local, j_local + 1, 'p')] = P[2];
            system_matrix[rowbase, map_var(i_local - 1, j_local, 'p')] = P[3];
            //system_matrix[rowbase, map_var(i_local, j_local - 1, 'p')] = P[4];
            //system_matrix[rowbase, map_var(i_local + 1, j_local + 1, 'p')] = P[5];

            system_matrix[rowbase, map_var(i_local, j_local, 'u')] = Q[0];
            //system_matrix[rowbase, map_var(i_local + 1, j_local, 'u')] = Q[1];
            system_matrix[rowbase, map_var(i_local, j_local + 1, 'u')] = Q[2];
            system_matrix[rowbase, map_var(i_local - 1, j_local, 'u')] = Q[3];
            //system_matrix[rowbase, map_var(i_local, j_local - 1, 'u')] = Q[4];
            //system_matrix[rowbase, map_var(i_local + 1, j_local + 1, 'u')] = Q[5];

            rhs[rowbase] -= p_boundary_conditions[global_i, BoundaryConditions.Direction.Negative_Y] * P[4];
            rhs[rowbase] -= p_boundary_conditions[global_j, BoundaryConditions.Direction.Positive_X] * P[1];
            rhs[rowbase] -= p_boundary_conditions[global_j+1, BoundaryConditions.Direction.Positive_X] * P[5];

            rhs[rowbase] -= u_boundary_conditions[global_i, BoundaryConditions.Direction.Negative_Y] * Q[4];
            rhs[rowbase] -= u_boundary_conditions[global_j, BoundaryConditions.Direction.Positive_X] * Q[1];
            rhs[rowbase] -= u_boundary_conditions[global_j+1, BoundaryConditions.Direction.Positive_X] * Q[5];

            //y-momentum
            system_matrix[rowbase + 1, map_var(i_local, j_local, 'p')] = R[0];
            //system_matrix[rowbase + 1, map_var(i_local + 1, j_local, 'p')] = R[1];
            system_matrix[rowbase + 1, map_var(i_local, j_local + 1, 'p')] = R[2];
            system_matrix[rowbase + 1, map_var(i_local - 1, j_local, 'p')] = R[3];
            //system_matrix[rowbase + 1, map_var(i_local, j_local - 1, 'p')] = R[4];
            //system_matrix[rowbase + 1, map_var(i_local + 1, j_local + 1, 'p')] = R[5];

            system_matrix[rowbase + 1, map_var(i_local, j_local, 'v')] = S[0];
            //system_matrix[rowbase + 1, map_var(i_local + 1, j_local, 'v')] = S[1];
            system_matrix[rowbase + 1, map_var(i_local, j_local + 1, 'v')] = S[2];
            system_matrix[rowbase + 1, map_var(i_local - 1, j_local, 'v')] = S[3];
            //system_matrix[rowbase + 1, map_var(i_local, j_local - 1, 'v')] = S[4];
            //system_matrix[rowbase + 1, map_var(i_local + 1, j_local + 1, 'v')] = S[5];

            rhs[rowbase + 1] -= p_boundary_conditions[global_i, BoundaryConditions.Direction.Negative_Y] * R[4];
            rhs[rowbase + 1] -= p_boundary_conditions[global_j, BoundaryConditions.Direction.Positive_X] * R[1];
            rhs[rowbase + 1] -= p_boundary_conditions[global_j+1, BoundaryConditions.Direction.Positive_X] * R[5];

            rhs[rowbase + 1] -= v_boundary_conditions[global_i, BoundaryConditions.Direction.Negative_Y] * S[4];
            rhs[rowbase + 1] -= v_boundary_conditions[global_j, BoundaryConditions.Direction.Positive_X] * S[1];
            rhs[rowbase + 1] -= v_boundary_conditions[global_j+1, BoundaryConditions.Direction.Positive_X] * S[5];

            //continuity
            system_matrix[rowbase + 2, map_var(i_local, j_local, 'u')] = T[0];
            //system_matrix[rowbase + 2, map_var(i_local + 1, j_local, 'u')] = T[1];
            system_matrix[rowbase + 2, map_var(i_local, j_local + 1, 'u')] = T[2];
            system_matrix[rowbase + 2, map_var(i_local - 1, j_local, 'u')] = T[3];
            //system_matrix[rowbase + 2, map_var(i_local, j_local - 1, 'u')] = T[4];
            //system_matrix[rowbase + 2, map_var(i_local + 1, j_local + 1, 'u')] = T[5];

            system_matrix[rowbase + 2, map_var(i_local, j_local, 'v')] = W[0];
            //system_matrix[rowbase + 2, map_var(i_local + 1, j_local, 'v')] = W[1];
            system_matrix[rowbase + 2, map_var(i_local, j_local + 1, 'v')] = W[2];
            system_matrix[rowbase + 2, map_var(i_local - 1, j_local, 'v')] = W[3];
            //system_matrix[rowbase + 2, map_var(i_local, j_local - 1, 'v')] = W[4];
            //system_matrix[rowbase + 2, map_var(i_local + 1, j_local + 1, 'v')] = W[5];

            rhs[rowbase + 2] -= u_boundary_conditions[global_i, BoundaryConditions.Direction.Negative_Y] * T[4];
            rhs[rowbase + 2] -= u_boundary_conditions[global_j, BoundaryConditions.Direction.Positive_X] * T[1];
            rhs[rowbase + 2] -= u_boundary_conditions[global_j+1, BoundaryConditions.Direction.Positive_X] * T[5];

            rhs[rowbase + 2] -= v_boundary_conditions[global_i, BoundaryConditions.Direction.Negative_Y] * W[4];
            rhs[rowbase + 2] -= v_boundary_conditions[global_j, BoundaryConditions.Direction.Positive_X] * W[1];
            rhs[rowbase + 2] -= v_boundary_conditions[global_j+1, BoundaryConditions.Direction.Positive_X] * W[5];


        }

        private void map_interior_node_eqs(int rowbase, int i_local, int j_local, double[] P, double[] Q, double [] R, double[] S, double[] T, double[] W)
        {
                //x-momentum
                system_matrix[rowbase, map_var(i_local, j_local, 'p')] = P[0];
                system_matrix[rowbase, map_var(i_local + 1, j_local, 'p')] = P[1];
                system_matrix[rowbase, map_var(i_local, j_local + 1, 'p')] = P[2];
                system_matrix[rowbase, map_var(i_local - 1, j_local, 'p')] = P[3];
                system_matrix[rowbase, map_var(i_local, j_local - 1, 'p')] = P[4];
                system_matrix[rowbase, map_var(i_local + 1, j_local + 1, 'p')] = P[5];

                system_matrix[rowbase, map_var(i_local, j_local, 'u')] = Q[0];
                system_matrix[rowbase, map_var(i_local + 1, j_local, 'u')] = Q[1];
                system_matrix[rowbase, map_var(i_local, j_local + 1, 'u')] = Q[2];
                system_matrix[rowbase, map_var(i_local - 1, j_local, 'u')] = Q[3];
                system_matrix[rowbase, map_var(i_local, j_local - 1, 'u')] = Q[4];
                system_matrix[rowbase, map_var(i_local + 1, j_local + 1, 'u')] = Q[5];

                //y-momentum
                system_matrix[rowbase + 1, map_var(i_local, j_local, 'p')] = R[0];
                system_matrix[rowbase + 1, map_var(i_local + 1, j_local, 'p')] = R[1];
                system_matrix[rowbase + 1, map_var(i_local, j_local + 1, 'p')] = R[2];
                system_matrix[rowbase + 1, map_var(i_local - 1, j_local, 'p')] = R[3];
                system_matrix[rowbase + 1, map_var(i_local, j_local - 1, 'p')] = R[4];
                system_matrix[rowbase + 1, map_var(i_local + 1, j_local + 1, 'p')] = R[5];

                system_matrix[rowbase + 1, map_var(i_local, j_local, 'v')] = S[0];
                system_matrix[rowbase + 1, map_var(i_local + 1, j_local, 'v')] = S[1];
                system_matrix[rowbase + 1, map_var(i_local, j_local + 1, 'v')] = S[2];
                system_matrix[rowbase + 1, map_var(i_local - 1, j_local, 'v')] = S[3];
                system_matrix[rowbase + 1, map_var(i_local, j_local - 1, 'v')] = S[4];
                system_matrix[rowbase + 1, map_var(i_local + 1, j_local + 1, 'v')] = S[5];

                //continuity
                system_matrix[rowbase + 2, map_var(i_local, j_local, 'u')] = T[0];
                system_matrix[rowbase + 2, map_var(i_local + 1, j_local, 'u')] = T[1];
                system_matrix[rowbase + 2, map_var(i_local, j_local + 1, 'u')] = T[2];
                system_matrix[rowbase + 2, map_var(i_local - 1, j_local, 'u')] = T[3];
                system_matrix[rowbase + 2, map_var(i_local, j_local - 1, 'u')] = T[4];
                system_matrix[rowbase + 2, map_var(i_local + 1, j_local + 1, 'u')] = T[5];

                system_matrix[rowbase + 2, map_var(i_local, j_local, 'v')] = W[0];
                system_matrix[rowbase + 2, map_var(i_local + 1, j_local, 'v')] = W[1];
                system_matrix[rowbase + 2, map_var(i_local, j_local + 1, 'v')] = W[2];
                system_matrix[rowbase + 2, map_var(i_local - 1, j_local, 'v')] = W[3];
                system_matrix[rowbase + 2, map_var(i_local, j_local - 1, 'v')] = W[4];
                system_matrix[rowbase + 2, map_var(i_local + 1, j_local + 1, 'v')] = W[5];

            rhs[rowbase] = (-1 / deltat) * prev_u[i_local+1, j_local+1].Value;
            rhs[rowbase+1] = (-1 / deltat) * prev_v[i_local+1, j_local+1].Value;
            rhs[rowbase+2] = 0;
        }

        private int map_var(int local_i, int local_j, char symbol)
        {
            int local_offset = 0;
            switch (symbol)
            {
                case 'p':
                    {
                        local_offset = 0;
                        break;
                    }
                case 'u':
                    {
                        local_offset = 1;
                        break;
                    }
                case 'v':
                    {
                        local_offset = 2;
                        break;
                    }
            }
            int block = (local_j) + (eq_nodes_y * (local_i));
            int y =  (3 * block) + local_offset;
            return y;
        }
        private int[] inverse_map_var(int index, char symbol)
        {
            int local_offset = 0;
            switch (symbol)
            {
                case 'p':
                    {
                        local_offset = 0;
                        break;
                    }
                case 'u':
                    {
                        local_offset = 1;
                        break;
                    }
                case 'v':
                    {
                        local_offset = 2;
                        break;
                    }
            }
            int block = (index - local_offset) / 3;
            int local_j = block % eq_nodes_y;
            int local_i = (block - local_j) / eq_nodes_y;
            int[] output = new int[2];
            output[0] = local_i + 1;
            output[1] = local_j + 1;
            return output;
        }
    }
}
