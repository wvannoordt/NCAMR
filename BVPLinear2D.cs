/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace m_435_NCAMR
{
    class BVPLinear2D
    {
        private bool console_output = false;
        private BoundaryConditions boundary_conditions;
        public BoundaryConditions Boundary_Conditions
        {
            get { return boundary_conditions; }
            set { boundary_conditions = value; }
        }
        LinearOperatorOrder2 lin_operator;
        public LinearOperatorOrder2 LinearOperator
        {
            get { return lin_operator; }
            set { lin_operator = value; }
        }
        NCGrid_Distribution discretization;
        public NCGrid_Distribution Discretization
        {
            get { return discretization; }
            set { discretization = value; }
        }
        public BVPLinear2D(BoundaryConditions _conditions, LinearOperatorOrder2 _operator, NCGrid_Distribution _discretization)
        { 
            lin_operator = _operator;
            discretization = _discretization;
            boundary_conditions = _conditions;
            bool[] all_necessary_compatibility_conditions =
            {
                boundary_conditions.Bounds.Xmax == discretization.Bounds.Xmax,
                boundary_conditions.Bounds.Xmin == discretization.Bounds.Xmin,
                boundary_conditions.Bounds.Ymax == discretization.Bounds.Ymax,
                boundary_conditions.Bounds.Ymin == discretization.Bounds.Ymin,
                boundary_conditions.Xcount == discretization.Xcount,
                boundary_conditions.Ycount == discretization.Ycount
            };
            bool compatible = true;
            foreach (bool i in all_necessary_compatibility_conditions)
            {
                compatible = compatible && i;
            }
            if (!compatible) { throw new Exception("Error: Boundary conditions and initial discretization are incompatible."); }
        }
        private enum BoundaryCase
        {
            PositiveXBoundary,
            NegativeXBoundary,
            PositiveYBoundary,
            NegativeYBoundary,
            LLCorner,
            LRCorner,
            ULCorner,
            URCorner,
            Body
        }
        public void EnableConsoleOutput()
        {
            console_output = true;
        }
        public void DisableConsoleOutput()
        {
            console_output = false;
        }
        public NCGrid_Distribution Solve()
        {
            int m = discretization.Xcount;
            int n = discretization.Ycount;
            int total_nodes = (m - 2) * (n - 2);
            NCGrid_Distribution dist = discretization.Clone();
            dist.ApplyBoundary(boundary_conditions);
            Matrix system_matrix = new Matrix(total_nodes, total_nodes);
            Matrix RHS = new Matrix(total_nodes, 1);
            int neg_y = 0;
            int pos_y = 0;
            int neg_x = 0;
            int pos_x = 0;
            int body = 0;
            int ll_corner = 0;
            int lr_corner = 0;
            int ur_corner = 0;
            int ul_corner = 0;
            if (console_output) { Console.WriteLine("Populating linear system..."); }
            for (int i = 1; i < m - 1; i++)
            {
                if (console_output && i % (m - 1) / 13 == 0)
                {
                    Console.WriteLine((100 * i / (m - 1)).ToString() + "%");
                }
                for (int j = 1; j < n - 1; j++)
                {
                    double rhs_here = 0;
                    //for now, assume zero forcing function.
                    BoundaryCase _case;
                    bool interior = !isCloseToBoundary(i, j, m, n, out _case);
                    int surplusi = 1;
                    int surplusj = 1;
                    Matrix b = dist.GetTaylorSystemCoeffs(i, j, surplusi, surplusj);
                    double uijterm = 0;
                    double ui1jterm = 0;
                    double uij1term = 0;
                    double ui_1jterm = 0;
                    double uij_1term = 0;
                    double usurplusterm = 0;
                    int row = (m - 2) * (i - 1) + j - 1;
                    for (int h = 0; h < 5; h++)
                    {
                        ui1jterm += b[h, 0] * lin_operator[h];
                        uij1term += b[h, 1] * lin_operator[h];
                        ui_1jterm += b[h, 2] * lin_operator[h];
                        uij_1term += b[h, 3] * lin_operator[h];
                        usurplusterm += b[h, 4] * lin_operator[h];
                        double temp_ij = 0;
                        for (int k = 0; k < 5; k++)
                        {
                            temp_ij += b[h, k];
                        }
                        uijterm += lin_operator[h] * temp_ij;
                    }
                    switch (_case)
                    {
                        case BoundaryCase.NegativeYBoundary:
                            {
                                rhs_here -= dist[i, j - 1].Value * uij_1term;
                                uij_1term = 0;
                                neg_y++;
                                break;
                            }
                        case BoundaryCase.PositiveYBoundary:
                            {
                                rhs_here -= dist[i, j + 1].Value * uij1term;
                                rhs_here -= dist[i + surplusi, j + surplusj].Value * usurplusterm;
                                usurplusterm = 0;
                                uij1term = 0;
                                pos_y++;
                                break;
                            }
                        case BoundaryCase.NegativeXBoundary:
                            {
                                rhs_here -= dist[i - 1, j].Value * ui_1jterm;
                                ui_1jterm = 0;
                                neg_x++;
                                break;
                            }
                        case BoundaryCase.PositiveXBoundary:
                            {
                                rhs_here -= dist[i + 1, j].Value * ui1jterm;
                                rhs_here -= dist[i + surplusi, j + surplusj].Value * usurplusterm;
                                usurplusterm = 0;
                                ui1jterm = 0;
                                pos_x++;
                                break;
                            }
                        case BoundaryCase.ULCorner:
                            {
                                rhs_here -= dist[i, j + 1].Value * uij1term;
                                rhs_here -= dist[i - 1, j].Value * ui_1jterm;
                                rhs_here -= dist[i + surplusi, j + surplusj].Value * usurplusterm;
                                uij1term = 0;
                                ui_1jterm = 0;
                                usurplusterm = 0;
                                ul_corner++;
                                break;
                            }
                        case BoundaryCase.LLCorner:
                            {
                                rhs_here -= dist[i - 1, j].Value * ui_1jterm;
                                rhs_here -= dist[i, j - 1].Value * uij_1term;
                                ui_1jterm = 0;
                                uij_1term = 0;
                                ll_corner++;
                                break;
                            }
                        case BoundaryCase.URCorner:
                            {
                                rhs_here -= dist[i + 1, j].Value * ui1jterm;
                                rhs_here -= dist[i, j + 1].Value * uij1term;
                                rhs_here -= dist[i + surplusi, j + surplusj].Value * usurplusterm;
                                ui1jterm = 0;
                                uij1term = 0;
                                usurplusterm = 0;
                                ur_corner++;
                                break;
                            }
                        case BoundaryCase.LRCorner:
                            {
                                rhs_here -= dist[i, j - 1].Value * uij_1term;
                                rhs_here -= dist[i + 1, j].Value * ui1jterm;
                                rhs_here -= dist[i + surplusj, j + surplusj].Value * usurplusterm;
                                uij_1term = 0;
                                ui1jterm = 0;
                                usurplusterm = 0;
                                lr_corner++;
                                break;
                            }
                        case BoundaryCase.Body:
                            {
                                body++;
                                break;
                            }
                    }
                    system_matrix[row, map2Dto1D(i - 1, j - 1, m - 1, n - 1)] = -1 * uijterm;
                    if (ui1jterm != 0)
                    {
                        system_matrix[row, map2Dto1D(i, j - 1, m - 1, n - 1)] = ui1jterm;
                    }
                    if (uij1term != 0)
                    {
                        system_matrix[row, map2Dto1D(i - 1, j, m - 1, n - 1)] = uij1term;
                    }
                    if (ui_1jterm != 0)
                    {
                        system_matrix[row, map2Dto1D(i - 2, j - 1, m - 1, n - 1)] = ui_1jterm;
                    }
                    if (uij_1term != 0)
                    {
                        system_matrix[row, map2Dto1D(i - 1, j - 2, m - 1, n - 1)] = uij_1term;
                    }
                    if (usurplusterm != 0)
                    {
                        system_matrix[row, map2Dto1D(i - 1 + surplusi, j - 1 + surplusj, m - 1, n - 1)] = usurplusterm;
                    }
                    RHS[row] = rhs_here;
                }
            }
            if (console_output)
            {
                Console.WriteLine("System populated.");
            }
            Matrix results = RunExteriorSolve(system_matrix, RHS);
            for (int i = 0; i < total_nodes; i++)
            {
                int offset = n - 2;
                int J = i % offset;
                int I = (i - J) / offset;
                dist.assign_value_at(I + 1, J + 1, results[i]);
            }
            return dist;
        }
        public NCGrid_Distribution SolveKaczMarzExt(int maxiterations, double tolerance, int radius)
        {
            int m = discretization.Xcount;
            int n = discretization.Ycount;
            int total_nodes = (m - 2) * (n - 2);
            NCGrid_Distribution dist = discretization.Clone();
            dist.ApplyBoundary(boundary_conditions);
            Matrix system_matrix = new Matrix(total_nodes, total_nodes);
            Matrix RHS = new Matrix(total_nodes, 1);
            int neg_y = 0;
            int pos_y = 0;
            int neg_x = 0;
            int pos_x = 0;
            int body = 0;
            int ll_corner = 0;
            int lr_corner = 0;
            int ur_corner = 0;
            int ul_corner = 0;
            if (console_output) { Console.WriteLine("Populating linear system..."); }
            for (int i = 1; i < m - 1; i++)
            {
                if (console_output && i % (m - 1) / 13 == 0)
                {
                    Console.WriteLine((100 * i / (m - 1)).ToString() + "%");
                }
                for (int j = 1; j < n - 1; j++)
                {
                    double rhs_here = 0;
                    //for now, assume zero forcing function.
                    BoundaryCase _case;
                    bool interior = !isCloseToBoundary(i, j, m, n, out _case);
                    int surplusi = 1;
                    int surplusj = 1;
                    Matrix b = dist.GetTaylorSystemCoeffs(i, j, surplusi, surplusj);
                    double uijterm = 0;
                    double ui1jterm = 0;
                    double uij1term = 0;
                    double ui_1jterm = 0;
                    double uij_1term = 0;
                    double usurplusterm = 0;
                    int row = (m - 2) * (i - 1) + j - 1;
                    for (int h = 0; h < 5; h++)
                    {
                        ui1jterm += b[h, 0] * lin_operator[h];
                        uij1term += b[h, 1] * lin_operator[h];
                        ui_1jterm += b[h, 2] * lin_operator[h];
                        uij_1term += b[h, 3] * lin_operator[h];
                        usurplusterm += b[h, 4] * lin_operator[h];
                        double temp_ij = 0;
                        for (int k = 0; k < 5; k++)
                        {
                            temp_ij += b[h, k];
                        }
                        uijterm += lin_operator[h] * temp_ij;
                    }
                    switch (_case)
                    {
                        case BoundaryCase.NegativeYBoundary:
                            {
                                rhs_here -= dist[i, j - 1].Value * uij_1term;
                                uij_1term = 0;
                                neg_y++;
                                break;
                            }
                        case BoundaryCase.PositiveYBoundary:
                            {
                                rhs_here -= dist[i, j + 1].Value * uij1term;
                                rhs_here -= dist[i + surplusi, j + surplusj].Value * usurplusterm;
                                usurplusterm = 0;
                                uij1term = 0;
                                pos_y++;
                                break;
                            }
                        case BoundaryCase.NegativeXBoundary:
                            {
                                rhs_here -= dist[i - 1, j].Value * ui_1jterm;
                                ui_1jterm = 0;
                                neg_x++;
                                break;
                            }
                        case BoundaryCase.PositiveXBoundary:
                            {
                                rhs_here -= dist[i + 1, j].Value * ui1jterm;
                                rhs_here -= dist[i + surplusi, j + surplusj].Value * usurplusterm;
                                usurplusterm = 0;
                                ui1jterm = 0;
                                pos_x++;
                                break;
                            }
                        case BoundaryCase.ULCorner:
                            {
                                rhs_here -= dist[i, j + 1].Value * uij1term;
                                rhs_here -= dist[i - 1, j].Value * ui_1jterm;
                                rhs_here -= dist[i + surplusi, j + surplusj].Value * usurplusterm;
                                uij1term = 0;
                                ui_1jterm = 0;
                                usurplusterm = 0;
                                ul_corner++;
                                break;
                            }
                        case BoundaryCase.LLCorner:
                            {
                                rhs_here -= dist[i - 1, j].Value * ui_1jterm;
                                rhs_here -= dist[i, j - 1].Value * uij_1term;
                                ui_1jterm = 0;
                                uij_1term = 0;
                                ll_corner++;
                                break;
                            }
                        case BoundaryCase.URCorner:
                            {
                                rhs_here -= dist[i + 1, j].Value * ui1jterm;
                                rhs_here -= dist[i, j + 1].Value * uij1term;
                                rhs_here -= dist[i + surplusi, j + surplusj].Value * usurplusterm;
                                ui1jterm = 0;
                                uij1term = 0;
                                usurplusterm = 0;
                                ur_corner++;
                                break;
                            }
                        case BoundaryCase.LRCorner:
                            {
                                rhs_here -= dist[i, j - 1].Value * uij_1term;
                                rhs_here -= dist[i + 1, j].Value * ui1jterm;
                                rhs_here -= dist[i + surplusj, j + surplusj].Value * usurplusterm;
                                uij_1term = 0;
                                ui1jterm = 0;
                                usurplusterm = 0;
                                lr_corner++;
                                break;
                            }
                        case BoundaryCase.Body:
                            {
                                body++;
                                break;
                            }
                    }
                    system_matrix[row, map2Dto1D(i - 1, j - 1, m - 1, n - 1)] = -1 * uijterm;
                    if (ui1jterm != 0)
                    {
                        system_matrix[row, map2Dto1D(i, j - 1, m - 1, n - 1)] = ui1jterm;
                    }
                    if (uij1term != 0)
                    {
                        system_matrix[row, map2Dto1D(i - 1, j, m - 1, n - 1)] = uij1term;
                    }
                    if (ui_1jterm != 0)
                    {
                        system_matrix[row, map2Dto1D(i - 2, j - 1, m - 1, n - 1)] = ui_1jterm;
                    }
                    if (uij_1term != 0)
                    {
                        system_matrix[row, map2Dto1D(i - 1, j - 2, m - 1, n - 1)] = uij_1term;
                    }
                    if (usurplusterm != 0)
                    {
                        system_matrix[row, map2Dto1D(i - 1 + surplusi, j - 1 + surplusj, m - 1, n - 1)] = usurplusterm;
                    }
                    RHS[row] = rhs_here;
                }
            }
            if (console_output)
            {
                Console.WriteLine("System populated.");
            }
            Matrix results = Matrix.ext_kaczmarz_radius(system_matrix, RHS, tolerance, maxiterations, radius);
            for (int i = 0; i < total_nodes; i++)
            {
                int offset = n - 2;
                int J = i % offset;
                int I = (i - J) / offset;
                dist.assign_value_at(I + 1, J + 1, results[i]);
            }
            return dist;
        }
        public NCGrid_Distribution Solve(Matrix.SystemSolvingScheme scheme)
        {
            int m = discretization.Xcount;
            int n = discretization.Ycount;
            int total_nodes = (m - 2) * (n - 2);
            NCGrid_Distribution dist = discretization.Clone();
            dist.ApplyBoundary(boundary_conditions);
            Matrix system_matrix = new Matrix(total_nodes, total_nodes);
            Matrix RHS = new Matrix(total_nodes, 1);
            int neg_y = 0;
            int pos_y = 0;
            int neg_x = 0;
            int pos_x = 0;
            int body = 0;
            int ll_corner = 0;
            int lr_corner = 0;
            int ur_corner = 0;
            int ul_corner = 0;
            if (console_output) { Console.WriteLine("Populating linear system..."); }
            for (int i = 1; i < m - 1; i++)
            {
                if (console_output && i % (m - 1) / 13 == 0)
                {
                    Console.WriteLine((100 * i / (m - 1)).ToString() + "%");
                }
                for (int j = 1; j < n - 1; j++)
                {
                    double rhs_here = 0;
                    //for now, assume zero forcing function.
                    BoundaryCase _case;
                    bool interior = !isCloseToBoundary(i, j, m, n, out _case);
                    int surplusi = 1;
                    int surplusj = 1;
                    Matrix b = dist.GetTaylorSystemCoeffs(i, j, surplusi, surplusj);
                    double uijterm = 0;
                    double ui1jterm = 0;
                    double uij1term = 0;
                    double ui_1jterm = 0;
                    double uij_1term = 0;
                    double usurplusterm = 0;
                    int row = (m - 2) * (i - 1) + j - 1;
                    for (int h = 0; h < 5; h++)
                    {
                        ui1jterm += b[h, 0] * lin_operator[h];
                        uij1term += b[h, 1] * lin_operator[h];
                        ui_1jterm += b[h, 2] * lin_operator[h];
                        uij_1term += b[h, 3] * lin_operator[h];
                        usurplusterm += b[h, 4] * lin_operator[h];
                        double temp_ij = 0;
                        for (int k = 0; k < 5; k++)
                        {
                            temp_ij += b[h, k];
                        }
                        uijterm += lin_operator[h] * temp_ij;
                    }
                    switch (_case)
                    {
                        case BoundaryCase.NegativeYBoundary:
                            {
                                rhs_here -= dist[i, j - 1].Value * uij_1term;
                                uij_1term = 0;
                                neg_y++;
                                break;
                            }
                        case BoundaryCase.PositiveYBoundary:
                            {
                                rhs_here -= dist[i, j + 1].Value * uij1term;
                                rhs_here -= dist[i + surplusi, j + surplusj].Value * usurplusterm;
                                usurplusterm = 0;
                                uij1term = 0;
                                pos_y++;
                                break;
                            }
                        case BoundaryCase.NegativeXBoundary:
                            {
                                rhs_here -= dist[i - 1, j].Value * ui_1jterm;
                                ui_1jterm = 0;
                                neg_x++;
                                break;
                            }
                        case BoundaryCase.PositiveXBoundary:
                            {
                                rhs_here -= dist[i + 1, j].Value * ui1jterm;
                                rhs_here -= dist[i + surplusi, j + surplusj].Value * usurplusterm;
                                usurplusterm = 0;
                                ui1jterm = 0;
                                pos_x++;
                                break;
                            }
                        case BoundaryCase.ULCorner:
                            {
                                rhs_here -= dist[i, j + 1].Value * uij1term;
                                rhs_here -= dist[i - 1, j].Value * ui_1jterm;
                                rhs_here -= dist[i + surplusi, j + surplusj].Value * usurplusterm;
                                uij1term = 0;
                                ui_1jterm = 0;
                                usurplusterm = 0;
                                ul_corner++;
                                break;
                            }
                        case BoundaryCase.LLCorner:
                            {
                                rhs_here -= dist[i - 1, j].Value * ui_1jterm;
                                rhs_here -= dist[i, j - 1].Value * uij_1term;
                                ui_1jterm = 0;
                                uij_1term = 0;
                                ll_corner++;
                                break;
                            }
                        case BoundaryCase.URCorner:
                            {
                                rhs_here -= dist[i + 1, j].Value * ui1jterm;
                                rhs_here -= dist[i, j + 1].Value * uij1term;
                                rhs_here -= dist[i + surplusi, j + surplusj].Value * usurplusterm;
                                ui1jterm = 0;
                                uij1term = 0;
                                usurplusterm = 0;
                                ur_corner++;
                                break;
                            }
                        case BoundaryCase.LRCorner:
                            {
                                rhs_here -= dist[i, j - 1].Value * uij_1term;
                                rhs_here -= dist[i + 1, j].Value * ui1jterm;
                                rhs_here -= dist[i + surplusj, j + surplusj].Value * usurplusterm;
                                uij_1term = 0;
                                ui1jterm = 0;
                                usurplusterm = 0;
                                lr_corner++;
                                break;
                            }
                        case BoundaryCase.Body:
                            {
                                body++;
                                break;
                            }
                    }
                    system_matrix[row, map2Dto1D(i - 1, j - 1, m - 1, n - 1)] = -1 * uijterm;
                    if (ui1jterm != 0)
                    {
                        system_matrix[row, map2Dto1D(i, j - 1, m - 1, n - 1)] = ui1jterm;
                    }
                    if (uij1term != 0)
                    {
                        system_matrix[row, map2Dto1D(i - 1, j, m - 1, n - 1)] = uij1term;
                    }
                    if (ui_1jterm != 0)
                    {
                        system_matrix[row, map2Dto1D(i - 2, j - 1, m - 1, n - 1)] = ui_1jterm;
                    }
                    if (uij_1term != 0)
                    {
                        system_matrix[row, map2Dto1D(i - 1, j - 2, m - 1, n - 1)] = uij_1term;
                    }
                    if (usurplusterm != 0)
                    {
                        system_matrix[row, map2Dto1D(i - 1 + surplusi, j - 1 + surplusj, m - 1, n - 1)] = usurplusterm;
                    }
                    RHS[row] = rhs_here;
                }
            }
            if (console_output)
            {
                Console.WriteLine("System populated.");
            }
            Matrix results = Matrix.solve_system(system_matrix, RHS, scheme);
            for (int i = 0; i < total_nodes; i++)
            {
                int offset = n - 2;
                int J = i % offset;
                int I = (i - J) / offset;
                dist.assign_value_at(I + 1, J + 1, results[i]);
            }
            return dist;
        }
        private int maxval(List<int> f)
        {
            int max = f[0];
            for (int i = 1; i < f.Count; i++)
            {
                if (f[i] > max) { max = f[i]; }
            }
            return max;
        }
        private void buildsystem(out Matrix system_matrix, out Matrix RHS)
        {
            int m = discretization.Xcount;
            int n = discretization.Ycount;
            int total_nodes = (m - 2) * (n - 2);
            NCGrid_Distribution dist = discretization.Clone();
            dist.ApplyBoundary(boundary_conditions);
            system_matrix = new Matrix(total_nodes, total_nodes);
            RHS = new Matrix(total_nodes, 1);
            int neg_y = 0;
            int pos_y = 0;
            int neg_x = 0;
            int pos_x = 0;
            int body = 0;
            int ll_corner = 0;
            int lr_corner = 0;
            int ur_corner = 0;
            int ul_corner = 0;
            if (console_output) { Console.WriteLine("Populating linear system..."); }
            for (int i = 1; i < m - 1; i++)
            {
                if (console_output && i % (m - 1) / 13 == 0)
                {
                    Console.WriteLine((100 * i / (m - 1)).ToString() + "%");
                }
                for (int j = 1; j < n - 1; j++)
                {
                    double rhs_here = 0;
                    //for now, assume zero forcing function.
                    BoundaryCase _case;
                    bool interior = !isCloseToBoundary(i, j, m, n, out _case);
                    int surplusi = 1;
                    int surplusj = 1;
                    Matrix b = dist.GetTaylorSystemCoeffs(i, j, surplusi, surplusj);
                    double uijterm = 0;
                    double ui1jterm = 0;
                    double uij1term = 0;
                    double ui_1jterm = 0;
                    double uij_1term = 0;
                    double usurplusterm = 0;
                    int row = (m - 2) * (i - 1) + j - 1;
                    for (int h = 0; h < 5; h++)
                    {
                        ui1jterm += b[h, 0] * lin_operator[h];
                        uij1term += b[h, 1] * lin_operator[h];
                        ui_1jterm += b[h, 2] * lin_operator[h];
                        uij_1term += b[h, 3] * lin_operator[h];
                        usurplusterm += b[h, 4] * lin_operator[h];
                        double temp_ij = 0;
                        for (int k = 0; k < 5; k++)
                        {
                            temp_ij += b[h, k];
                        }
                        uijterm += lin_operator[h] * temp_ij;
                    }
                    switch (_case)
                    {
                        case BoundaryCase.NegativeYBoundary:
                            {
                                rhs_here -= dist[i, j - 1].Value * uij_1term;
                                uij_1term = 0;
                                neg_y++;
                                break;
                            }
                        case BoundaryCase.PositiveYBoundary:
                            {
                                rhs_here -= dist[i, j + 1].Value * uij1term;
                                rhs_here -= dist[i + surplusi, j + surplusj].Value * usurplusterm;
                                usurplusterm = 0;
                                uij1term = 0;
                                pos_y++;
                                break;
                            }
                        case BoundaryCase.NegativeXBoundary:
                            {
                                rhs_here -= dist[i - 1, j].Value * ui_1jterm;
                                ui_1jterm = 0;
                                neg_x++;
                                break;
                            }
                        case BoundaryCase.PositiveXBoundary:
                            {
                                rhs_here -= dist[i + 1, j].Value * ui1jterm;
                                rhs_here -= dist[i + surplusi, j + surplusj].Value * usurplusterm;
                                usurplusterm = 0;
                                ui1jterm = 0;
                                pos_x++;
                                break;
                            }
                        case BoundaryCase.ULCorner:
                            {
                                rhs_here -= dist[i, j + 1].Value * uij1term;
                                rhs_here -= dist[i - 1, j].Value * ui_1jterm;
                                rhs_here -= dist[i + surplusi, j + surplusj].Value * usurplusterm;
                                uij1term = 0;
                                ui_1jterm = 0;
                                usurplusterm = 0;
                                ul_corner++;
                                break;
                            }
                        case BoundaryCase.LLCorner:
                            {
                                rhs_here -= dist[i - 1, j].Value * ui_1jterm;
                                rhs_here -= dist[i, j - 1].Value * uij_1term;
                                ui_1jterm = 0;
                                uij_1term = 0;
                                ll_corner++;
                                break;
                            }
                        case BoundaryCase.URCorner:
                            {
                                rhs_here -= dist[i + 1, j].Value * ui1jterm;
                                rhs_here -= dist[i, j + 1].Value * uij1term;
                                rhs_here -= dist[i + surplusi, j + surplusj].Value * usurplusterm;
                                ui1jterm = 0;
                                uij1term = 0;
                                usurplusterm = 0;
                                ur_corner++;
                                break;
                            }
                        case BoundaryCase.LRCorner:
                            {
                                rhs_here -= dist[i, j - 1].Value * uij_1term;
                                rhs_here -= dist[i + 1, j].Value * ui1jterm;
                                rhs_here -= dist[i + surplusj, j + surplusj].Value * usurplusterm;
                                uij_1term = 0;
                                ui1jterm = 0;
                                usurplusterm = 0;
                                lr_corner++;
                                break;
                            }
                        case BoundaryCase.Body:
                            {
                                body++;
                                break;
                            }
                    }
                    system_matrix[row, map2Dto1D(i - 1, j - 1, m - 1, n - 1)] = -1 * uijterm;
                    if (ui1jterm != 0)
                    {
                        system_matrix[row, map2Dto1D(i, j - 1, m - 1, n - 1)] = ui1jterm;
                    }
                    if (uij1term != 0)
                    {
                        system_matrix[row, map2Dto1D(i - 1, j, m - 1, n - 1)] = uij1term;
                    }
                    if (ui_1jterm != 0)
                    {
                        system_matrix[row, map2Dto1D(i - 2, j - 1, m - 1, n - 1)] = ui_1jterm;
                    }
                    if (uij_1term != 0)
                    {
                        system_matrix[row, map2Dto1D(i - 1, j - 2, m - 1, n - 1)] = uij_1term;
                    }
                    if (usurplusterm != 0)
                    {
                        system_matrix[row, map2Dto1D(i - 1 + surplusi, j - 1 + surplusj, m - 1, n - 1)] = usurplusterm;
                    }
                    RHS[row] = rhs_here;
                }
            }
            if (console_output)
            {
                Console.WriteLine("System populated.");
            }
        }
        public NCGrid_Distribution SolveSRDD()
        {
            Matrix system_matrix, RHS;
            NCGrid_Distribution dist = discretization.Clone();
            dist.ApplyBoundary(boundary_conditions);
            buildsystem(out system_matrix, out RHS);
            Matrix results = RunExteriorSolve(system_matrix, RHS);
            for (int i = 0; i < RHS.Rows; i++)
            {
                int offset = discretization.Ycount - 2;
                int J = i % offset;
                int I = (i - J) / offset;
                dist.assign_value_at(I + 1, J + 1, results[i]);
            }
            return dist;
        }
        private Matrix RunExteriorSolve(Matrix system_matrix, Matrix RHS)
        {
            system_matrix.exportToFile(@"C:\Users\Will\Desktop\Folders\MATH435\repo\matrices\sys.csv");
            Matrix system = Matrix.ConcatenateRight(system_matrix, RHS);
            Matrix result = new Matrix(system.Rows, 1);
            int m = system.Rows;
            for (int diagrow = 0; diagrow < m; diagrow++)
            {
                if (console_output && diagrow % (m/17) == 0)
                {
                    Console.WriteLine(diagrow.ToString() + " of " + m.ToString() + " rows normalized (" + ((100*diagrow)/m).ToString() + "%).");
                }
                double scl1 = 1 / system[diagrow, diagrow];
                system.RowScale(diagrow, scl1);
                for (int i = 0; i < m; i++)
                {
                    if (i != diagrow && system[i, diagrow] != 0)
                    {
                        double scl2 = -1 * system[i, diagrow];
                        system.RowAdd(i, diagrow, scl2);
                    }
                }
                //foreach (int i in targets[diagrow])
                //{
                //    double scl2 = -1 * system[i, diagrow];
                //    system.RowAdd(i, diagrow, scl2);
                //}
            }
            for (int i = 0; i < m; i++)
            {
                result[i] = system[i, m];
            }
            return result;
        }
        public NCGrid_Distribution solve_iterative_morph(int max_morph_count, double size)
        {
            NCGrid_Distribution soln = Solve();
            for (int i = 0; i < max_morph_count - 1; i++)
            {
                soln.QuickSketch("soln-" + i.ToString());
                soln.ApplyMeshMorphGA(size);
                discretization = soln.Clone();
                soln = Solve();
            }
            return soln;

        }
        public NCGrid_Distribution SolveSlow()
        {
            int m = discretization.Xcount;
            int n = discretization.Ycount;
            int total_nodes = (m - 2) * (n - 2);
            //NCGrid_Distribution dist = new NCGrid_Distribution(boundary_conditions.Bounds, discretization.Xcount, discretization.Ycount);
            NCGrid_Distribution dist = discretization.Clone();
            dist.ApplyBoundary(boundary_conditions);
            Matrix system_matrix = new Matrix(total_nodes,total_nodes);
            Matrix RHS = new Matrix(total_nodes, 1);
            int neg_y = 0;
            int pos_y = 0;
            int neg_x = 0;
            int pos_x = 0;
            int body = 0;
            int ll_corner = 0;
            int lr_corner = 0;
            int ur_corner = 0;
            int ul_corner = 0;
            List<string> indices = new List<string>();
            if (console_output) { Console.WriteLine("Populating linear system..."); }
            for (int i = 1; i < m-1; i++)
            {
                if (console_output && i % (m-1)/13 == 0)
                {
                    Console.WriteLine((100 * i / (m - 1)).ToString() + "%");
                }
                for (int j = 1; j < n-1; j++)
                {
                    double rhs_here = 0;
                    indices.Add(i.ToString() + "," + j.ToString());
                    //for now, assume zero forcing function.
                    BoundaryCase _case;
                    bool interior = !isCloseToBoundary(i, j, m, n, out _case);
                    int surplusi = 1;
                    int surplusj = 1;
                    Matrix b = dist.GetTaylorSystemCoeffs(i, j, surplusi, surplusj);
                    double uijterm = 0;
                    double ui1jterm = 0;
                    double uij1term = 0;
                    double ui_1jterm = 0;
                    double uij_1term = 0;
                    double usurplusterm = 0;
                    int row = (m - 2) * (i - 1) + j - 1;
                    for (int h = 0; h < 5; h++)
                    {
                        ui1jterm += b[h, 0] * lin_operator[h];
                        uij1term += b[h, 1] * lin_operator[h];
                        ui_1jterm += b[h, 2] * lin_operator[h];
                        uij_1term += b[h, 3] * lin_operator[h];
                        usurplusterm += b[h, 4] * lin_operator[h];
                        double temp_ij = 0;
                        for (int k = 0; k < 5; k++)
                        {
                            temp_ij += b[h, k];
                        }
                        uijterm += lin_operator[h] * temp_ij;
                    }
                    switch (_case)
                    {
                        case BoundaryCase.NegativeYBoundary:
                            {
                                rhs_here -= dist[i, j - 1].Value * uij_1term;
                                uij_1term = 0;
                                neg_y++;
                                break;
                            }
                        case BoundaryCase.PositiveYBoundary:
                            {
                                rhs_here -= dist[i, j + 1].Value * uij1term;
                                rhs_here -= dist[i + surplusi, j + surplusj].Value * usurplusterm;
                                usurplusterm = 0;
                                uij1term = 0;
                                pos_y++;
                                break;
                            }
                        case BoundaryCase.NegativeXBoundary:
                            {
                                rhs_here -= dist[i - 1, j].Value * ui_1jterm;
                                ui_1jterm = 0;
                                neg_x++;
                                break;
                            }
                        case BoundaryCase.PositiveXBoundary:
                            {
                                rhs_here -= dist[i + 1, j].Value * ui1jterm;
                                rhs_here -= dist[i + surplusi, j + surplusj].Value * usurplusterm;
                                usurplusterm = 0;
                                ui1jterm = 0;
                                pos_x++;
                                break;
                            }
                        case BoundaryCase.ULCorner:
                            {
                                rhs_here -= dist[i, j + 1].Value * uij1term;
                                rhs_here -= dist[i - 1, j].Value * ui_1jterm;
                                rhs_here -= dist[i + surplusi, j + surplusj].Value * usurplusterm;
                                uij1term = 0;
                                ui_1jterm = 0;
                                usurplusterm = 0;
                                ul_corner++;
                                break;
                            }
                        case BoundaryCase.LLCorner:
                            {
                                rhs_here -= dist[i - 1, j].Value * ui_1jterm;
                                rhs_here -= dist[i, j - 1].Value * uij_1term;
                                ui_1jterm = 0;
                                uij_1term = 0;
                                ll_corner++;
                                break;
                            }
                        case BoundaryCase.URCorner:
                            {
                                rhs_here -= dist[i + 1, j].Value * ui1jterm;
                                rhs_here -= dist[i, j + 1].Value * uij1term;
                                rhs_here -= dist[i + surplusi, j + surplusj].Value * usurplusterm;
                                ui1jterm = 0;
                                uij1term = 0;
                                usurplusterm = 0;
                                ur_corner++;
                                break;
                            }
                        case BoundaryCase.LRCorner:
                            {
                                rhs_here -= dist[i, j - 1].Value * uij_1term;
                                rhs_here -= dist[i + 1, j].Value * ui1jterm;
                                rhs_here -= dist[i + surplusj, j + surplusj].Value * usurplusterm;
                                uij_1term = 0;
                                ui1jterm = 0;
                                usurplusterm = 0;
                                lr_corner++;
                                break;
                            }
                        case BoundaryCase.Body:
                            {
                                body++;
                                break;
                            }
                    }
                    system_matrix[row, map2Dto1D(i - 1, j - 1, m - 1, n - 1)] = -1*uijterm;
                    if (ui1jterm != 0)
                    {
                        system_matrix[row, map2Dto1D(i, j - 1, m - 1, n - 1)] = ui1jterm;
                    }
                    if (uij1term != 0)
                    {
                        system_matrix[row, map2Dto1D(i-1, j, m - 1, n - 1)] = uij1term;
                    }
                    if (ui_1jterm != 0)
                    {
                        system_matrix[row, map2Dto1D(i-2, j - 1, m - 1, n - 1)] = ui_1jterm;
                    }
                    if (uij_1term != 0)
                    {
                        system_matrix[row, map2Dto1D(i - 1, j - 2, m - 1, n - 1)] = uij_1term;
                    }
                    if (usurplusterm != 0)
                    {
                        system_matrix[row, map2Dto1D(i - 1 + surplusi, j - 1 + surplusj, m - 1, n - 1)] = usurplusterm;
                    }
                    RHS[row] = rhs_here;
                }
            }
            if (console_output)
            {
                Console.WriteLine("System populated.");
            }
            Matrix results = Matrix.solve_system(system_matrix, RHS, Matrix.SystemSolvingScheme.Basic, true, true);
            for (int i = 0; i < total_nodes; i++)
            {
                int offset = n - 2;
                int J = i % offset;
                int I = (i-J)/offset;
                dist.assign_value_at(I+1,J+1, results[i]);
            }
            return dist;
        }
        private bool isCloseToBoundary(int i, int j, int xcount, int ycount, out BoundaryCase _case)
        {
            bool onTop = j == ycount - 2;
            bool onBottom = j == 1;
            bool onLeft = i == 1;
            bool onRight = i == xcount - 2;
            _case = BoundaryCase.Body;
            if (onBottom && !onRight && !onLeft) { _case = BoundaryCase.NegativeYBoundary; }
            if (onBottom && onRight) { _case = BoundaryCase.LRCorner; }
            if (onBottom && onLeft) { _case = BoundaryCase.LLCorner; }
            if (!onBottom&&!onTop&&onLeft) { _case = BoundaryCase.NegativeXBoundary; }
            if (onTop && !onRight && !onLeft) { _case = BoundaryCase.PositiveYBoundary; }
            if (onTop && onRight) { _case = BoundaryCase.URCorner; }
            if (onTop && onLeft) { _case = BoundaryCase.ULCorner; }
            if (!onBottom && !onTop && onRight) { _case = BoundaryCase.PositiveXBoundary; }
            return onTop || onBottom || onLeft || onRight;
        }
        private int map2Dto1D(int i, int j, int icount, int jcount)
        {
            return i * (jcount-1) + j;
        }
        private int[] map1Dto2D(int q, int icount, int jcount)
        {
            int j = q & jcount;
            int i = (q - j) / jcount;
            int[] output = { i, j };
            return output;
        }
    }
}
