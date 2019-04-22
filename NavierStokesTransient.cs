/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Threading.Tasks;

namespace m_435_NCAMR
{
    class NavierStokesTransient
    {
        int time_steps;
        NCGrid_Distribution current_ustar, current_vstar, initial_u, initial_v;
        BoundaryConditions p_boundary, u_boundary, v_boundary;
        bool enable_console_output = true;
        double deltat;
        FluidProperties properties;
        public bool ConsoleOutputEnabled
        {
            get { return enable_console_output; }
            set { enable_console_output = value; }
        }
        public NavierStokesTransient(int timesteps, NCGrid_Distribution u_initial, NCGrid_Distribution v_initial, BoundaryConditions p, BoundaryConditions u, BoundaryConditions v, FluidProperties props, double dt)
        {
            deltat = dt;
            properties = props;
            current_ustar = new NCGrid_Distribution(u_initial.Bounds, u_initial.Xcount, u_initial.Ycount);
            current_vstar = new NCGrid_Distribution(v_initial.Bounds, v_initial.Xcount, v_initial.Ycount);
            current_ustar.SetConstant(1);
            current_vstar.SetConstant(1);
            time_steps = timesteps;
            p_boundary = p;
            u_boundary = u;
            v_boundary = v;
            initial_u = u_initial;
            initial_v = v_initial;
        }
        public void Solve()
        {
            Solve("solution");
        }
        public void Solve(string title) // return p,u,v
        {
            string path = Paths.NavSolutionsRepo + "//" + title;
            if (Directory.Exists(path))
            {
                int n = 0;
                while (Directory.Exists(path + n.ToString()))
                {
                    n++;
                }
                path += n.ToString();
            }
            NCGrid_Distribution[] last;
            Directory.CreateDirectory(path);
            Directory.CreateDirectory(path + "\\p");
            Directory.CreateDirectory(path + "\\u");
            Directory.CreateDirectory(path + "\\v");
            double t = 0;
            NavierStokesProblemSingleFrame frame = new NavierStokesProblemSingleFrame(current_ustar, current_vstar, initial_u, initial_v, p_boundary, u_boundary, v_boundary, properties, deltat);
            NCGrid_Distribution[] sols = frame.Solve();
            sols[0].WriteToFile(path + "\\p\\" + title + "-pressure-" + bufferint(0) + ".dist", false, true, true);
            sols[1].WriteToFile(path + "\\u\\" + title + "-x-vel-" + bufferint(0) + ".dist", false, true, true);
            sols[2].WriteToFile(path + "\\v\\" + title + "-y-vel-" + bufferint(0) + ".dist", false, true, true);
            last = sols;
            for (int i = 0; i < time_steps; i++)
            {
                t += deltat;
                NavierStokesProblemSingleFrame newframe = new NavierStokesProblemSingleFrame(last[1], last[2], last[1], last[2], u_boundary, v_boundary, p_boundary, properties, deltat);
                NCGrid_Distribution[] solsnew = newframe.Solve();
                solsnew[0].WriteToFile(path + "\\p\\" + title + "-pressure-" + bufferint(i+1) + ".dist", false, true, true);
                solsnew[1].WriteToFile(path + "\\u\\" + title + "-x-vel-" + bufferint(i+1) + ".dist", false, true, true);
                solsnew[2].WriteToFile(path + "\\v\\" + title + "-y-vel-" + bufferint(i+1) + ".dist", false, true, true);
                if (enable_console_output)
                {
                    Console.ForegroundColor = ConsoleColor.Green;
                    Console.WriteLine("Step " + i.ToString() + " of " + time_steps.ToString() + " complete." );
                    Console.ForegroundColor = ConsoleColor.Gray;
                }
                last[0] = solsnew[0].Clone();
                last[1] = solsnew[1].Clone();
                last[2] = solsnew[2].Clone();
            }
        }
        private string bufferint(int t)
        {
            int numzeros = time_steps.ToString().Length - t.ToString().Length;
            string output = string.Empty;
            for (int i = 0; i < numzeros; i++)
            {
                output += "0";
            }
            output += t.ToString();
            return output;
        }
    }
}
