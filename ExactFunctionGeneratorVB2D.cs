/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using System.Threading;
using System.IO;

namespace m_435_NCAMR
{
    static class ExactFunctionGeneratorVB2D
    {
        private static string CSCRIPT = @"C:\Windows\SysWOW64\cscript.exe";
        private static string script_template = @"C:\Users\Will\Desktop\Folders\MATH435\repo\visual-basic-templates\function-generator-template.vbs";
        private static string script_disc_template = @"C:\Users\Will\Desktop\Folders\MATH435\repo\visual-basic-templates\function-generator-discretization-template.vbs";
        public static void GenerateFunction(string vb_syntax_eval, RBounds2D bounds, int xcount, int ycount, string name)
        {
            string[] vbsraw = File.ReadAllLines(script_template);
            vbsraw[0] = String.Format(vbsraw[0], ycount);
            vbsraw[1] = String.Format(vbsraw[1], xcount);
            vbsraw[2] = String.Format(vbsraw[2], name);
            vbsraw[3] = String.Format(vbsraw[3], bounds.Xmin);
            vbsraw[4] = String.Format(vbsraw[4], bounds.Xmax);
            vbsraw[5] = String.Format(vbsraw[5], bounds.Ymin);
            vbsraw[6] = String.Format(vbsraw[6], bounds.Ymax);
            vbsraw[14] = String.Format(vbsraw[14], (2 + xcount * ycount));
            vbsraw[25] = String.Format(vbsraw[25], vb_syntax_eval);
            vbsraw[50] = String.Format(vbsraw[50], (2 + xcount * ycount));
            string newname = @"C:\Users\Will\Desktop\Folders\MATH435\repo\visual-basic-templates\TEMPORARY.vbs";
            File.WriteAllLines(newname, vbsraw);
            run_script(newname);
        }
        private static void run_script(string scriptname_absolute)
        {
            Process vb_script = Process.Start(CSCRIPT, scriptname_absolute);
            while (!vb_script.HasExited)
            {
                Thread.Sleep(10);
            }
            File.Delete(scriptname_absolute);
        }
        public static void quickPlot(string vb_eval, string title, RBounds2D bounds)
        {
            int n = 100;
            GenerateFunction(vb_eval, bounds, n, n, title);
            NCGrid_Distribution p = NCGrid_Distribution.from_file(title, false);
            File.Delete(Paths.DistributionRepo + "\\" + title + ".dist");
            DistributionSketchSettings S = DistributionSketchSettings.Fancy();
            S.HasHeatmap = true;
            S.HasInfo = false;
            S.Mode = SketchMode.DISTRIBUTION_ONLY;
            DistributionSketch2D sk = new DistributionSketch2D(p, S);
            sk.CreateSketch(false);
            sk.SaveImage(title, false);
        }
        public static void quickPlot(string vb_eval, string title)
        {
            quickPlot(vb_eval, title, RBounds2D.Square(10));
        }
        public static NCGrid_Distribution GenerateFunctionToGrid(string vb_eval, NCGrid_Distribution discretization)
        {
            string absolute_temp_path = Paths.DistributionRepo + "\\temp.dist";
            GenerateFunction(vb_eval, discretization, "temp");
            NCGrid_Distribution output = NCGrid_Distribution.from_file("temp", false);
            File.Delete(absolute_temp_path);
            return output;
        }
        public static void GenerateFunction(string vb_syntax_eval, NCGrid_Distribution discretization, string name)
        {
            discretization.WriteToFile(@"C:\Users\Will\Desktop\Folders\MATH435\repo\visual-basic-templates\TEMPORARY.dist", false, true, true);
            string[] vbsraw = File.ReadAllLines(script_disc_template);
            int fileconst = discretization.Xcount * discretization.Ycount + 2;
            vbsraw[0] = String.Format(vbsraw[0], discretization.Ycount);
            vbsraw[1] = String.Format(vbsraw[1], discretization.Xcount);
            vbsraw[16] = String.Format(vbsraw[16], fileconst);
            vbsraw[13] = String.Format(vbsraw[13], fileconst);
            vbsraw[2] = String.Format(vbsraw[2], name);
            vbsraw[28] = String.Format(vbsraw[28], vb_syntax_eval);
            string newname = @"C:\Users\Will\Desktop\Folders\MATH435\repo\visual-basic-templates\TEMPORARY.vbs";
            File.WriteAllLines(newname, vbsraw);
            run_script(newname);
            File.Delete(@"C:\Users\Will\Desktop\Folders\MATH435\repo\visual-basic-templates\TEMPORARY.dist");
        }
    }
}
