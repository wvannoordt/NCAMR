/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Drawing;
using AForge.Video.VFW;
using _3DSimple;


namespace m_435_NCAMR
{
    class FollowedDataSet : DataSet
    {
        public FollowedDataSet(string title, bool using_full_path, bool enable_console_output) : base(title, using_full_path, enable_console_output)
        {

        }
        public new static FollowedDataSet FirstAvailable(bool enable_console_output)
        {
            string[] allnames = Directory.GetDirectories(default_source);
            if (allnames.Length != 0)
            {
                return new FollowedDataSet(allnames[0], true, enable_console_output);
            }
            else
            {
                throw new Exception("No datasets available.");
            }
        }
        public override void create_animation(string output_title, bool enable_console_output, bool using_full_path, bool allow_overwrite, bool use_floating_color_scale, DistributionSketchSettings S)
        {
            string absolute_video_output_path = video_default_repo + "\\" + output_title + ".avi";
            if (using_full_path) { absolute_video_output_path = output_title; }
            if (File.Exists(absolute_video_output_path) && !allow_overwrite) { throw new Exception("Error: File exists and overwrite permission not granted."); }
            Directory.CreateDirectory(AbsoluteDataSetLocation + "\\images");
            string[] allfiles = Directory.GetFiles(absolute_set_location);
            List<string> distrib_files = new List<string>();
            int q = 0;
            foreach (string i in allfiles)
            {
                if (i.EndsWith(".dist"))
                {
                    distrib_files.Add(i);
                    NCGrid_Distribution n = NCGrid_Distribution.from_file(i, true);
                    NCGrid_Distribution following_grid = new NCGrid_Distribution(n.Bounds, n.Xcount, n.Ycount);
                    int border = 1;
                    for (int ii = border; ii < following_grid.Xcount- border; ii++)
                    {
                        for (int jj = border; jj < following_grid.Ycount- border; jj++)
                        {
                            double deltax = (n.Xmax - n.Xmin) / (n.Xcount - 1);
                            double deltay = (n.Ymax - n.Ymin) / (n.Ycount - 1);
                            double[] stencil = 
                            {
                                n[ii-1,jj-1].Value,
                                n[ii,jj-1].Value,
                                n[ii+1,jj-1].Value,
                                n[ii-1,jj].Value,
                                n[ii,jj].Value,
                                n[ii+1,jj].Value,
                                n[ii-1,jj+1].Value,
                                n[ii,jj+1].Value,
                                n[ii+1,jj+1].Value,
                            };
                            double ux = (stencil[5] - stencil[3]) / (2 * deltax);
                            double uy = (stencil[7] - stencil[1]) / (2 * deltay);
                            double uxx = (stencil[5]+stencil[3] - (2*stencil[4])) / (deltax * deltax);
                            double uyy = (stencil[7] + stencil[1] - (2 * stencil[4])) / (deltay * deltay);
                            double uxy = (stencil[8]+stencil[0]-stencil[6]-stencil[2])/(4*deltax*deltay);
                            double dx = 0.06 * ((uxx*ux)+(uxy*uy));
                            double dy = 0.06 * ((uyy*uy)+(ux*uxy));
                            Vector3 move = new Vector3(dx, dy, 0);
                            while ((following_grid[ii, jj] + move).X > n.Bounds.Xmax || (following_grid[ii, jj] + move).X < n.Bounds.Xmin || (following_grid[ii, jj] + move).Y > n.Bounds.Ymax || (following_grid[ii, jj] + move).Y < n.Bounds.Ymin)
                            {
                                dx = 0.5 * dx;
                                dy = 0.5 * dy;
                                move = new Vector3(dx, dy, 0);
                            }
                            following_grid[ii, jj] = following_grid[ii, jj] + move;
                        }
                    }
                    HybridDistributionSketch2D sk = new HybridDistributionSketch2D(n, S);
                    if (!use_floating_color_scale)
                    {
                        sk.override_colormap_limits(info);
                    }
                    sk.set_superimposed_grid(following_grid);
                    sk.CreateSketch(true);
                    sk.SaveImage(AbsoluteDataSetLocation + "\\images\\render" + bufferint(q++, 5) + ".bmp", true);
                }
            }
            AVIWriter V = new AVIWriter("wmv3");
            string[] filenames = Directory.GetFiles(AbsoluteDataSetLocation + "\\images");
            List<string> bmpnames = new List<string>();
            foreach (string i in filenames)
            {
                if (i.EndsWith(".bmp")) { bmpnames.Add(i); }
            }
            if (bmpnames.Count != 0)
            {
                Bitmap first = (Bitmap)Bitmap.FromFile(bmpnames[0]);
                V.Open(absolute_video_output_path, first.Width, first.Height);
                V.FrameRate = 50;
                V.AddFrame(first);
                int ct = bmpnames.Count;
                for (int i = 1; i < bmpnames.Count; i++)
                {
                    Bitmap butt = (Bitmap)Bitmap.FromFile(bmpnames[i]);
                    if (enable_console_output) { Console.WriteLine(i.ToString() + " of " + ct.ToString() + " frames stacked" + "(" + (i * 100 / ct).ToString() + "%)"); }
                    V.AddFrame(butt);
                    butt.Dispose();
                }
                V.Close();
                first.Dispose();
                if (enable_console_output) { Console.WriteLine("AVI successfully created on " + DateTime.Now.ToString() + " in directory " + absolute_video_output_path); }
                string[] imagenames = Directory.GetFiles(AbsoluteDataSetLocation + "\\images");
                foreach (string g in imagenames)
                {
                    File.Delete(g);
                }
                Directory.Delete(AbsoluteDataSetLocation + "\\images");
            }
            else
            {
                Console.WriteLine("No bitmap images found in the current directory.");
            }
        }
    }
}
