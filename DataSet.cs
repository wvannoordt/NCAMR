/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Drawing;
using AForge.Video.VFW;

namespace m_435_NCAMR
{
    class DataSet
    {
        protected string absolute_set_location;
        protected const string video_default_repo = @"C:\Users\Will\Desktop\Folders\MATH435\repo\animations";
        protected const string default_source = @"C:\Users\Will\Desktop\Folders\MATH435\repo\datasets";
        protected DataSetInfo info;
        public string AbsoluteDataSetLocation
        {
            get { return absolute_set_location; }
            set { absolute_set_location = value; }
        }
        protected static string bufferint(int n, int size)
        {
            string l = string.Empty;
            for (int i = 0; i < size + 1 - n.ToString().Length; i++)
            {
                l += "0";
            }
            l += n.ToString();
            return l;
        }
        public virtual void create_animation(string output_title, bool enable_console_output, bool using_full_path, bool allow_overwrite, bool use_floating_color_scale, DistributionSketchSettings S)
        {
            string absolute_video_output_path = video_default_repo + "\\" + output_title + ".avi";
            if (using_full_path) { absolute_video_output_path = output_title; }
            if (File.Exists(absolute_video_output_path)&&!allow_overwrite) { throw new Exception("Error: File exists and overwrite permission not granted."); }
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
                    DistributionSketch2D sk = new DistributionSketch2D(n, S);
                    if (!use_floating_color_scale)
                    {
                        sk.override_colormap_limits(info);
                    }
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
        public DataSet(string title, bool using_full_path, bool enable_console_output)
        {
            absolute_set_location = default_source + "\\" + title;
            if (using_full_path) { absolute_set_location = title; }
            info = new DataSetInfo(Directory.GetFiles(absolute_set_location), enable_console_output);
        }
        public static DataSet FirstAvailable(bool enable_console_output)
        {
            string[] allnames = Directory.GetDirectories(default_source);
            if (allnames.Length != 0)
            {
                return new DataSet(allnames[0], true, enable_console_output);
            }
            else
            {
                throw new Exception("No datasets available.");
            }
        }
    }
}
