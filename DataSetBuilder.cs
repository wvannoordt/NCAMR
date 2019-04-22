/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO.Compression;
using System.IO;

namespace m_435_NCAMR
{
    static class DataSetBuilder
    {
        private static string info_file_tite = "info.inf";
        private static string default_repo = @"C:\Users\Will\Desktop\Folders\MATH435\repo\datasets";
        private static string default_source = @"C:\Users\Will\Desktop\Folders\MATH435\repo\distributions";
        public static void build_from_default_repo(string output_title, bool using_full_path, bool enable_console_output)
        {
            string absolute_source_path = default_source;
            string absolute_output_path = default_repo + "\\" + output_title;
            if (using_full_path) { absolute_source_path = output_title; }
            string[] all_absolute_names = Directory.GetFiles(absolute_source_path);
            //will need interpolation/regression code.
            DataSetInfo info = new DataSetInfo(all_absolute_names, enable_console_output);
            Directory.CreateDirectory(absolute_output_path);
            int n = 0;
            foreach(string i in all_absolute_names)
            {
                if (enable_console_output) { Console.WriteLine("Processing "+ (++n).ToString() + " of " + all_absolute_names.Length.ToString() + " files." ); }
                File.Copy(i, absolute_output_path + "\\" + rip_file_title(i));
                File.Delete(i);
            }
            info.write_to_info_file(absolute_output_path + "\\" + info_file_tite);
        }
        private static string rip_file_title(string fullpath)
        {
            return fullpath.Split('\\').Last();
        }
        public static void build_by_keyword(string keyword, string path)
        {

        }
    }
}
