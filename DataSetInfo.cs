/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace m_435_NCAMR
{
    class DataSetInfo
    {
        private int distrib_count;
        private double universal_max, universal_min;
        public int DistributionCount
        {
            get { return distrib_count; }
        }
        public double UniversalMax
        {
            get { return universal_max; }
        }
        public double UniversalMin
        {
            get { return universal_min; }
        }
        public DataSetInfo(string[] fullfilenames, bool enable_console_output)
        {
            DateTime then = DateTime.Now;
            universal_min = double.MaxValue;
            universal_max = double.MinValue;
            List<string> all_dist = new List<string>();
            foreach (string p in fullfilenames)
            {
                if (p.EndsWith(".dist"))
                {
                    all_dist.Add(p);
                }
            }
            distrib_count = all_dist.Count;
            for (int i = 0; i < distrib_count; i++)
            {
                using (StreamReader r = new StreamReader(all_dist[i]))
                {
                    string first_line = r.ReadLine();
                    string[] spt = first_line.Split(':');
                    double tempmax, tempmin;
                    if (!(double.TryParse(spt[6], out tempmax)&& double.TryParse(spt[8], out tempmin))) { throw new Exception("File improperly formatted."); }
                    if (tempmax > universal_max) { universal_max = tempmax; }
                    if (tempmin < universal_min) { universal_min = tempmin; }
                }
            }
            DateTime now = DateTime.Now;
            if (enable_console_output)
            {
                Console.WriteLine("DataSetInfo generated in " + (now - then).Milliseconds.ToString() + " milliseconds");
                Console.WriteLine(this.ToString());
            }
        }
        public override string ToString()
        {
            //incomplete: will need interpolation coefficients
            return "dist_count:" + distrib_count + ":max_value:" + universal_max.ToString() + ":min_value:" + universal_min.ToString();
        }
        public void write_to_info_file(string fullpath)
        {
            string[] contents = { this.ToString() };
            File.WriteAllLines(fullpath, contents);
        }
    }
}
