/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace m_435_NCAMR
{
    static class RepoManagement
    {
        public static void ClearRepo(params string[] abs_repo_dists)
        {
            List<string[]> arrs = new List<string[]>();
            foreach (string i in abs_repo_dists)
            {
                arrs.Add(Directory.GetFiles(i));
            }
            foreach (string[] repnames in arrs)
            {
                foreach (string y in repnames)
                {
                    File.Delete(y);
                }
            }
        }
    }
}
