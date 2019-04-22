/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

namespace m_435_NCAMR
{
    static class OpenImage
    {
        //needs work
        static string prg = @"C:\Program Files (x86)\Windows Live\Photo Gallery\WLXPhotoGallery";
        static string img_repo = @"C:\Users\Will\Desktop\Folders\MATH435\repo\images";
        public static void Show(string title)
        {
            string fullpath = img_repo + "\\" + title + ".bmp";
            Process.Start(prg, fullpath);
        }
        public static void Show(string path, bool using_full)
        {
            Process.Start(prg, path);
        }
    }
}
