/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;
using _3DSimple;
using ColorMapping;

namespace m_435_NCAMR
{
    class DistributionSketch3DBasic
    {
        protected const string default_repo = @"C:\Users\Will\Desktop\Folders\MATH435\repo\images";
        private NCGrid_Distribution subject;
        Bitmap output;
        public NCGrid_Distribution Subject
        {
            get { return subject; }
            set { subject = value; }
        }
        public DistributionSketch3DBasic(NCGrid_Distribution sub)
        {
            subject = sub;
            output = new Bitmap(2000, 2000);
        }
        public void SaveImage(string title, bool using_full_path)
        {
            string path = default_repo + "\\" + title + ".bmp";
            if (using_full_path) { path = title; }
            output.Save(path);
        }
        public void GenerateSketch()
        {
            using (var g = Graphics.FromImage(output))
            {
                Point og = new Point(output.Width / 2, output.Height / 2);
                Point3 light = new Point3(subject.Bounds.Xmin+0.5*(subject.Bounds.Xmax-subject.Bounds.Xmin), subject.Bounds.Ymin + 0.5 * (subject.Bounds.Ymax - subject.Bounds.Ymin), 1.5*subject.MaxValue);
                double camZ = 45;
                double camXY = 45;
                double luster = 0.9;
                double zoom = 100;
                double dx = (subject.Bounds.Xmax - subject.Bounds.Xmin) / (subject.Xcount - 1);
                double dy = (subject.Bounds.Ymax - subject.Bounds.Ymin) / (subject.Ycount - 1);
                g.Clear(Color.DarkGray);
                for (int i = 0; i < subject.Xcount-1; i++)
                {
                    for (int j = 0; j < subject.Ycount-1; j++)
                    {
                        double k = Math.PI / 180;
                        Point3 P00 = subject[i, j].toPoint();
                        Point3 P01 = subject[i, j + 1].toPoint();
                        Point3 P10 = subject[i + 1, j].toPoint();
                        Point3 P11 = subject[i + 1, j + 1].toPoint();
                        Vector3 V1 = new Vector3(P00, P11);
                        Vector3 V2 = new Vector3(P10, P01);
                        float cosZ = (float)Math.Cos(camZ * k);
                        float cosXY = (float)Math.Cos(camXY * k);
                        float sinZ = (float)Math.Sin(camZ * k);
                        float sinXY = (float)Math.Sin(camXY * k);
                        Vector3 vw = new Vector3(cosZ * cosXY, cosZ * sinXY, sinZ);
                        Vector3 Normal = V1 % V2;
                        double zavg = 0.25 * (P00.Z + P10.Z + P01.Z + P11.Z);
                        Point a00 = transform(P00, camZ, camXY, zoom, og);
                        Point a01 = transform(P01, camZ, camXY, zoom, og);
                        Point a10 = transform(P10, camZ, camXY, zoom, og);
                        Point a11 = transform(P11, camZ, camXY, zoom, og);
                        Point[] pts = new Point[4];
                        pts[0] = a00;
                        pts[1] = a01;
                        pts[2] = a11;
                        pts[3] = a10;
                        Vector3 liteVector = new Vector3(P00, light);
                        double ux2 = ((P10.Z - P00.Z) / dx) * ((P10.Z - P00.Z) / dx);
                        double uy2 = ((P01.Z - P00.Z) / dy) * ((P01.Z - P00.Z) / dy);
                        double strain = 0.5 * ((Math.Sqrt(1 + ux2) - 1) + (Math.Sqrt(1 + uy2) - 1));
                        SolidBrush B = new SolidBrush(ShadeSingle(Color.Red, liteVector, Normal, camXY, camZ, luster));
                        g.FillPolygon(B, pts);
                    }
                }
                g.Dispose();
            }
        }
        Point transform(Point3 objectLocation, double cameraZ, double cameraXY, double zoom, Point form_origin)
        {
            cameraXY = 90 - cameraXY;
            double conv = Math.PI / 180;
            double cos_z = Math.Cos(cameraZ * conv);
            double sin_z = Math.Sin(cameraZ * conv);
            double cos_xy = Math.Cos(cameraXY * conv);
            double sin_xy = Math.Sin(cameraXY * conv);
            int px = (int)(((objectLocation.Y * sin_xy) - (objectLocation.X * cos_xy)) * zoom);
            int py = (int)(((objectLocation.Z * cos_z) - (objectLocation.X * sin_xy * sin_z) - (objectLocation.Y * sin_z * cos_xy)) * zoom);
            return new Point(form_origin.X + px, form_origin.Y - py);
        }
        private Color ShadeSingle(Color C, Vector3 light, Vector3 Normal, double camXY, double camZ, double luster)
        {
            double lus_thresh = luster * 0.1;
            if (lus_thresh > 1) { lus_thresh = 1; }
            Vector3 N = Normal.unit();
            if (N.K < 0) { N = -1 * N; }
            Vector3 L = light.unit();
            double k = Math.PI / 180;
            float cosZ = (float)Math.Cos(camZ * k);
            float cosXY = (float)Math.Cos(camXY * k);
            float sinZ = (float)Math.Sin(camZ * k);
            float sinXY = (float)Math.Sin(camXY * k);
            Vector3 vw = new Vector3(cosZ * cosXY, cosZ * sinXY, sinZ);
            double mu = ((2 * (float)(L * N) * N) - L) * vw;
            if (mu < lus_thresh) { return ColorMap.MapGradient(Color.FromArgb(255, 10, 10, 16), C, mu, -1, lus_thresh); }
            int r = (int)(lus_thresh * 255);
            return ColorMap.MapGradient(C, Color.FromArgb(255, r, r, r), mu, lus_thresh, 1);
        }
    }
}
