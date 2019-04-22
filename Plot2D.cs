/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;
using ColorMapping;
using System.IO;

namespace m_435_NCAMR
{
    class Plot2D
    {
        protected const string default_repo = @"C:\Users\Will\Desktop\Folders\MATH435\repo\images";
        protected Plot2DSettings settings;
        protected float override_max;
        protected float override_min;
        protected bool use_overriden_limits = false;
        protected List<Color> HReferenceColors = new List<Color>();
        protected List<double> HReferenceValues = new List<double>();
        protected List<Color> VReferenceColors = new List<Color>();
        protected List<double> VReferenceValues = new List<double>();
        protected List<double[]> xsurplus = new List<double[]>();
        protected List<double[]> ysurplus = new List<double[]>();
        protected List<Color> colorsurplus = new List<Color>();
        protected string[] overwrite_x_labels = null;
        protected string[] overwrite_y_labels = null;

        //parameters
        protected float MAIN_SUBJECT_LOWER_LEFT_X_2D;
        protected float MAIN_SUBJECT_LOWER_LEFT_Y_2D;
        protected float MAIN_SUBJECT_UPPER_RIGHT_X_2D;
        protected float MAIN_SUBJECT_UPPER_RIGHT_Y_2D;

        //constants
        protected const float FIGURETITLE_Y = 0.04f;
        protected const float FIGURETITLE_X = 0.1f;
        protected const float SQUARE_BORDER_THIN = 0.07f;
        protected const float SQUARE_BORDER_THICK = 0.130f;
        protected const float GRIDLINE_THICKNESS = 3.3f;
        protected const float X_AXIS_TITLE_X = 0.5f;
        protected const float Y_AXIS_TITLE_Y = 0.5f;
        protected const float Y_AXIS_TITLE_X = 0.1f * SQUARE_BORDER_THICK;
        protected const float X_AXIS_TITLE_Y = 0.5f * SQUARE_BORDER_THICK;
        protected const float LABEL_OFFSET_Y = -0.056f;

        protected const float HEATBOX_BORDER_LEFT = 0.84f;
        protected const float HEATBOX_BORDER_RIGHT = 0.97f;
        protected static float HEATBOX_BORDER_BOTTOM = 0.35f;
        protected float HEATBOX_BORDER_TOP = 1 - HEATBOX_BORDER_BOTTOM;
        protected const float HEATBOX_LL_X = 0.85f;
        protected float HEATBOX_UR_X = 0.96f;
        protected static float HEATBOX_LL_Y = 0.36f;
        protected float HEATBOX_UR_Y = 1 - HEATBOX_LL_Y;

        protected float HEATMAP_LABELS_X = 0.91f;
        protected static float HEATMAP_LABELS_YU = 0.33f;
        protected static float HEATMAP_LABELS_YL = 0.64f;

        protected float X_LABEL_BASE_Y = 0.94f * SQUARE_BORDER_THICK;
        protected float X_LABEL_BASE_X = 1.0f * SQUARE_BORDER_THICK;
        protected float Y_LABEL_BASE_Y = 1.0f * SQUARE_BORDER_THICK;
        protected float Y_LABEL_BASE_X = 0.44f * SQUARE_BORDER_THICK;
        protected float VALUE_TEXT_PROPORTION = 40f / 2000f;
        protected float TITLE_TEXT_PROPORTION = 65f / 2000f;
        protected const float THIN_PEN = 0.8f;
        protected const float THICK_PEN = 2.4f;
        protected Pen THICK_GRID_PEN = new Pen(Color.Black, THICK_PEN);
        protected Pen THIN_GRID_PEN = new Pen(Color.Black, THIN_PEN);
        protected Color BACKGROUND_COLOR = Color.White;
        protected Color TEXTCOLOR = Color.Black;
        protected Color DEFAULTCURVECOLOR = Color.Black;
        RBounds2D sketch_bounds;
        double[] xvals;
        double[] yvals;

        //calculated relatives
        protected float grid_xmin;
        protected float grid_xmax;
        protected float grid_ymin;
        protected float grid_ymax;
        protected float pix_per_unit_x;
        protected float pix_per_unit_y;
        protected static Color[] default_colors = {Color.Blue, Color.Red, Color.Orange, Color.Green, Color.Green, Color.Yellow, Color.Gray, Color.DarkGoldenrod, Color.DarkGreen};
        protected Bitmap canvas;
        protected int x_pixel_count, y_pixel_count;
        public Bitmap Canvas
        {
            get { return canvas; }
        }
        public Plot2D(string function_VBSYNTAX, Plot2DSettings S)
        {
            settings = S;
            throw new Exception("Error: Function not yet implemented");
            canvas = new Bitmap(S.ImageWidth, S.ImageHeight);
            x_pixel_count = canvas.Width;
            y_pixel_count = canvas.Height;
            get_positions();

        }
        public Plot2D(double[] x, double[] y, Plot2DSettings S)
        {
            settings = S;
            if (x.Length != y.Length) { throw new Exception("Error: Array dimensions are inconsistent."); }
            xvals = x;
            yvals = y;
            canvas = new Bitmap(S.ImageWidth, S.ImageHeight);
            x_pixel_count = canvas.Width;
            y_pixel_count = canvas.Height;
            double ymax = double.NegativeInfinity;
            double ymin = double.PositiveInfinity;
            double xmax = double.NegativeInfinity;
            double xmin = double.PositiveInfinity;
            for (int i = 0; i < x.Length; i++)
            {
                if (x[i] > xmax) { xmax = x[i]; }
                if (x[i] < xmin) { xmin = x[i]; }
                if (y[i] > ymax) { ymax = y[i]; }
                if (y[i] < ymin) { ymin = y[i]; }
            }
            double deltax = xmax - xmin;
            double deltay = ymax - ymin;
            double xbar = (xmax + xmin) * 0.5;
            double ybar = (ymax + ymin) * 0.5;
            double xscl = 1.0;
            double yscl = 1.0;
            double newdeltax = xscl * deltax;
            double newdeltay = yscl * deltay;
            xmax = xbar + 0.5 * newdeltax;
            xmin = xbar - 0.5 * newdeltax;
            ymax = ybar + 0.5 * newdeltay;
            ymin = ybar - 0.5 * newdeltay;
            sketch_bounds = new RBounds2D(xmin, xmax, ymin, ymax);
            get_positions();
        }
        public void ScaleCurveY(double scale_factor)
        {
            for (int i = 0; i < yvals.Length; i++)
            {
                yvals[i] = yvals[i] * scale_factor;
                sketch_bounds.Ymax *= scale_factor;
                sketch_bounds.Ymin *= scale_factor;
            }
        }
        public virtual void CreateSketch(bool allow_console_output)
        {
            using (Graphics drawer = Graphics.FromImage(canvas))
            {
                drawer.Clear(BACKGROUND_COLOR);
                draw_references(drawer);
                draw_curve(drawer);
                draw_fig_title(drawer);
                if (settings.HasAxisTitles) { draw_axis_titles(drawer); }
                if (settings.HasAxisValues) { draw_axis_values(drawer); }
                if (settings.HasGridlines) { draw_gridlines(drawer); }
            }
        }
        public void AddHorizontalReference(double val, Color C)
        {
            HReferenceValues.Add(val);
            HReferenceColors.Add(C);
        }
        public void AddVerticalReference(double val, Color C)
        {
            VReferenceValues.Add(val);
            VReferenceColors.Add(C);
        }
        protected void draw_references(Graphics drawer)
        {
            for (int i = 0; i < HReferenceValues.Count; i++)
            {
                if (is_between(HReferenceValues[i], sketch_bounds.Ymin, sketch_bounds.Ymax))
                {
                    Pen currentpen = new Pen(HReferenceColors[i], 6.0f);
                    drawer.DrawLine(currentpen, mapXY(sketch_bounds.Xmin, HReferenceValues[i]), mapXY(sketch_bounds.Xmax, HReferenceValues[i]));
                }
            }
            for (int i = 0; i < VReferenceValues.Count; i++)
            {
                if (is_between(VReferenceValues[i], sketch_bounds.Xmin, sketch_bounds.Xmax))
                {
                    Pen currentpen = new Pen(VReferenceColors[i], 6.0f);
                    drawer.DrawLine(currentpen, mapXY(VReferenceValues[i], sketch_bounds.Ymin), mapXY(VReferenceValues[i], sketch_bounds.Ymax));
                }
            }
        }
        protected bool is_between(double val, double lower, double upper)
        {
            return (val < upper) && (val > lower);
        }
        protected void draw_curve(Graphics drawer)
        {
            Pen curvepen = new Pen(Color.Black, 4.5f);
            List<PointF> curve = new List<PointF>();
            for (int i = 0; i < xvals.Length; i++)
            {
                curve.Add(mapXY(xvals[i], yvals[i]));
            }
            for (int i = 0; i < xsurplus.Count; i++)
            {
                List<PointF> curcurve = new List<PointF>();
                Pen curcurvepen = new Pen(colorsurplus[i], 4.5f);
                for (int j = 0; j < xsurplus[i].Length; j++)
                {
                    curcurve.Add(mapXY(xsurplus[i][j], ysurplus[i][j]));
                }
                drawer.DrawLines(curcurvepen, curcurve.ToArray());
            }
            drawer.DrawLines(curvepen, curve.ToArray());
        }
        public void AddCurve(double[] xs, double[] ys, Color C)
        {
            xsurplus.Add(xs);
            ysurplus.Add(ys);
            colorsurplus.Add(C);
        }
        protected void get_positions()
        {
            if (settings.HasAxisValues || settings.HasAxisTitles)
            {
                MAIN_SUBJECT_LOWER_LEFT_X_2D = SQUARE_BORDER_THICK + settings.WidthOffset;
                MAIN_SUBJECT_LOWER_LEFT_Y_2D = SQUARE_BORDER_THICK;
                MAIN_SUBJECT_UPPER_RIGHT_X_2D = 1 - SQUARE_BORDER_THICK + settings.WidthOffset;
                MAIN_SUBJECT_UPPER_RIGHT_Y_2D = 1 - SQUARE_BORDER_THICK;
            }
            else
            {
                MAIN_SUBJECT_LOWER_LEFT_X_2D = SQUARE_BORDER_THIN;
                MAIN_SUBJECT_LOWER_LEFT_Y_2D = SQUARE_BORDER_THIN;
                MAIN_SUBJECT_UPPER_RIGHT_X_2D = 1 - SQUARE_BORDER_THIN;
                MAIN_SUBJECT_UPPER_RIGHT_Y_2D = 1 - SQUARE_BORDER_THIN;
            }

            grid_xmin = MAIN_SUBJECT_LOWER_LEFT_X_2D * (float)x_pixel_count;
            grid_xmax = MAIN_SUBJECT_UPPER_RIGHT_X_2D * (float)x_pixel_count;
            grid_ymin = MAIN_SUBJECT_LOWER_LEFT_Y_2D * (float)y_pixel_count;
            grid_ymax = MAIN_SUBJECT_UPPER_RIGHT_Y_2D * (float)y_pixel_count;
            pix_per_unit_x = (grid_xmax - grid_xmin) / (float)(sketch_bounds.Xmax - sketch_bounds.Xmin);
            pix_per_unit_y = (grid_ymax - grid_ymin) / (float)(sketch_bounds.Ymax - sketch_bounds.Ymin);
        }
        public void OverrideLabels(string[] xlabels, string[] ylabels)
        {
            overwrite_x_labels = xlabels;
            overwrite_y_labels = ylabels;
        }
        protected void draw_axis_values(Graphics drawer)
        {
            Brush textbrush = new SolidBrush(TEXTCOLOR);
            float xlabel_abs_pos_x = ((float)x_pixel_count * (X_LABEL_BASE_X+settings.WidthOffset));
            float xlabel_abs_increment = (((1 - 2 * (Y_LABEL_BASE_Y)) * (float)x_pixel_count) / (float)(settings.XLabelCount - 1));
            float ylabel_abs_pos_y = ((float)y_pixel_count * (Y_LABEL_BASE_X));
            float ylabel_abs_increment = (((1 - 2 * Y_LABEL_BASE_Y) * (float)y_pixel_count) / (float)(settings.YLabelCount - 1));
            float xlabel_abs_pos_y = (X_LABEL_BASE_Y * (float)x_pixel_count);
            float ylabel_abs_pos_x = ((Y_LABEL_BASE_X - LABEL_OFFSET_Y) * (float)y_pixel_count);
            if (overwrite_x_labels != null && overwrite_y_labels != null)
            {
                double delta_x = (sketch_bounds.Xmax - sketch_bounds.Xmin) / (overwrite_x_labels.Length - 1);
                double delta_y = (sketch_bounds.Ymax - sketch_bounds.Ymin) / (overwrite_y_labels.Length - 1);
                xlabel_abs_increment = (((1 - 2 * (Y_LABEL_BASE_Y)) * (float)x_pixel_count) / (float)(overwrite_x_labels.Length - 1));
                ylabel_abs_increment = (((1 - 2 * Y_LABEL_BASE_Y) * (float)y_pixel_count) / (float)(overwrite_y_labels.Length - 1));
                for (int i = 0; i < overwrite_x_labels.Length; i++)
                {
                    drawer.DrawString(overwrite_x_labels[i], new Font("arial", settings.FontSize), textbrush, new PointF(xlabel_abs_pos_x + i * xlabel_abs_increment, settings.ImageHeight - xlabel_abs_pos_y));
                }
                for (int j = 0; j < overwrite_y_labels.Length; j++)
                {
                    drawer.DrawString(overwrite_y_labels[j], new Font("arial", settings.FontSize), textbrush, new PointF(ylabel_abs_pos_y, ylabel_abs_pos_x + j * ylabel_abs_increment));
                }
            }
            else
            {
                double delta_x = (sketch_bounds.Xmax - sketch_bounds.Xmin) / (settings.XLabelCount - 1);
                double delta_y = (sketch_bounds.Ymax - sketch_bounds.Ymin) / (settings.YLabelCount - 1);
                for (int i = 0; i < settings.XLabelCount; i++)
                {
                    drawer.DrawString(truncate_string((sketch_bounds.Xmin + i * delta_x).ToString(), 4), new Font("arial", settings.FontSize), textbrush, new PointF(xlabel_abs_pos_x + i * xlabel_abs_increment, settings.ImageHeight - xlabel_abs_pos_y));
                }
                for (int j = 0; j < settings.YLabelCount; j++)
                {
                    drawer.DrawString(truncate_string((sketch_bounds.Ymax - j * delta_y).ToString(), 4), new Font("arial", settings.FontSize), textbrush, new PointF(ylabel_abs_pos_y, ylabel_abs_pos_x + j * ylabel_abs_increment));
                }
            }
            
        }
        protected void draw_gridlines(Graphics drawer)
        {
            Pen gridlinepen = new Pen(Color.FromArgb(130, 0, 0, 0), GRIDLINE_THICKNESS);
            NCAMRNode[,] gridpoints = new NCAMRNode[settings.XGridlineCount, settings.YGridlineCount];
            double coordinate_dx = (sketch_bounds.Xmax - sketch_bounds.Xmin) / (settings.XGridlineCount - 1);
            double coordinate_dy = (sketch_bounds.Ymax - sketch_bounds.Ymin) / (settings.YGridlineCount - 1);
            for (int i = 0; i < settings.XGridlineCount; i++)
            {
                for (int j = 0; j < settings.YGridlineCount; j++)
                {
                    gridpoints[i, j] = new NCAMRNode(sketch_bounds.Xmin + i * coordinate_dx, sketch_bounds.Ymin + j * coordinate_dy, 0);
                }
            }
            for (int i = 0; i < settings.XGridlineCount; i++)
            {
                for (int j = 0; j < settings.YGridlineCount - 1; j++)
                {
                    //thick grid pen temporary
                    drawer.DrawLine(gridlinepen, mapGridNode(gridpoints[i, j]), mapGridNode(gridpoints[i, j + 1]));
                }
            }
            for (int i = 0; i < settings.YGridlineCount; i++)
            {
                for (int j = 0; j < settings.XGridlineCount - 1; j++)
                {
                    //thick grid pen temporary
                    drawer.DrawLine(gridlinepen, mapGridNode(gridpoints[j + 1, i]), mapGridNode(gridpoints[j, i]));
                }
            }
        }
        protected void draw_fig_title(Graphics drawer)
        {
            float absx = settings.ImageWidth * FIGURETITLE_X;
            float absy = settings.ImageHeight * FIGURETITLE_Y;
            Brush textbrush = new SolidBrush(TEXTCOLOR);
            drawer.DrawString(settings.FigureTitle, new Font("arial", 1.9f * VALUE_TEXT_PROPORTION * settings.ImageHeight), textbrush, new PointF(absx, absy));
        }
        protected void draw_axis_titles(Graphics drawer)
        {
            float x_title_abs_pos_x = X_AXIS_TITLE_X * settings.ImageWidth;
            float x_title_abs_pos_y = X_AXIS_TITLE_Y * settings.ImageWidth;
            float y_title_abs_pos_x = Y_AXIS_TITLE_X * settings.ImageHeight;
            float y_title_abs_pos_y = Y_AXIS_TITLE_Y * settings.ImageHeight;
            Brush textbrush = new SolidBrush(TEXTCOLOR);
            drawer.DrawString(settings.HorizontalTitle, new Font("arial", TITLE_TEXT_PROPORTION * settings.ImageHeight), textbrush, new PointF(x_title_abs_pos_x, settings.ImageHeight - x_title_abs_pos_y));
            drawer.DrawString(settings.VerticalTitle, new Font("arial", TITLE_TEXT_PROPORTION * settings.ImageHeight), textbrush, new PointF(y_title_abs_pos_x, y_title_abs_pos_y));
        }
        protected PointF mapXY(double x, double y)
        {
            return new PointF(grid_xmin + pix_per_unit_x * (float)(x - sketch_bounds.Xmin), grid_ymin + pix_per_unit_y * (float)((sketch_bounds.Ymax - y)));
        }
        protected PointF mapGridNode(NCAMRNode toMap)
        {
            return new PointF(grid_xmin + pix_per_unit_x * (float)(toMap.X - sketch_bounds.Xmin), grid_ymin + pix_per_unit_y * (float)((sketch_bounds.Ymax - toMap.Y)));
        }
        public static Plot2D ReadCSV(string path, Plot2DSettings S, params int[] yindeces)
        {
            string[] filestuff = File.ReadAllLines(path);
            int stuffcount = filestuff[0].Split(',').Length-1;
            int[] index = yindeces;
            if (yindeces.Length == 0)
            {
                index = new int[stuffcount];
                for (int i = 0; i < stuffcount; i++)
                {
                    index[i] = i + 1;
                }
            }
            List<double>[] arrs = new List<double>[index.Length + 1];
            for (int i = 0; i < arrs.Length; i++)
            {
                arrs[i] = new List<double>();
            }
            for (int i = 0; i < filestuff.Length; i++)
            {
                string[] spt = filestuff[i].Split(',');
                double x = Convert.ToDouble(spt[0]);
                List<double> curdouble = new List<double>();
                foreach (int f in index)
                {
                    curdouble.Add(Convert.ToDouble(spt[f]));
                }
                arrs[0].Add(x);
                for (int q = 1; q <= curdouble.Count; q++)
                {
                    arrs[q].Add(curdouble[q - 1]);
                }
            }
            double[] basex = arrs[0].ToArray();
            double[] basey = arrs[1].ToArray();
            Plot2D plot = new Plot2D(basex, basey, S);
            for (int i = 2; i < arrs.Length; i++)
            {
                plot.AddCurve(basex, arrs[i].ToArray(), default_colors[i - 2]);
            }
            return plot;
        }
        public void SaveImage(string title, bool using_full_path)
        {
            string path = default_repo + "\\" + title + ".png";
            if (using_full_path) { path = title; }
            Image imgout = (Image)canvas;
            imgout.Save(path);
        }
        public void override_colormap_limits(float new_min, float new_max)
        {
            use_overriden_limits = true;
            override_max = new_max;
            override_min = new_min;
        }
        public void override_colormap_limits(DataSetInfo D)
        {
            use_overriden_limits = true;
            override_max = (float)D.UniversalMax;
            override_min = (float)D.UniversalMin;
        }
        public void reset_colormap_limits()
        {
            use_overriden_limits = false;
        }
        public void SetSketchSettings(Plot2DSettings S)
        {
            settings = S;
            get_positions();
        }
        protected string truncate_string(string input, int number)
        {
            if (sketch_bounds.Xmax-sketch_bounds.Xmin > 0.1 && input.Contains("E-")) { return "0"; }
            if (input.Length <= number) return input;
            else
            {
                string output = string.Empty;
                for (int i = 0; i < number; i++)
                {
                    if (!(i == number - 1 && input[i] == '.'))
                    {
                        output += input[i];
                    }
                }
                return output;
            }
        }
        public void OverrideBounds(RBounds2D newbounds)
        {
            sketch_bounds = newbounds;
            get_positions();

        }
    }
}
