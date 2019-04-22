/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;
using ColorMapping;

namespace m_435_NCAMR
{
    enum SketchMode
    {
        EMPTY,
        GRID_DISTRIBUTION,
        GRID_ONLY,
        DISTRIBUTION_ONLY
    }
    class DistributionSketch2D
    {
        protected const string default_repo = @"C:\Users\Will\Desktop\Folders\MATH435\repo\images";
        protected DistributionSketchSettings settings;
        protected float override_max;
        protected float override_min;
        protected bool use_overriden_limits = false;

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
        protected const float GRIDLINE_THICKNESS = 2.3f;
        protected const float X_AXIS_TITLE_X = 0.5f;
        protected const float Y_AXIS_TITLE_Y = 0.5f;
        protected float Y_AXIS_TITLE_X = 0.31f * SQUARE_BORDER_THICK;
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

        protected float HEATMAP_LABELS_X = 0.87f;
        protected static float HEATMAP_LABELS_YU = 0.33f;
        protected static float HEATMAP_LABELS_YL = 0.64f;

        protected float X_LABEL_BASE_Y = 0.94f * SQUARE_BORDER_THICK;
        protected float X_LABEL_BASE_X = 1.0f * SQUARE_BORDER_THICK;
        protected float Y_LABEL_BASE_Y = 1.0f * SQUARE_BORDER_THICK;
        protected float Y_LABEL_BASE_X = 0.54f * SQUARE_BORDER_THICK;
        protected float VALUE_TEXT_PROPORTION = 40f / 2000f;
        protected float TITLE_TEXT_PROPORTION = 65f / 2000f;
        protected const float THIN_PEN = 1.5f;
        protected const float THICK_PEN = 2.4f;
        protected Pen THICK_GRID_PEN = new Pen(Color.Black, THICK_PEN);
        protected Pen THIN_GRID_PEN = new Pen(Color.Black, THIN_PEN);
        protected Color BACKGROUND_COLOR = Color.White;
        protected Color TEXTCOLOR = Color.Black;

        //calculated relatives
        protected float grid_xmin;
        protected float grid_xmax;
        protected float grid_ymin;
        protected float grid_ymax;
        protected float pix_per_unit_x;
        protected float pix_per_unit_y;

        protected Bitmap canvas;
        protected int x_pixel_count, y_pixel_count;
        protected NCGrid_Distribution subject;
        public SketchMode Sketch_Mode
        {
            get { return settings.Mode; }
            set { settings.Mode = value; }
        }
        public Bitmap Canvas
        {
            get { return canvas; }
        }
        public DistributionSketch2D(NCGrid_Distribution _subject, DistributionSketchSettings S)
        {
            settings = S;
            subject = _subject;
            canvas = new Bitmap(S.ImageWidth, S.ImageHeight);
            x_pixel_count = canvas.Width;
            y_pixel_count = canvas.Height;
            get_positions();

        }
        public virtual void CreateSketch(bool allow_console_output)
        {
            using ( Graphics drawer = Graphics.FromImage(canvas))
            {
                drawer.Clear(BACKGROUND_COLOR);
                bool no_dt = false;
                DateTime then = DateTime.Now;
                switch (settings.Mode)
                {
                    case SketchMode.EMPTY: { if (allow_console_output) { Console.WriteLine("Warning: Sketchmode set to \"EMPTY\". No drawing operations executed."); no_dt = true; } break; }
                    case SketchMode.GRID_ONLY:
                        {
                            draw_grid(drawer);
                            break;
                        }
                    case SketchMode.DISTRIBUTION_ONLY:
                        {
                            draw_distrib(drawer);
                            break;
                        }
                    case SketchMode.GRID_DISTRIBUTION:
                        {
                            draw_distrib(drawer);
                            draw_grid(drawer);
                            break;
                        }
                }
                if (settings.HasAxisValues) { draw_axis_values(drawer); }
                if (settings.HasAxisTitles) { draw_axis_titles(drawer); }
                if (settings.HasGridlines) { draw_gridlines(drawer); }
                if (settings.HasHeatmap) { draw_heatmap_box(drawer); }
                DateTime now = DateTime.Now;
                draw_fig_title(drawer);
                if (allow_console_output && !no_dt) { Console.WriteLine("Sketch generated in " + (now - then).Milliseconds.ToString() + " milliseconds."); }
            }
        }
        protected void get_positions()
        {
            if (settings.HasAxisValues || settings.HasAxisTitles)
            {
                MAIN_SUBJECT_LOWER_LEFT_X_2D = SQUARE_BORDER_THICK;
                MAIN_SUBJECT_LOWER_LEFT_Y_2D = SQUARE_BORDER_THICK;
                MAIN_SUBJECT_UPPER_RIGHT_X_2D = 1 - SQUARE_BORDER_THICK;
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
            pix_per_unit_x = (grid_xmax - grid_xmin) / (float)(subject.Xmax - subject.Xmin);
            pix_per_unit_y = (grid_ymax - grid_ymin) / (float)(subject.Ymax - subject.Ymin);
        }
        protected void draw_axis_values(Graphics drawer)
        {
            Brush textbrush = new SolidBrush(TEXTCOLOR);
            float xlabel_abs_pos_x = ((float)x_pixel_count * X_LABEL_BASE_X);
            float xlabel_abs_increment = (((1 - 2 * (Y_LABEL_BASE_Y)) * (float)x_pixel_count) / (float)(settings.XLabelCount-1));
            float ylabel_abs_pos_y = ((float)y_pixel_count * (Y_LABEL_BASE_X));
            float ylabel_abs_increment = (((1 - 2 * Y_LABEL_BASE_Y) * (float)y_pixel_count) / (float)(settings.YLabelCount-1));
            float xlabel_abs_pos_y = (X_LABEL_BASE_Y * (float)x_pixel_count);
            float ylabel_abs_pos_x = ((Y_LABEL_BASE_X-LABEL_OFFSET_Y) * (float)y_pixel_count);
            double delta_x = (subject.Xmax - subject.Xmin) / (settings.XLabelCount - 1);
            double delta_y = (subject.Ymax - subject.Ymin) / (settings.YLabelCount - 1);
            for (int i = 0; i < settings.XLabelCount; i++)
            {
                drawer.DrawString(truncate_string((subject.Xmin + i * delta_x).ToString(),4), new Font("arial", VALUE_TEXT_PROPORTION * settings.ImageHeight), textbrush, new PointF(xlabel_abs_pos_x + i * xlabel_abs_increment, settings.ImageHeight - xlabel_abs_pos_y));
            }
            for (int j = 0;  j< settings.YLabelCount; j++)
            {
                drawer.DrawString(truncate_string((subject.Ymax - j * delta_y).ToString(),4), new Font("arial", VALUE_TEXT_PROPORTION * settings.ImageHeight), textbrush, new PointF(ylabel_abs_pos_y, ylabel_abs_pos_x + j*ylabel_abs_increment));
            }
        }
        protected void draw_heatmap_box(Graphics drawer)
        {
            float abspos_hb_ll_x = HEATBOX_LL_X * settings.ImageWidth;
            float abspos_hb_ll_y = HEATBOX_LL_Y * settings.ImageHeight;
            float abspos_hb_ur_x = HEATBOX_UR_X * settings.ImageWidth;
            float abspos_hb_ur_y = HEATBOX_UR_Y * settings.ImageHeight;
            float abspos_border_upper = HEATBOX_BORDER_TOP * settings.ImageHeight;
            float abspos_border_lower = HEATBOX_BORDER_BOTTOM * settings.ImageHeight;
            float abspos_border_left = HEATBOX_BORDER_LEFT * settings.ImageWidth;
            float abspos_border_right = HEATBOX_BORDER_RIGHT * settings.ImageWidth;
            PointF[] border =
            {
                new PointF(abspos_border_left, abspos_border_upper),
                new PointF(abspos_border_left, abspos_border_lower),
                new PointF(abspos_border_right, abspos_border_lower),
                new PointF(abspos_border_right, abspos_border_upper)
            };
            Brush borderbrush = new SolidBrush(Color.White);
            drawer.FillPolygon(borderbrush, border);
            int heatmap_pixel_count = (int)Math.Abs((abspos_hb_ll_y - abspos_hb_ur_y));
            for (int i = 0; i <= heatmap_pixel_count; i++)
            {
                Pen drawingpen;
                if (settings.UseCustomColorMap) { drawingpen = new Pen(ColorMap.MapColorMultiple(settings.ColorMap, (double)i, 0, (double)heatmap_pixel_count),1.0f); }
                else { drawingpen = new Pen(ColorMap.MapRainbowColor((float)i, (float)heatmap_pixel_count, 0),1.0f); }
                PointF left = new PointF(abspos_hb_ll_x, abspos_hb_ur_y - i);
                PointF right = new PointF(abspos_hb_ur_x, abspos_hb_ur_y - i);
                drawer.DrawLine(drawingpen, left, right);
            }
            float abs_upper_label_y = HEATMAP_LABELS_YU * settings.ImageHeight;
            float abs_lower_label_y = HEATMAP_LABELS_YL * settings.ImageHeight;
            float abs_label_x = HEATMAP_LABELS_X * settings.ImageWidth;
            string upper_label, lower_label;
            if (use_overriden_limits)
            {
                upper_label = truncate_string(override_max.ToString(),4);
                lower_label = truncate_string(override_min.ToString(), 4);
            }
            else
            {
                upper_label = truncate_string(subject.MaxValue.ToString(), 4);
                lower_label = truncate_string(subject.MinValue.ToString(), 4);
            }
            Brush textbrush = new SolidBrush(TEXTCOLOR);
            drawer.DrawString(upper_label, new Font("arial", VALUE_TEXT_PROPORTION * settings.ImageHeight), textbrush, new PointF(abs_label_x, abs_upper_label_y));
            drawer.DrawString(lower_label, new Font("arial", VALUE_TEXT_PROPORTION * settings.ImageHeight), textbrush, new PointF(abs_label_x, abs_lower_label_y));
        }
        protected void draw_gridlines(Graphics drawer)
        {
            Pen gridlinepen = new Pen(Color.FromArgb(130, 0, 0, 0), GRIDLINE_THICKNESS);
            NCAMRNode[,] gridpoints = new NCAMRNode[settings.XGridlineCount, settings.YGridlineCount];
            double coordinate_dx = (subject.Xmax - subject.Xmin) / (settings.XGridlineCount - 1);
            double coordinate_dy = (subject.Ymax - subject.Ymin) / (settings.YGridlineCount - 1);
            for (int i = 0; i < settings.XGridlineCount; i++)
            {
                for (int j = 0; j < settings.YGridlineCount; j++)
                {
                    gridpoints[i, j] = new NCAMRNode(subject.Xmin + i * coordinate_dx, subject.Ymin + j * coordinate_dy, 0);
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
                    drawer.DrawLine(gridlinepen, mapGridNode(gridpoints[j+1, i]), mapGridNode(gridpoints[j, i]));
                }
            }
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
        protected void draw_grid(Graphics drawer)
        {
            
            for (int i = 0; i < subject.Xcount - 1; i++)
            {
                drawer.DrawLine(THICK_GRID_PEN, mapGridNode(subject[i, subject.Ycount - 1]), mapGridNode(subject[i + 1, subject.Ycount - 1]));
                drawer.DrawLine(THICK_GRID_PEN, mapGridNode(subject[i, 0]), mapGridNode(subject[i + 1, 0]));
            }
            for (int j = 0; j < subject.Ycount - 1; j++)
            {
                drawer.DrawLine(THICK_GRID_PEN, mapGridNode(subject[subject.Xcount - 1, j]), mapGridNode(subject[subject.Xcount - 1, j + 1]));
                drawer.DrawLine(THICK_GRID_PEN, mapGridNode(subject[0, j]), mapGridNode(subject[0, j + 1]));
            }
            for (int i = 0; i < subject.Xcount - 1; i++)
            {
                for (int j = 0; j < subject.Xcount - 1; j++)
                {
                    drawer.DrawLine(THIN_GRID_PEN, mapGridNode(subject[i, j]), mapGridNode(subject[i + 1, j]));
                    drawer.DrawLine(THIN_GRID_PEN, mapGridNode(subject[i, j]), mapGridNode(subject[i, j+1]));
                }
            }
        }
        public void update_label_y_pos(double d)
        {
            Y_AXIS_TITLE_X = (float)d;
        }
        protected void draw_distrib(Graphics drawer)
        {
            double min = subject.MinValue;
            double max = subject.MaxValue;
            if (use_overriden_limits)
            {
                max = override_max;
                min = override_min;
            }
            for (int i = 0; i < subject.Xcount - 1; i++)
            {
                for (int j = 0; j < subject.Ycount - 1; j++)
                {
                    NCAMRNode[] nodes = { subject[i, j], subject[i+1, j], subject[i+1, j+1], subject[i, j+1] };
                    double zavg = 0.25d * (nodes[0].Value+ nodes[1].Value+ nodes[2].Value+ nodes[3].Value);
                    if (Math.Abs(max - min) < 0.001d)
                    {
                        max = 1;
                        min = 0;
                        zavg = 0;
                    }
                    SolidBrush B;
                    if (settings.UseCustomColorMap)
                    {
                        B = new SolidBrush(ColorMap.MapColorMultiple(settings.ColorMap, zavg, min, max));
                    }
                    else
                    {
                        B = new SolidBrush(ColorMap.MapRainbowColor((float)zavg, (float)max, (float)min));
                    }
                    PointF[] pts = new PointF[4];
                    for (int k = 0; k < 4; k++)
                    {
                        pts[k] = mapGridNode(nodes[k]);
                    }
                    drawer.FillPolygon(B, pts);
                }
            }
        }
        protected PointF mapGridNode(NCAMRNode toMap)
        {
            return new PointF(grid_xmin + pix_per_unit_x * (float)(toMap.X - subject.Xmin), grid_ymin + pix_per_unit_y * (float)((subject.Ymax - toMap.Y)));
        }
        public void SaveImage(string title, bool using_full_path)
        {
            string path = default_repo + "\\" + title + ".bmp";
            if (using_full_path) { path = title; }
            canvas.Save(path);
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
        public void SetSketchSettings(DistributionSketchSettings S)
        {
            settings = S;
            get_positions();
        }
        protected void draw_fig_title(Graphics drawer)
        {
            float absx = settings.ImageWidth * FIGURETITLE_X;
            float absy = settings.ImageHeight * FIGURETITLE_Y;
            Brush textbrush = new SolidBrush(TEXTCOLOR);
            drawer.DrawString(settings.FigureTitle, new Font("arial", 1.9f*VALUE_TEXT_PROPORTION * settings.ImageHeight), textbrush, new PointF(absx, absy));
        }
        protected string truncate_string(string input, int number)
        {
            string supplementary = string.Empty;
            if (input.Contains("E"))
            {
                int idx = input.LastIndexOf('E');
                for (int i = idx; i < input.Length; i++)
                {
                    supplementary += input[i];
                }
            }
            if (input.Length <= number) return input;
            else
            {
                string output = string.Empty;
                for (int i = 0; i < number;  i++)
                {
                    output += input[i];
                }
                return output+supplementary;
            }
        }
    }
}
