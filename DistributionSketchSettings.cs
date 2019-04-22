/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;

namespace m_435_NCAMR
{
    class DistributionSketchSettings
    {
        private string fig_title = string.Empty;
        private SketchMode _mode;
        private const int DEFAULTIMAGESIZE = 2000;
        private const int DEFAULTLABELCOUNT = 5;
        private bool has_axis_values, has_info, has_axis_titles, use_special_colormap, has_heatmap, has_cartesian_gridlines;
        private int x_cart_gridlines, y_cart_gridlines;
        private int x_label_count, y_label_count, image_width, image_height;
        private List<Color> special_map;
        private string horizontal_title, vertical_title;
        public string FigureTitle
        {
            get { return fig_title; }
        }
        public bool HasGridlines
        {
            get { return has_cartesian_gridlines; }
            set { has_cartesian_gridlines = value; }
        }
        public int XGridlineCount
        {
            get { return x_cart_gridlines; }
            set { x_cart_gridlines = value; }
        }
        public int YGridlineCount
        {
            get { return y_cart_gridlines; }
            set { y_cart_gridlines = value; }
        }
        public string HorizontalTitle
        {
            get { return horizontal_title; }
            set { horizontal_title = value; }
        }
        public string VerticalTitle
        {
            get { return vertical_title; }
            set { vertical_title = value; }
        }
        public List<Color> ColorMap
        {
            get { return special_map; }
        }
        public int ImageWidth
        {
            get { return image_width; }
            set { image_width = value; }
        }
        public int ImageHeight
        {
            get { return image_height; }
            set { image_height = value; }
        }
        public SketchMode Mode
        {
            get { return _mode; }
            set { _mode = value; }
        }
        public bool HasAxisValues
        {
            get { return has_axis_values; }
            set { has_axis_values = value; }
        }
        public bool HasHeatmap
        {
            get { return has_heatmap; }
            set { has_heatmap = value; }
        }
        public bool HasInfo
        {
            get { return has_info; }
            set { has_info = value; }
        }
        public bool HasAxisTitles
        {
            get { return has_axis_titles; }
            set { has_axis_titles = value; }
        }
        public bool UseCustomColorMap
        {
            get { return use_special_colormap; }
        }
        public int XLabelCount
        {
            get { return x_label_count; }
            set { x_label_count = value; }
        }
        public int YLabelCount
        {
            get { return y_label_count; }
            set { y_label_count = value; }
        }
        public DistributionSketchSettings(SketchMode skmode, int imagesize)
        {
            _mode = skmode;
            has_axis_values = false;
            has_axis_titles = false;
            has_heatmap = false;
            use_special_colormap = false;
            has_info = false;
            has_cartesian_gridlines = false;
            x_cart_gridlines = 10;
            y_cart_gridlines = 10;
            x_label_count = 5;
            y_label_count = 5;
            image_height = imagesize;
            image_width = imagesize;
            horizontal_title = "x";
            vertical_title = "y";
        }
        public DistributionSketchSettings(SketchMode skmode, bool hasaxisvalues, bool hasaxistitles, bool hasheatmap, bool hasinfo, bool hascartesians, int xgridcount, int ygridcount, int xlabelcount, int ylabelcount, int imagesize)
        {
            _mode = skmode;
            has_axis_values = hasaxisvalues;
            has_axis_titles = hasaxistitles;
            has_heatmap = hasheatmap;
            use_special_colormap = false;
            has_info = hasinfo;
            has_cartesian_gridlines = false;
            x_cart_gridlines = xgridcount;
            y_cart_gridlines = ygridcount;
            x_label_count = xlabelcount;
            y_label_count = ylabelcount;
            image_height = imagesize;
            image_width = imagesize;
            horizontal_title = "x";
            vertical_title = "y";
        }
        public void SetColormap(List<Color> colors)
        {
            use_special_colormap = true;
            special_map = colors;
        }
        public void revert_simple()
        {
            has_axis_titles = false;
            has_axis_values = false;
            has_cartesian_gridlines = false;
            has_heatmap = false;
            has_info = false;
        }
        public void SetColormap(DefaultColorMaps t)
        {
            use_special_colormap = true;
            List<Color> the_colors = new List<Color>();
            switch(t)
            {
                case DefaultColorMaps.AMERICAN:
                    {
                        the_colors.Add(Color.Blue);
                        the_colors.Add(Color.LightBlue);
                        the_colors.Add(Color.White);
                        the_colors.Add(Color.Pink);
                        the_colors.Add(Color.Red);
                        break;
                    }
                case DefaultColorMaps.RWBB:
                    {
                        the_colors.Add(Color.Blue);
                        the_colors.Add(Color.DarkBlue);
                        the_colors.Add(Color.Black);
                        the_colors.Add(Color.Gray);
                        the_colors.Add(Color.LightGray);
                        the_colors.Add(Color.White);
                        the_colors.Add(Color.LightPink);
                        the_colors.Add(Color.Pink);
                        the_colors.Add(Color.Red);
                        break;
                    }
                case DefaultColorMaps.GREYSCALE:
                    {
                        the_colors.Add(Color.Black);
                        the_colors.Add(Color.DarkGray);
                        the_colors.Add(Color.Gray);
                        the_colors.Add(Color.LightGray);
                        the_colors.Add(Color.White);
                        break;
                    }
                case DefaultColorMaps.RED:
                    {
                        the_colors.Add(Color.Red);
                        the_colors.Add(Color.OrangeRed);
                        the_colors.Add(Color.Orange);
                        the_colors.Add(Color.Yellow);
                        the_colors.Add(Color.LightYellow);
                        break;
                    }
                case DefaultColorMaps.REV_AMERICAN:
                    {
                        the_colors.Add(Color.Red);
                        the_colors.Add(Color.Pink);
                        the_colors.Add(Color.White);
                        the_colors.Add(Color.LightBlue);
                        the_colors.Add(Color.Blue);
                        break;
                    }
            }
            special_map = the_colors;
        }
        public void RevokeColorMap()
        {
            use_special_colormap = false;
        }
        public enum DefaultColorMaps
        {
            AMERICAN, RWBB, GREYSCALE, RED,  REV_AMERICAN
        }
        public void SetFigureTitle(string figtitle)
        {
            fig_title = figtitle;
        }
        public static DistributionSketchSettings Basic()
        {
            return new DistributionSketchSettings(SketchMode.DISTRIBUTION_ONLY, false, false, false, false, false, DEFAULTLABELCOUNT, DEFAULTLABELCOUNT, 10, 10, DEFAULTIMAGESIZE);
        }
        public static DistributionSketchSettings Fancy()
        {
            return new DistributionSketchSettings(SketchMode.GRID_DISTRIBUTION, true, true, true, true, true, DEFAULTLABELCOUNT, DEFAULTLABELCOUNT, 10, 10, DEFAULTIMAGESIZE);
        }

    }
}
