/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace m_435_NCAMR
{
    class Plot2DSettings
    {
        private const int DEFAULTIMAGESIZE = 2000;
        private float width_offset = 0.08f;
        private string fig_title;
        private float fontsize = 40f;
        private const int DEFAULTLABELCOUNT = 5;
        private bool has_axis_values, has_axis_titles, has_cartesian_gridlines;
        private int x_cart_gridlines, y_cart_gridlines;
        private int x_label_count, y_label_count, image_width, image_height;
        private string horizontal_title, vertical_title;
        public float WidthOffset
        {
            get { return width_offset; }
        }
        public void setWidthOffset(float to_set)
        {
            width_offset = to_set;
        }
        public float FontSize
        {
            get { return fontsize; }
        }
        public void setFontSize(float size)
        {
            fontsize = size;
        }
        public string FigureTitle
        {
            get { return fig_title; }
        }
        public void SetFigureTitle(string title)
        {
            fig_title = title;
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
        public bool HasAxisValues
        {
            get { return has_axis_values; }
            set { has_axis_values = value; }
        }
        public bool HasAxisTitles
        {
            get { return has_axis_titles; }
            set { has_axis_titles = value; }
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
        public Plot2DSettings(int imagesize)
        {
            has_axis_values = false;
            has_axis_titles = false;
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
        public Plot2DSettings(bool hasaxisvalues, bool hasaxistitles, bool hascartesians, int xgridcount, int ygridcount, int xlabelcount, int ylabelcount, int imagesize)
        {
            has_axis_values = hasaxisvalues;
            has_axis_titles = hasaxistitles;
            has_cartesian_gridlines = hascartesians;
            x_cart_gridlines = xgridcount;
            y_cart_gridlines = ygridcount;
            x_label_count = xlabelcount;
            y_label_count = ylabelcount;
            image_height = imagesize;
            image_width = imagesize;
            horizontal_title = "x";
            vertical_title = "y";
        }
        public void revert_simple()
        {
            has_axis_titles = false;
            has_axis_values = false;
            has_cartesian_gridlines = false;
        }
        public static Plot2DSettings Fancy()
        {
            return new Plot2DSettings(true, true, true, 10, 10, DEFAULTLABELCOUNT, DEFAULTLABELCOUNT, DEFAULTIMAGESIZE);
        }
    }
}
