/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;

namespace m_435_NCAMR
{
    class HybridDistributionSketch2D : DistributionSketch2D
    {
        private NCGrid_Distribution superimpose_grid;
        private bool has_superimposed;
        public HybridDistributionSketch2D(NCGrid_Distribution subject, DistributionSketchSettings S) : base(subject, S)
        {
            has_superimposed = false;
        }
        public void set_superimposed_grid(NCGrid_Distribution _super)
        {
            has_superimposed = true;
            superimpose_grid = _super;
        }
        public void revoke_superimposed_grid()
        {
            has_superimposed = false;
        }
        public override void CreateSketch(bool allow_console_output)
        {
            settings.Mode = SketchMode.DISTRIBUTION_ONLY;
            base.CreateSketch(allow_console_output);
            if (has_superimposed)
            {
                using (Graphics drawer = Graphics.FromImage(canvas))
                {
                    for (int i = 0; i < superimpose_grid.Xcount - 1; i++)
                    {
                        drawer.DrawLine(THICK_GRID_PEN, mapGridNode(superimpose_grid[i, superimpose_grid.Ycount - 1]), mapGridNode(superimpose_grid[i + 1, superimpose_grid.Ycount - 1]));
                        drawer.DrawLine(THICK_GRID_PEN, mapGridNode(superimpose_grid[i, 0]), mapGridNode(superimpose_grid[i + 1, 0]));
                    }
                    for (int j = 0; j < superimpose_grid.Ycount - 1; j++)
                    {
                        drawer.DrawLine(THICK_GRID_PEN, mapGridNode(superimpose_grid[superimpose_grid.Xcount - 1, j]), mapGridNode(superimpose_grid[superimpose_grid.Xcount - 1, j + 1]));
                        drawer.DrawLine(THICK_GRID_PEN, mapGridNode(superimpose_grid[0, j]), mapGridNode(superimpose_grid[0, j + 1]));
                    }
                    for (int i = 0; i < superimpose_grid.Xcount - 1; i++)
                    {
                        for (int j = 0; j < superimpose_grid.Xcount - 1; j++)
                        {
                            drawer.DrawLine(THIN_GRID_PEN, mapGridNode(superimpose_grid[i, j]), mapGridNode(superimpose_grid[i + 1, j]));
                            drawer.DrawLine(THIN_GRID_PEN, mapGridNode(superimpose_grid[i, j]), mapGridNode(superimpose_grid[i, j + 1]));
                        }
                    }
                if (settings.HasHeatmap)
                {
                        draw_heatmap_box(drawer);
                }
            }
            }
        }
        
    }
}
