/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using _3DSimple;

namespace m_435_NCAMR
{
    class NCAMRNode
    {
        double node_value, x, y;
        public NCAMRNode(double _x, double _y, double _val)
        {
            node_value = _val;
            x = _x;
            y = _y;
        }
        public double Value
        {
            get { return node_value; }
            set { node_value = value; }
        }
        public double X
        {
            get { return x; }
            set { x = value; }
        }
        public double Y
        {
            get { return y; }
            set { y = value; }
        }
        public Point3 toPoint()
        {
            return new Point3(x, y, node_value);
        }
        public static NCAMRNode operator +(NCAMRNode root, Vector3 increment)
        {
            return new NCAMRNode(root.X + increment.I, root.Y + increment.J, root.Value);
        }
        public static NCAMRNode operator -(NCAMRNode root, Vector3 increment)
        {
            return new NCAMRNode(root.X - increment.I, root.Y - increment.J, root.Value);
        }
        public NCAMRNode Clone()
        {
            NCAMRNode b = new NCAMRNode(x, y, node_value);
            return b;
        }
        public override string ToString()
        {
            return "node:" + x.ToString() + ":" + y.ToString() + ":" + node_value.ToString();
        }
        public static NCAMRNode fromstring(string node_string)
        {
            double _x, _y, _val;
            string[] ar = node_string.Split(':');
            if(!(double.TryParse(ar[1], out _x) && double.TryParse(ar[2], out _y) && double.TryParse(ar[3], out _val))) { throw new ArgumentException("String improperly formatted."); }
            else { return new NCAMRNode(_x, _y, _val); }
        }
    }
}
