/* Copyright (c) 2018 William van Noordt */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using _3DSimple;
using System.IO;

namespace m_435_NCAMR
{
    class Matrix
    {
        private const double JUST_ABOUT_ZERO = 1E-10;
        static int zero_count_parameter = 3;
        private int rows, columns;
        private double[] contents;
        private double determinant;
        private bool isSquare;
        private Matrix _inverse = null;
        private Matrix _ref = null;
        private Matrix _rref = null;
        private Matrix _transpose = null;

        public Matrix Inverse
        {
            get
            {
                if (this.isSquare)
                {
                    return _inverse ?? internalInverse(this);
                }
                else
                {
                    throw new ArgumentException("Error: Matrix is not square.");
                }
            }
        }
        public int Dimension
        {
            get { return rows * columns; }
        }
        public bool isSingular
        {
            get { return this.Determinant == 0; }
        }
        public double Determinant
        {
            get
            {
                if (!isSquare) { throw new ArgumentException("Matrix not square."); }
                else
                {
                    if (determinant == double.PositiveInfinity)
                    {
                        determinant = basicDet(this);
                        return determinant;
                    }
                    else
                    {
                        return determinant;
                    }
                }
            }
        }
        public int Rows
        {
            get { return rows; }
        }
        public int Columns
        {
            get { return columns; }
        }
        public bool IsSquare
        {
            get { return isSquare; }
        }
        public enum SystemSolvingScheme
        {
            Basic,
            SortToDiagonal,
            Kaczmarz
        }
        public Matrix(int _rows, int _columns)
        {
            isSquare = _rows == _columns;
            determinant = double.PositiveInfinity;
            contents = new double[_rows * _columns];
            rows = _rows;
            columns = _columns;
            for (int i = 0; i < _rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    contents[j + (i * columns)] = 0;
                }
            }
        }
        public Matrix(int rank)
        {
            isSquare = true;
            determinant = double.PositiveInfinity;
            contents = new double[rank * rank];
            rows = rank;
            columns = rank;
            for (int i = 0; i < rank; i++)
            {
                for (int j = 0; j < rank; j++)
                {
                    contents[j + (i * columns)] = 0;
                }
            }
        }
        public Matrix(int _rows, int _columns, params double[] cts)
        {
            isSquare = _rows == _columns;
            determinant = double.PositiveInfinity;
            contents = new double[_rows * _columns];
            rows = _rows;
            columns = _columns;
            for (int i = 0; i < _rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    contents[j + (i * columns)] = cts[j + (i * columns)];
                }
            }
        }
        public static Matrix eliminateEntry(Matrix A, int row, int col)
        {
            Matrix outMat = new Matrix(A.rows - 1, A.Columns - 1);
            int colOffset = 0;
            for (int i = 0; i < A.Columns - 1; i++)
            {
                int rowOffset = 0;
                if (i == col) { colOffset = 1; }
                for (int j = 0; j < A.Rows - 1; j++)
                {
                    if (j == row) { rowOffset = 1; }
                    outMat[i, j] = A[i + colOffset, j + rowOffset];
                }
            }
            return outMat;
        }
        /// <summary>
        /// "angle" in degrees
        /// </summary>
        public static Matrix RotationMatrix3D(double angle, Vector3 axis)
        {
            Vector3 u = axis.unit();
            double aangle = angle;
            double cos = Math.Cos(aangle);
            double sin = Math.Sin(aangle);
            double onemns = 1 - cos;
            double[] C = new double[9];
            C[0] = cos + (u.I * u.I * onemns);
            C[1] = (u.I * u.J * onemns) - (u.K * sin);
            C[2] = (u.I * u.K * onemns) + (u.J * sin);
            C[3] = (u.J * u.I * onemns) + (u.K * sin);
            C[4] = cos + (u.J * u.J * onemns);
            C[5] = (u.J * u.K * onemns) - (u.I * sin);
            C[6] = (u.K * u.I * onemns) - (u.J * sin);
            C[7] = (u.K * u.J * onemns) + (u.I * sin);
            C[8] = cos + (u.K * u.K * onemns);
            return new Matrix(3, 3, C);
        }
        public double this[int row, int col]
        {
            get { int index = col + (row * columns); return contents[index]; }
            set
            {
                if (row >= rows) { throw new IndexOutOfRangeException("Row index out of range: exceeded maximum value."); }
                if (col >= columns) { throw new IndexOutOfRangeException("Column index out of range: exceeded maximum value."); }
                if (row < 0) { throw new IndexOutOfRangeException("Row index out of range: exceeded minimum value."); }
                if (col < 0) { throw new IndexOutOfRangeException("Column index out of range: exceeded minimum value."); }
                double inter = value;
                int index = col + (row * columns);
                contents[index] = inter;
                clearDependent();
            }
        }
        public static Matrix ConcatenateRight(Matrix A, Matrix B)
        {
            if (A.rows != B.rows)
            {
                throw new Exception("Error: Inconsistent dimensions.");
            }
            int rs = A.rows;
            int firstcols = A.columns;
            int secondcols = B.columns;
            int totalcols = firstcols + secondcols;
            Matrix output = new Matrix(rs, totalcols);
            for (int i = 0; i < rs; i++)
            {
                for (int j = 0; j < firstcols; j++)
                {
                    output[i, j] = A[i, j];
                }
                for (int j = firstcols; j < totalcols; j++)
                {
                    output[i, j] = B[i, j-firstcols];
                }
            }
            return output;
        }
        public double this[int index]
        {
            get
            {
                if (rows == 1 || columns == 1)
                {
                    return contents[index];
                }
                else
                {
                    throw new Exception("Error: Matrix is not a row or column vector.");
                }
            }
            set
            {
                if (rows == 1 || columns == 1)
                {
                    contents[index] = value;
                    clearDependent();
                }
                else
                {
                    throw new Exception("Error: Matrix is not a row or column vector.");
                }
            }
        }
        private double basicDet(Matrix A)
        {
            if (A.Rows == 1 && A.Columns == 1)
            {
                return A[0, 0];
            }
            else
            {
                double acc = 0;
                double sgn = 1;
                for (int i = 0; i < A.columns; i++)
                {
                    acc += A[i, 0] * sgn * basicDet(Matrix.eliminateEntry(A, 0, i));
                    sgn *= -1;
                }
                return acc;
            }
        }
        public static Matrix Identity(int m)
        {
            Matrix I = new Matrix(m);
            for (int i = 0; i < m; i++)
            {
                I[i, i] = 1;
            }
            return I;
        }
        public Matrix Clone()
        {
            Matrix X = new Matrix(rows, columns);
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    X[i, j] = this[i, j];
                }
            }
            return X;
        }
        private Matrix internalInverse(Matrix A)
        {
            //will need modification to switch around rows if first diagonal entry is zero.
            int m = A.rows;
            Matrix invs = Matrix.Identity(m);
            Matrix clone = this.Clone();
            for (int diagrow = 0; diagrow < m; diagrow++)
            {
                double scl1 = 1 / clone[diagrow, diagrow];
                invs.RowScale(diagrow, scl1);
                clone.RowScale(diagrow, scl1);
                for (int i = 0; i < m; i++)
                {
                    if (i != diagrow && clone[i, diagrow] != 0)
                    {
                        double scl2 = -1 * clone[i, diagrow];
                        invs.RowAdd(i, diagrow, scl2);
                        clone.RowAdd(i, diagrow, scl2);
                    }
                }
            }
            _inverse = invs;
            return invs;
        }
        private Matrix internalInverse(Matrix A, bool doLatexOutput, string path, bool approximmateZero)
        {
            List<string> contents = new List<string>();
            int m = A.rows;
            contents.Add("Let $M = " + this.Latexoutput(approximmateZero) + "$. Then, the following operations are used to invert $M$ to find $M^{-1}$:");
            Matrix clone = this;
            Matrix invs = Matrix.Identity(m);
            contents.Add("\\begin{equation*}");
            contents.Add(clone.Latexoutput(approximmateZero));
            contents.Add("\\longleftrightarrow");
            contents.Add(invs.Latexoutput(approximmateZero));
            contents.Add("\\end{equation*}");
            for (int diagrow = 0; diagrow < m; diagrow++)
            {
                double scl1 = 1 / clone[diagrow, diagrow];
                invs.RowScale(diagrow, scl1);
                clone.RowScale(diagrow, scl1);
                contents.Add("\\begin{equation*}");
                contents.Add(clone.Latexoutput(approximmateZero));
                contents.Add("\\longleftrightarrow");
                contents.Add(invs.Latexoutput(approximmateZero));
                contents.Add("\\end{equation*}");
                for (int i = 0; i < m; i++)
                {
                    if (i != diagrow)
                    {
                        double scl2 = -1 * clone[i, diagrow];
                        invs.RowAdd(i, diagrow, scl2);
                        clone.RowAdd(i, diagrow, scl2);
                    }
                }
                contents.Add("\\begin{equation*}");
                contents.Add(clone.Latexoutput(approximmateZero));
                contents.Add("\\longleftrightarrow");
                contents.Add(invs.Latexoutput(approximmateZero));
                contents.Add("\\end{equation*}");
            }
            _inverse = invs;
            contents.Add("Therefore, $M^{-1} = " + _inverse.Latexoutput(approximmateZero) + "$.");
            File.WriteAllLines(path, contents.ToArray());
            Console.WriteLine("Latex file output.");
            return invs;
        }
        public void RowScale(int row_zer_bas, double scale)
        {
            for (int j = 0; j < columns; j++)
            {
                this[row_zer_bas, j] *= scale;
            }
        }
        private static bool isZero(double D)
        {
            return Math.Abs(D) < JUST_ABOUT_ZERO;
        }
        public void RowAdd(int row_add_to, int row_add_from, double prescale)
        {
            for (int j = 0; j < columns; j++)
            {
                this[row_add_to, j] += (prescale * this[row_add_from, j]);
            }
        }
        public static Matrix operator +(Matrix A, Matrix B)
        {
            if (A.Rows != B.Rows || A.Columns != B.Columns)
            {
                throw new ArgumentException("Error: matrix dimensions must agree.");
            }
            else
            {
                Matrix C = new Matrix(A.rows, B.Columns);
                for (int i = 0; i < A.Rows; i++)
                {
                    for (int j = 0; j < B.columns; j++)
                    {
                        C[i, j] = A[i, j] + B[i, j];
                    }
                }
                return C;
            }
        }
        private static double vec_inner_prod(Matrix A, Matrix B)
        {
            double acc = 0;
            for (int i = 0; i < A.Dimension; i++)
            {
                acc += A[i] * B[i];
            }
            return acc;
        }
        public static Matrix operator -(Matrix A, Matrix B)
        {
            if (A.Rows != B.Rows || A.Columns != B.Columns)
            {
                throw new ArgumentException("Error: matrix dimensions must agree.");
            }
            else
            {
                Matrix C = new Matrix(A.rows, B.Columns);
                for (int i = 0; i < A.Rows; i++)
                {
                    for (int j = 0; j < B.columns; j++)
                    {
                        C[i, j] = A[i, j] - B[i, j];
                    }
                }
                return C;
            }
        }
        public static Matrix operator *(Matrix A, Matrix B)
        {
            if (A.Columns != B.Rows)
            {
                throw new ArgumentException("Error: Matrix dimensions are incompatible with multiplication");
            }
            else
            {
                Matrix C = new Matrix(A.Rows, B.Columns);
                for (int i = 0; i < C.Rows; i++)
                {
                    for (int j = 0; j < C.Columns; j++)
                    {
                        double acc = 0;
                        for (int y = 0; y < B.Rows; y++)
                        {
                            acc += A[i, y] * B[y, j];
                        }
                        C[i, j] = acc;
                    }
                }
                return C;
            }
        }
        private void clearDependent()
        {
            determinant = double.PositiveInfinity;
            _inverse = null;
            _ref = null;
            _rref = null;
            _transpose = null;
        }
        /// <summary>
        /// Solves an N by N system of equations, with right-hand N-vecctor.
        /// </summary>
        /// <param name="coeffs">Must be a square matrix</param>
        /// <param name="rhs">Right hand side</param>
        /// <returns></returns>
        public static Matrix solve_system(Matrix coeffs, Matrix rhs)
        {
            return solve_system(coeffs, rhs, SystemSolvingScheme.Basic, true, false);
        }
        /// <summary>
        /// Solves an N by N system of equations, with right-hand N-vecctor.
        /// </summary>
        /// <param name="coeffs">Must be a square matrix</param>
        /// <param name="rhs">Right hand side</param>
        /// <param name="scheme">Indicates the scheme used to solve the system</param>
        /// <returns></returns>
        public static Matrix solve_system(Matrix coeffs, Matrix rhs, SystemSolvingScheme scheme)
        {
            return solve_system(coeffs, rhs, scheme, true, false);
        }
        /// <summary>
        /// Solves an N by N system of equations, with right-hand N-vecctor.
        /// </summary>
        /// <param name="coeffs">Must be a square matrix</param>
        /// <param name="rhs">Right hand side</param>
        /// <param name="suppress_determinant_validation">Used to supress determinant validation, lower computing time, but risks error</param>
        /// <returns></returns>
        public static Matrix solve_system(Matrix coeffs, Matrix rhs, bool suppress_determinant_validation)
        {
            return solve_system(coeffs, rhs, SystemSolvingScheme.Basic, suppress_determinant_validation, false);
        }
        /// <summary>
        /// Solves an N by N system of equations, with right-hand N-vecctor.
        /// </summary>
        /// <param name="coeffs">Must be a square matrix</param>
        /// <param name="rhs">Right hand side</param>
        /// <param name="scheme">Indicates the scheme used to solve the system</param>
        /// <param name="suppress_determinant_validation">Used to supress determinant validation, lower computing time, but risks error</param>
        /// <returns></returns>
        public static Matrix solve_system(Matrix coeffs, Matrix rhs, SystemSolvingScheme scheme, bool suppress_determinant_validation, bool enable_console_output)
        {
            if (!coeffs.IsSquare || (rhs.Columns != 1) || (rhs.Rows != coeffs.Rows))
            {
                throw new ArgumentException("System has incompatible dimensions.");
            }
            if (!suppress_determinant_validation)
            {
                if (coeffs.isSingular)
                {
                    throw new ArgumentException("Coefficient matrix is singular.");
                }
            }
            switch (scheme)
            {
                case SystemSolvingScheme.Basic:
                    {
                        return coeffs.Inverse * rhs;
                    }
                case SystemSolvingScheme.SortToDiagonal:
                    {
                        return solve_using_sort_diagonal(coeffs, rhs, true);
                    }
                case SystemSolvingScheme.Kaczmarz:
                    {
                        return kaczmarz_solve(coeffs, rhs, 100, 1E-5);
                    }
                default:
                    {
                        throw new ArgumentException("Invalid scheme");
                    }
            }

        }
        private static Matrix solve_using_sort_diagonal(Matrix coeffs, Matrix rhs, bool enable_console_output)
        {
            int rws = coeffs.rows;
            int cols = coeffs.columns + 1;
            Matrix augment = new Matrix(rws, cols + 1);
            if (enable_console_output) { Console.WriteLine("Augmenting..."); }
            for (int i = 0; i < rws; i++)
            {
                if (enable_console_output && i%(rws/13) == 0)
                {
                    Console.WriteLine(((100 * i) / rws).ToString() + "%");
                }
                for (int j = 0; j < cols - 1; j++)
                {
                    augment[i, j] = coeffs[i, j];
                }
                augment[i, cols - 1] = rhs[i];
                augment[i, cols] = 0;
            }
            for (int cur_row = 0; cur_row < augment.rows; cur_row++)
            {
                bool zero = true;
                for (int cur_col = 0; cur_col < augment.columns && zero; cur_col++)
                {
                    zero = isZero(augment[cur_row, cur_col]);
                    if (zero)
                    {
                        augment[cur_row, cols] = cur_col + 1;
                    }
                }
            }
            sort_by_last(augment);
            int m = augment.rows;
            int n = augment.columns - 1;
            if (enable_console_output) { Console.WriteLine("Solving..."); }
            for (int diagrow = 0; diagrow < m; diagrow++)
            {
                if (enable_console_output && diagrow%(m/13) == 0)
                {
                    Console.WriteLine(((100*diagrow)/m).ToString()+"%");
                }
                double scl1 = 1 / augment[diagrow, diagrow];
                augment.RowScale(diagrow, scl1);
                int zerocount = 0;
                for (int i = 0; i < m && zerocount < zero_count_parameter; i++)
                {
                    if (Math.Abs(augment[i, diagrow]) < 1E-10) { zerocount++; }
                    if (i != diagrow)
                    {
                        double scl2 = -1 * augment[i, diagrow];
                        augment.RowAdd(i, diagrow, scl2);
                    }
                }
            }
            double[] result = new double[augment.rows];
            for (int i = 0; i < augment.rows; i++)
            {
                result[i] = augment[i, augment.columns - 2];
            }
            return new Matrix(augment.rows, 1, result);
        }
        public void RowSwap(int row, int other_row)
        {
            for (int i = 0; i < this.columns; i++)
            {
                double med = this[other_row, i];
                this[other_row, i] = this[row, i];
                this[row, i] = med;
            }
        }
        public double sumrow(int row)
        {
            double sum = 0;
            for (int i = 0; i < columns; i++)
            {
                sum += this[row, i];
            }
            return sum;
        }
        public void exportToFile(string path)
        {
            string filetype = path.Split('.')[1];
            switch (filetype)
            {
                case "csv":
                    {
                        List<string> filestuff = new List<string>();
                        for (int i = 0; i < rows; i++)
                        {
                            string currentline = string.Empty;
                            currentline += this[i, 0].ToString();
                            for (int j = 1; j < columns; j++)
                            {
                                currentline += "," + this[i, j].ToString();
                            }
                            filestuff.Add(currentline);
                        }
                        File.WriteAllLines(path, filestuff.ToArray());
                        break;
                    }
                case "txt":
                    {
                        File.WriteAllLines(path, this.ToStrings());
                        break;
                    }
                default:
                    {
                        throw new Exception("Error: Invalid filetype.");
                    }
            }
        }
        public string[] ToStrings()
        {
            List<string> sts = new List<string>();
            int shelfsize = 8;
            int trunc = 6;
            for (int i = 0; i < this.Rows; i++)
            {
                string line = string.Empty;
                for (int j = 0; j < this.Columns; j++)
                {
                    string entry = this[i, j].ToString();
                    string tr_entry = entry.Substring(0, Math.Min(trunc, entry.Length));
                    line += tr_entry.PadRight(shelfsize);
                }
                sts.Add(line);
            }
            return sts.ToArray();
        }
        public string[] ToStrings(int _trunc, int _shelfsize)
        {
            List<string> sts = new List<string>();
            int shelfsize = _shelfsize;
            int trunc = _trunc;
            for (int i = 0; i < this.Rows; i++)
            {
                string line = string.Empty;
                for (int j = 0; j < this.Columns; j++)
                {
                    string entry = this[i, j].ToString();
                    string tr_entry = entry.Substring(0, Math.Min(trunc, entry.Length));
                    line += tr_entry.PadRight(shelfsize);
                }
                sts.Add(line);
            }
            return sts.ToArray();
        }
        private static void sort_by_last(Matrix M)
        {
            int last = M.rows - 1;
            quicksort(M, 0, last);
        }
        private static void quicksort(Matrix M, int lo, int hi)
        {
            if (lo < hi)
            {
                int pi = partition(M, lo, hi);

                quicksort(M, lo, pi - 1);
                quicksort(M, pi + 1, hi);
            }
        }
        private static int partition(Matrix M, int lo, int hi)
        {
            int lst = M.columns - 1;
            int pivot = (int)M[hi, lst];
            int i = (lo - 1);
            for (int j = lo; j <= hi - 1; j++)
            {
                if (M[j, lst] <= pivot)
                {
                    i++;
                    M.RowSwap(i, j);
                }
            }
            M.RowSwap(i + 1, hi);
            return (i + 1);
        }
        public string Latexoutput(bool approximate_zeros)
        {
            string output = "\\begin{bmatrix}";
            for (int j = 0; j < rows; j++)
            {
                output += ((float)(this[j, 0])).ToString();
                for (int i = 1; i < columns; i++)
                {
                    if (approximate_zeros && Math.Abs(this[j, i]) < 1E-9)
                    {
                        output += "&" + "0";
                    }
                    else
                    {
                        output += "&" + ((float)(this[j, i])).ToString();
                    }
                }
                output += "\\" + "\\";
            }
            return output + "\\end{bmatrix}";
        }
        public static Matrix ext_kaczmarz_radius(Matrix A, Matrix b, double tolerance, int maxiterations, int radius)
        {
            Matrix xk = Matrix.constvec(b.Rows, 1);
            int n = b.Rows;
            int k = 0;
            Matrix xprev = xk.Clone();
            while (k <= maxiterations)
            {
                double complete_residual = Matrix.VecNorm(A * xk - b);
                for (int i = 0; i < xk.Dimension; i++)
                {
                    Matrix ai = A.RowVector(i);
                    double prod = vec_prod_radius(ai, xk, radius, i);
                    double top = prod - b[i];
                    double mag = Matrix.VecNorm(ai);
                    double bottom = mag * mag;
                    double scale = top / bottom;
                    Matrix subtract = scale * ai;
                    xprev = xk.Clone();
                    xk = xk - subtract.Transpose();
                }
                double residual = Matrix.VecNorm(xk - xprev) / Matrix.VecNorm(xk);
                Console.WriteLine(complete_residual.ToString() + "," + residual.ToString());
                if (residual < tolerance)
                {
                    break;
                }
                k++;
            }
            return xk;
        }
        private static double vec_prod_radius(Matrix a, Matrix b, int assumedradius, int center)
        {
            int bottom = Math.Max(0, center - assumedradius);
            int top = Math.Min(b.Dimension - 1, center + assumedradius);
            double acc = 0;
            for (int i = bottom; i < top; i++)
            {
                acc += a[i] * b[i];
            }
            return acc;
        }
        public static Matrix FromCsv(string path)
        {
            if (!path.EndsWith(".csv")) { throw new Exception("Invalid file type."); }
            List<double[]> cts = new List<double[]>();
            string[] stuff = File.ReadAllLines(path);
            int rows = stuff.Length;
            int columns = stuff[0].Split(',').Length;
            Matrix M = new Matrix(rows, columns);
            for (int i = 0; i < rows; i++)
            {
                string[] line = stuff[i].Split(',');
                for (int j = 0; j < columns; j++)
                {
                    double entry;
                    if (!double.TryParse(line[j], out entry))
                    {
                        throw new Exception("Error: invalid entry at " + i.ToString() + ", " + j.ToString());
                    }
                    else
                    {
                        M[i, j] = entry;
                    }
                }
            }
            return M;
        }
        public Matrix WriteInverseToFile(string path)
        {
            return internalInverse(this, true, path, true);
        }
        public Matrix Transpose()
        {
            Matrix output = new Matrix(this.columns, this.rows);
            for (int i = 0; i < output.rows; i++)
            {
                for (int j = 0; j < output.Columns; j++)
                {
                    output[i, j] = this[j, i];
                }
            }
            return output;
        }
        private static Matrix kaczmarz_solve(Matrix A, Matrix b, int maxiterations, double tolerance)
        {
            //initial guess
            List<double> xs = new List<double>();
            List<double> ys = new List<double>();
            List<double> zs = new List<double>();
            Matrix xk = Matrix.constvec(b.Rows, 1);
            int n = b.Rows;
            int k = 0;
            int ct = 0;
            Matrix xprev = xk.Clone();
            while (k <= maxiterations)
            {
                double complete_residual = Matrix.VecNorm(A * xk - b);
                double residual = Matrix.VecNorm(xk - xprev) / Matrix.VecNorm(xk);
                for (int i = 0; i < xk.Dimension; i++)
                {
                    xs.Add(k);
                    ct++;
                    Matrix ai = A.RowVector(i);
                    double prod = Matrix.vec_inner_prod(ai, xk);
                    double top = prod - b[i];
                    double mag = Matrix.VecNorm(ai);
                    double bottom = mag * mag;
                    double scale = top / bottom;
                    Matrix subtract = scale * ai;
                    xprev = xk.Clone();
                    xk = xk - subtract.Transpose();
                    ys.Add(complete_residual);
                    zs.Add(5*residual);
                    if (residual < tolerance && k > 2)
                    {
                        return xk;
                    }
                }
                Console.WriteLine(complete_residual.ToString() + "," + residual.ToString());
                k++;
            }
            return xk;
            
        }
        public static Matrix operator*(double scale, Matrix A)
        {
            Matrix output = new Matrix(A.rows, A.columns);
            for (int i = 0; i < A.rows; i++)
            {
                for (int j = 0; j < A.columns; j++)
                {
                    output[i, j] = scale * A[i, j];
                }
            }
            return output;
        }
        public static double VecNorm(Matrix vectormatrix)
        {
            double acc = 0;
            for (int i = 0; i < vectormatrix.Dimension; i++)
            {
                acc += vectormatrix[i] * vectormatrix[i];
            }
            return Math.Sqrt(acc);
        }
        public Matrix ColumnVector(int i)
        {
            Matrix output = new Matrix(rows, 1);
            for (int z = 0; z < rows; z++)
            {
                output[z] = this[z, i];
            }
            return output;
        }
        public Matrix RowVector(int i)
        {
            Matrix output = new Matrix(1, columns);
            for (int z = 0; z < columns; z++)
            {
                output[z] = this[i, z];
            }
            return output;
        }
        private static Matrix constvec(int rowdimension, double val)
        {
            Matrix z = new Matrix(rowdimension, 1);
            for (int i = 0; i < rowdimension; i++)
            {
                z[i] = val;
            }
            return z;
        }
    }
}
