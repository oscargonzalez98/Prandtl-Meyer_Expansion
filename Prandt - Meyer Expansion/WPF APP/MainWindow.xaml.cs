using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using LIBRERIA_CLASES;
using Microsoft.Win32;
using System.IO;

using Excel = Microsoft.Office.Interop.Excel;
//using ExcelLibrary.CompoundDocumentFormat;
//using ExcelLibrary.SpreadSheet;

namespace WPF_APP
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        Tablas tablas = new Tablas();

        public List<Celda> listPolygons = new List<Celda>();

        double L = 65;
        int rows = 41;
        int columns = 120;

        Matriz matrix;


        public MainWindow()
        {
            InitializeComponent();
        }

        private void bt_Calculate_Click(object sender, RoutedEventArgs e)
        {
            canvas.Children.Clear();

            double M1 = Convert.ToDouble(tb_M.Text.Replace(Convert.ToChar("."), Convert.ToChar(",")));
            double P1 = Convert.ToDouble(tb_P.Text.Replace(Convert.ToChar("."), Convert.ToChar(",")));
            double rho1 = Convert.ToDouble(tb_rho.Text.Replace(Convert.ToChar("."), Convert.ToChar(",")));
            double T1 = Convert.ToDouble(tb_T.Text.Replace(Convert.ToChar("."), Convert.ToChar(",")));
            double theta = Convert.ToDouble(tb_theta.Text.Replace(Convert.ToChar("."), Convert.ToChar(",")));
            double gamma = Convert.ToDouble(tb_Gamma.Text.Replace(Convert.ToChar("."), Convert.ToChar(",")));
            double R = Convert.ToDouble(tb_R.Text.Replace(Convert.ToChar("."), Convert.ToChar(",")));

            matrix = new Matriz(rows, columns);
            double U1 = M1 * Math.Sqrt(R * gamma * T1);
            matrix.SetInitialConditions(M1, U1, 0, P1, rho1, T1, theta, gamma, R);
            matrix.CalculateSteps(theta, gamma, R);

            matrix.CalculatePolygons();

            double x = 0;
            for (int j = 0; j < columns && x < L; j++)
            {
                x = matrix.Matrix[0, j].x;
                for (int i = 0; i < rows - 1; i++)
                {
                    canvas.Children.Add(matrix.Matrix[i, j].poligono);

                    // añadimos la cela a una lista.
                    listPolygons.Add(matrix.Matrix[i, j]);
                }
            }
        }

        private void Polygon_MouseLeave(object sender, MouseEventArgs e)
        {

        }

        private void bt_RestartSimulation_Click(object sender, RoutedEventArgs e)
        {
            canvas.Children.Clear();
        }

        private void Polygon_MouseEnter(object sender, MouseEventArgs e)
        {
            var i = canvas.Children.IndexOf(sender as Polygon);

            tb_M1.Content = "M = " + Math.Round(listPolygons[i].M,4).ToString();
            tb_P1.Content = "P = " + Math.Round(listPolygons[i].P, 4).ToString(); ;
            tb_rho1.Content = "ρ = " + Math.Round(listPolygons[i].rho, 4).ToString(); ;
            tb_T1.Content = "T = " + Math.Round(listPolygons[i].T, 4).ToString(); ;
            tb_U1.Content = "U = " + Math.Round(listPolygons[i].U, 4).ToString(); ;
            tb_V1.Content = "V = " + Math.Round(listPolygons[i].V, 4).ToString(); ;
            tb_X1.Content = "X = " + Math.Round(listPolygons[i].x, 4).ToString(); ;
            tb_Y1.Content = "Y = " + Math.Round(listPolygons[i].y, 4).ToString(); ;
        }

        private void bt_M_Click(object sender, RoutedEventArgs e)
        {
            canvas.Children.Clear();

            // Buscamos valor maximo y minimo
            double max_M = listPolygons.Max(r => r.M);
            double min_M = listPolygons.Min(r => r.M);

            //double blanco = Convert.ToDouble(Convert.ToInt32("ffffff", 16));
            double blanco = Convert.ToDouble(Convert.ToInt32("ffffff", 16));
            double rojo = Convert.ToDouble(Convert.ToInt32("ff0000", 16));

            double m = (rojo - blanco) / (max_M - min_M);
            double n = blanco - m * min_M;

            // Recorremos todos los polígonos

            double x = 0;
            for (int j = 0; j < columns && x < L; j++)
            {
                x = matrix.Matrix[0, j].x;
                for (int i = 0; i < rows - 1; i++)
                {
                    // Ahora calculamos su color
                    double opacity = (matrix.Matrix[i, j].M - min_M) / max_M;

                    matrix.Matrix[i, j].poligono.Stroke = Brushes.LightGray;
                    matrix.Matrix[i, j].poligono.StrokeThickness = 1;

                    int color = Convert.ToInt32(Math.Floor(m * matrix.Matrix[i, j].M + n));
                    string color_hex = String.Concat("#", color.ToString("X"));
                    SolidColorBrush brush = new SolidColorBrush((Color)ColorConverter.ConvertFromString(color_hex));
                    matrix.Matrix[i, j].poligono.Fill = brush;

                    // añadimos los polígonos al canvas
                    canvas.Children.Add(matrix.Matrix[i, j].poligono);

                    // Añadimos eventos para poder ver el valor de las variables de este poligono/celfda si pasamos el mouse por encima
                    matrix.Matrix[i, j].poligono.MouseEnter += Polygon_MouseEnter;
                    matrix.Matrix[i, j].poligono.MouseLeave += Polygon_MouseLeave;
                }
            }

            // Actualizamos la gráficoa para esta magnitud

            List<double> axisX = new List<double>();
            List<double> axisY1 = new List<double>();
            List<double> axisY2 = new List<double>();

            x = 0;
            for (int j = 0; j < columns && x < L; j++)
            {
                List<double> listaAverage1 = new List<double>();

                x = matrix.Matrix[0, j].x;
                for (int i = 0; i < rows - 1; i++)
                {
                    listaAverage1.Add(matrix.Matrix[i, j].M);
                }
                axisY2.Add(matrix.Matrix[0, j].M);
                axisX.Add(x);
                axisY1.Add(listaAverage1.Average());
            }

            var X = axisX.ToArray();
            var Y1 = axisY1.ToArray();
            var Y2 = axisY2.ToArray();

            plot.plt.Clear();
            plot.plt.PlotScatter(X, Y1, markerSize:0, lineWidth: 3, label: "Average");
            plot.plt.PlotScatter(X, Y2, markerSize:0, lineWidth: 3, label: "Boundary");
            plot.plt.Legend();
        }

        private void bt_P_Click(object sender, RoutedEventArgs e)
        {
            canvas.Children.Clear();

            // Buscamos valor maximo y minimo
            double max_M = listPolygons.Max(r => r.P);
            double min_M = listPolygons.Min(r => r.P);

            //double blanco = Convert.ToDouble(Convert.ToInt32("ffffff", 16));
            double blanco = Convert.ToDouble(Convert.ToInt32("ffffff", 16));
            double rojo = Convert.ToDouble(Convert.ToInt32("ff0000", 16));

            double m = (rojo - blanco) / (max_M - min_M);
            double n = blanco - m * min_M;

            // Recorremos todos los polígonos

            double x = 0;
            for (int j = 0; j < columns && x < L; j++)
            {
                x = matrix.Matrix[0, j].x;
                for (int i = 0; i < rows - 1; i++)
                {
                    // Ahora calculamos su color
                    double opacity = (matrix.Matrix[i, j].P - min_M) / max_M;

                    matrix.Matrix[i, j].poligono.Stroke = Brushes.LightGray;
                    matrix.Matrix[i, j].poligono.StrokeThickness = 1;

                    int color = Convert.ToInt32(Math.Floor(m * matrix.Matrix[i, j].P + n));
                    string color_hex = String.Concat("#", color.ToString("X"));
                    SolidColorBrush brush = new SolidColorBrush((Color)ColorConverter.ConvertFromString(color_hex));
                    matrix.Matrix[i, j].poligono.Fill = brush;

                    // añadimos los polígonos al canvas
                    canvas.Children.Add(matrix.Matrix[i, j].poligono);

                    // Añadimos eventos para poder ver el valor de las variables de este poligono/celfda si pasamos el mouse por encima
                    matrix.Matrix[i, j].poligono.MouseEnter += Polygon_MouseEnter;
                    matrix.Matrix[i, j].poligono.MouseLeave += Polygon_MouseLeave;
                }
            }

            List<double> axisX = new List<double>();
            List<double> axisY1 = new List<double>();
            List<double> axisY2 = new List<double>();

            x = 0;
            for (int j = 0; j < columns && x < L; j++)
            {
                List<double> listaAverage1 = new List<double>();

                x = matrix.Matrix[0, j].x;
                for (int i = 0; i < rows - 1; i++)
                {
                    listaAverage1.Add(matrix.Matrix[i, j].P);
                }
                axisY2.Add(matrix.Matrix[0, j].P);
                axisX.Add(x);
                axisY1.Add(listaAverage1.Average());
            }

            var X = axisX.ToArray();
            var Y1 = axisY1.ToArray();
            var Y2 = axisY2.ToArray();

            plot.plt.Clear();
            plot.plt.PlotScatter(X, Y1, markerSize: 0, lineWidth: 3, label: "Average");
            plot.plt.PlotScatter(X, Y2, markerSize: 0, lineWidth: 3, label: "Boundary");
            plot.plt.Legend();
        }

        private void bt_Rho_Click(object sender, RoutedEventArgs e)
        {
            canvas.Children.Clear();

            // Buscamos valor maximo y minimo
            double max_M = listPolygons.Max(r => r.rho);
            double min_M = listPolygons.Min(r => r.rho);

            //double blanco = Convert.ToDouble(Convert.ToInt32("ffffff", 16));
            double blanco = Convert.ToDouble(Convert.ToInt32("ffffff", 16));
            double rojo = Convert.ToDouble(Convert.ToInt32("ff0000", 16));

            double m = (rojo - blanco) / (max_M - min_M);
            double n = blanco - m * min_M;

            // Recorremos todos los polígonos

            double x = 0;
            for (int j = 0; j < columns && x < L; j++)
            {
                x = matrix.Matrix[0, j].x;
                for (int i = 0; i < rows - 1; i++)
                {
                    // Ahora calculamos su color
                    double opacity = (matrix.Matrix[i, j].rho - min_M) / max_M;

                    matrix.Matrix[i, j].poligono.Stroke = Brushes.LightGray;
                    matrix.Matrix[i, j].poligono.StrokeThickness = 1;

                    int color = Convert.ToInt32(Math.Floor(m * matrix.Matrix[i, j].rho + n));
                    string color_hex = String.Concat("#", color.ToString("X"));
                    SolidColorBrush brush = new SolidColorBrush((Color)ColorConverter.ConvertFromString(color_hex));
                    matrix.Matrix[i, j].poligono.Fill = brush;

                    // añadimos los polígonos al canvas
                    canvas.Children.Add(matrix.Matrix[i, j].poligono);

                    // Añadimos eventos para poder ver el valor de las variables de este poligono/celfda si pasamos el mouse por encima
                    matrix.Matrix[i, j].poligono.MouseEnter += Polygon_MouseEnter;
                    matrix.Matrix[i, j].poligono.MouseLeave += Polygon_MouseLeave;
                }
            }

            List<double> axisX = new List<double>();
            List<double> axisY1 = new List<double>();
            List<double> axisY2 = new List<double>();

            x = 0;
            for (int j = 0; j < columns && x < L; j++)
            {
                List<double> listaAverage1 = new List<double>();

                x = matrix.Matrix[0, j].x;
                for (int i = 0; i < rows - 1; i++)
                {
                    listaAverage1.Add(matrix.Matrix[i, j].rho);
                }
                axisY2.Add(matrix.Matrix[0, j].rho);
                axisX.Add(x);
                axisY1.Add(listaAverage1.Average());
            }

            var X = axisX.ToArray();
            var Y1 = axisY1.ToArray();
            var Y2 = axisY2.ToArray();

            plot.plt.Clear();
            plot.plt.PlotScatter(X, Y1, markerSize: 0, lineWidth: 3, label: "Average");
            plot.plt.PlotScatter(X, Y2, markerSize: 0, lineWidth: 3, label: "Boundary");
            plot.plt.Legend();
        }

        private void bt_T_Click(object sender, RoutedEventArgs e)
        {
            canvas.Children.Clear();

            // Buscamos valor maximo y minimo
            double max_M = listPolygons.Max(r => r.T);
            double min_M = listPolygons.Min(r => r.T);

            //double blanco = Convert.ToDouble(Convert.ToInt32("ffffff", 16));
            double blanco = Convert.ToDouble(Convert.ToInt32("ffffff", 16));
            double rojo = Convert.ToDouble(Convert.ToInt32("ff0000", 16));

            double m = (rojo - blanco) / (max_M - min_M);
            double n = blanco - m * min_M;

            // Recorremos todos los polígonos

            double x = 0;
            for (int j = 0; j < columns && x < L; j++)
            {
                x = matrix.Matrix[0, j].x;
                for (int i = 0; i < rows - 1; i++)
                {
                    // Ahora calculamos su color
                    double opacity = (matrix.Matrix[i, j].T - min_M) / max_M;

                    matrix.Matrix[i, j].poligono.Stroke = Brushes.LightGray;
                    matrix.Matrix[i, j].poligono.StrokeThickness = 1;

                    int color = Convert.ToInt32(Math.Floor(m * matrix.Matrix[i, j].T + n));
                    string color_hex = String.Concat("#", color.ToString("X"));
                    SolidColorBrush brush = new SolidColorBrush((Color)ColorConverter.ConvertFromString(color_hex));
                    matrix.Matrix[i, j].poligono.Fill = brush;

                    // añadimos los polígonos al canvas
                    canvas.Children.Add(matrix.Matrix[i, j].poligono);

                    // Añadimos eventos para poder ver el valor de las variables de este poligono/celfda si pasamos el mouse por encima
                    matrix.Matrix[i, j].poligono.MouseEnter += Polygon_MouseEnter;
                    matrix.Matrix[i, j].poligono.MouseLeave += Polygon_MouseLeave;
                }
            }

            List<double> axisX = new List<double>();
            List<double> axisY1 = new List<double>();
            List<double> axisY2 = new List<double>();

            x = 0;
            for (int j = 0; j < columns && x < L; j++)
            {
                List<double> listaAverage1 = new List<double>();

                x = matrix.Matrix[0, j].x;
                for (int i = 0; i < rows - 1; i++)
                {
                    listaAverage1.Add(matrix.Matrix[i, j].T);
                }
                axisY2.Add(matrix.Matrix[0, j].T);
                axisX.Add(x);
                axisY1.Add(listaAverage1.Average());
            }

            var X = axisX.ToArray();
            var Y1 = axisY1.ToArray();
            var Y2 = axisY2.ToArray();

            plot.plt.Clear();
            plot.plt.PlotScatter(X, Y1, markerSize: 0, lineWidth: 3, label: "Average");
            plot.plt.PlotScatter(X, Y2, markerSize: 0, lineWidth: 3, label: "Boundary");
            plot.plt.Legend();
        }

        private void bt_U_Click(object sender, RoutedEventArgs e)
        {
            canvas.Children.Clear();

            // Buscamos valor maximo y minimo
            double max_M = listPolygons.Max(r => r.U);
            double min_M = listPolygons.Min(r => r.U);

            //double blanco = Convert.ToDouble(Convert.ToInt32("ffffff", 16));
            double blanco = Convert.ToDouble(Convert.ToInt32("ffffff", 16));
            double rojo = Convert.ToDouble(Convert.ToInt32("ff0000", 16));

            double m = (rojo - blanco) / (max_M - min_M);
            double n = blanco - m * min_M;

            // Recorremos todos los polígonos

            double x = 0;
            for (int j = 0; j < columns && x < L; j++)
            {
                x = matrix.Matrix[0, j].x;
                for (int i = 0; i < rows - 1; i++)
                {
                    // Ahora calculamos su color
                    double opacity = (matrix.Matrix[i, j].U - min_M) / max_M;

                    matrix.Matrix[i, j].poligono.Stroke = Brushes.LightGray;
                    matrix.Matrix[i, j].poligono.StrokeThickness = 1;

                    int color = Convert.ToInt32(Math.Floor(m * matrix.Matrix[i, j].U + n));
                    string color_hex = String.Concat("#", color.ToString("X"));
                    SolidColorBrush brush = new SolidColorBrush((Color)ColorConverter.ConvertFromString(color_hex));
                    matrix.Matrix[i, j].poligono.Fill = brush;

                    // añadimos los polígonos al canvas
                    canvas.Children.Add(matrix.Matrix[i, j].poligono);

                    // Añadimos eventos para poder ver el valor de las variables de este poligono/celfda si pasamos el mouse por encima
                    matrix.Matrix[i, j].poligono.MouseEnter += Polygon_MouseEnter;
                    matrix.Matrix[i, j].poligono.MouseLeave += Polygon_MouseLeave;
                }
            }

            List<double> axisX = new List<double>();
            List<double> axisY1 = new List<double>();
            List<double> axisY2 = new List<double>();

            x = 0;
            for (int j = 0; j < columns && x < L; j++)
            {
                List<double> listaAverage1 = new List<double>();

                x = matrix.Matrix[0, j].x;
                for (int i = 0; i < rows - 1; i++)
                {
                    listaAverage1.Add(matrix.Matrix[i, j].U);
                }
                axisY2.Add(matrix.Matrix[0, j].U);
                axisX.Add(x);
                axisY1.Add(listaAverage1.Average());
            }

            var X = axisX.ToArray();
            var Y1 = axisY1.ToArray();
            var Y2 = axisY2.ToArray();

            plot.plt.Clear();
            plot.plt.PlotScatter(X, Y1, markerSize: 0, lineWidth: 3, label: "Average");
            plot.plt.PlotScatter(X, Y2, markerSize: 0, lineWidth: 3, label: "Boundary");
            plot.plt.Legend();
        }

        private void bt_V_Click(object sender, RoutedEventArgs e)
        {
            canvas.Children.Clear();

            // Buscamos valor maximo y minimo
            double max_M = listPolygons.Max(r => r.V);
            double min_M = listPolygons.Min(r => r.V);

            //double blanco = Convert.ToDouble(Convert.ToInt32("ffffff", 16));
            double blanco = Convert.ToDouble(Convert.ToInt32("ffffff", 16));
            double rojo = Convert.ToDouble(Convert.ToInt32("ff0000", 16));

            double m = (rojo - blanco) / (max_M - min_M);
            double n = blanco - m * min_M;

            // Recorremos todos los polígonos

            double x = 0;
            for (int j = 0; j < columns && x < L; j++)
            {
                x = matrix.Matrix[0, j].x;
                for (int i = 0; i < rows - 1; i++)
                {
                    // Ahora calculamos su color
                    double opacity = (matrix.Matrix[i, j].V - min_M) / max_M;

                    matrix.Matrix[i, j].poligono.Stroke = Brushes.LightGray;
                    matrix.Matrix[i, j].poligono.StrokeThickness = 1;

                    int color = Convert.ToInt32(Math.Floor(m * matrix.Matrix[i, j].V + n));
                    string color_hex = String.Concat("#", color.ToString("X"));
                    SolidColorBrush brush = new SolidColorBrush((Color)ColorConverter.ConvertFromString(color_hex));
                    matrix.Matrix[i, j].poligono.Fill = brush;

                    // añadimos los polígonos al canvas
                    canvas.Children.Add(matrix.Matrix[i, j].poligono);

                    // Añadimos eventos para poder ver el valor de las variables de este poligono/celfda si pasamos el mouse por encima
                    matrix.Matrix[i, j].poligono.MouseEnter += Polygon_MouseEnter;
                    matrix.Matrix[i, j].poligono.MouseLeave += Polygon_MouseLeave;
                }
            }

            List<double> axisX = new List<double>();
            List<double> axisY1 = new List<double>();
            List<double> axisY2 = new List<double>();

            x = 0;
            for (int j = 0; j < columns && x < L; j++)
            {
                List<double> listaAverage1 = new List<double>();

                x = matrix.Matrix[0, j].x;
                for (int i = 0; i < rows - 1; i++)
                {
                    listaAverage1.Add(matrix.Matrix[i, j].V);
                }
                axisY2.Add(matrix.Matrix[0, j].V);
                axisX.Add(x);
                axisY1.Add(listaAverage1.Average());
            }

            var X = axisX.ToArray();
            var Y1 = axisY1.ToArray();
            var Y2 = axisY2.ToArray();

            plot.plt.Clear();
            plot.plt.PlotScatter(X, Y1, markerSize: 0, lineWidth: 3, label: "Average");
            plot.plt.PlotScatter(X, Y2, markerSize: 0, lineWidth: 3, label: "Boundary");
            plot.plt.Legend();
        }

        private void btn_Tables_Click(object sender, RoutedEventArgs e)
        {
            ShowTables(tablas);
        }

        public void ShowTables (Tablas tables)
        {
            this.Visibility = Visibility.Hidden;
            tablas.mainwindow = this;
            tablas.matrix = matrix;

            if (tables != null)
            {
                tables.Visibility = Visibility.Visible;
            }
            else
            {
                tablas.Show();
            }
        }

        private void Window_Loaded()
        {

        }
    }
}
