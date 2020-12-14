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

namespace WPF_APP
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public List<Celda> listPolygons = new List<Celda>();

        public MainWindow()
        {
            InitializeComponent();
        }

        private void bt_Calculate_Click(object sender, RoutedEventArgs e)
        {
            int rows = 41;
            int columns = 120;

            double M1 = Convert.ToDouble(tb_M.Text.Replace(Convert.ToChar("."), Convert.ToChar(",")));
            double P1 = Convert.ToDouble(tb_P.Text.Replace(Convert.ToChar("."), Convert.ToChar(",")));
            double rho1 = Convert.ToDouble(tb_rho.Text.Replace(Convert.ToChar("."), Convert.ToChar(",")));
            double T1 = Convert.ToDouble(tb_T.Text.Replace(Convert.ToChar("."), Convert.ToChar(",")));
            double theta = Convert.ToDouble(tb_theta.Text.Replace(Convert.ToChar("."), Convert.ToChar(",")));
            double gamma = Convert.ToDouble(tb_Gamma.Text.Replace(Convert.ToChar("."), Convert.ToChar(",")));
            double R = Convert.ToDouble(tb_R.Text.Replace(Convert.ToChar("."), Convert.ToChar(",")));

            Matriz matrix = new Matriz(rows, columns);
            double U1 = M1 * Math.Sqrt(R * gamma * T1);
            matrix.SetInitialConditions(M1, U1, 0, P1, rho1, T1, theta, gamma, R);
            matrix.CalculateSteps(theta, gamma, R);

            matrix.CalculatePolygons();

            for (int j = 0; j < columns; j++)
            {
                for (int i = 0; i < rows - 1; i++)
                {
                    canvas.Children.Add(matrix.Matrix[i, j].poligono);

                    matrix.Matrix[i, j].poligono.MouseEnter += Polygon_MouseEnter;
                    matrix.Matrix[i, j].poligono.MouseLeave += Polygon_MouseLeave;

                    listPolygons.Add(matrix.Matrix[i, j]);
                }
            }
        }

        private void Polygon_MouseLeave(object sender, MouseEventArgs e)
        {

        }

        private void Polygon_MouseEnter(object sender, MouseEventArgs e)
        {
            var i = canvas.Children.IndexOf(sender as Polygon);

            tb_M1.Content = "M = " + listPolygons[i].M.ToString();
            tb_P1.Content = "P = " + listPolygons[i].P.ToString();
            tb_rho1.Content = "ρ = " + listPolygons[i].rho.ToString();
            tb_T1.Content = "T = " + listPolygons[i].T.ToString();
            tb_U1.Content = "U = " + listPolygons[i].U.ToString();
            tb_V1.Content = "V = " + listPolygons[i].V.ToString();
            tb_X1.Content = "X = " + listPolygons[i].x.ToString();
            tb_Y1.Content = "Y = " + listPolygons[i].y.ToString();
        }
    }
}
