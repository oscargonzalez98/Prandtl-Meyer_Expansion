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

            Matriz Matrix = new Matriz(rows, columns);
            double U1 = M1 * Math.Sqrt(R * gamma * T1);
            Matrix.SetInitialConditions(M1, U1, 0, P1, rho1, T1, theta, gamma, R);
            Matrix.CalculateSteps();
        }
    }
}
