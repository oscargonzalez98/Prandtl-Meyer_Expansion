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
            int rows = 40;
            int columns = 120;

            double M = Convert.ToDouble(tb_M.Text.Replace(Convert.ToChar("."), Convert.ToChar(",")));
            double U = Convert.ToDouble(tb_U.Text.Replace(Convert.ToChar("."), Convert.ToChar(",")));
            double P = Convert.ToDouble(tb_P.Text.Replace(Convert.ToChar("."), Convert.ToChar(",")));
            double rho = Convert.ToDouble(tb_rho.Text.Replace(Convert.ToChar("."), Convert.ToChar(",")));
            double T = Convert.ToDouble(tb_T.Text.Replace(Convert.ToChar("."), Convert.ToChar(",")));
            double theta = Convert.ToDouble(tb_theta.Text.Replace(Convert.ToChar("."), Convert.ToChar(",")));


        }
    }
}
