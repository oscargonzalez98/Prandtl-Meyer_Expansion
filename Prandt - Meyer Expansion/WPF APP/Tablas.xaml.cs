using DocumentFormat.OpenXml.Spreadsheet;
using LIBRERIA_CLASES;
using System;
using System.Collections.Generic;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Shapes;

namespace WPF_APP
{
    /// <summary>
    /// Lógica de interacción para Tablas.xaml
    /// </summary>
    public partial class Tablas : Window
    {
        public MainWindow mainwindow;
        List<string> it = new List<string>();
        public Matriz matrix;

        public Tablas()
        {
            InitializeComponent();

            DataGridTextColumn textColumn = new DataGridTextColumn();
            textColumn.Header = "Y";
            textColumn.Binding = new Binding("Y");
            datagrid.Columns.Add(textColumn);

            textColumn = new DataGridTextColumn();
            textColumn.Header = "M";
            textColumn.Binding = new Binding("M");
            datagrid.Columns.Add(textColumn);

            textColumn = new DataGridTextColumn();
            textColumn.Header = "P";
            textColumn.Binding = new Binding("P");
            datagrid.Columns.Add(textColumn);

            textColumn = new DataGridTextColumn();
            textColumn.Header = "ρ";
            textColumn.Binding = new Binding("ρ");
            datagrid.Columns.Add(textColumn);

            textColumn = new DataGridTextColumn();
            textColumn.Header = "T";
            textColumn.Binding = new Binding("T");
            datagrid.Columns.Add(textColumn);

            textColumn = new DataGridTextColumn();
            textColumn.Header = "U";
            textColumn.Binding = new Binding("U");
            datagrid.Columns.Add(textColumn);

            textColumn = new DataGridTextColumn();
            textColumn.Header = "V";
            textColumn.Binding = new Binding("V");
            datagrid.Columns.Add(textColumn);
        }

        private void btn_Simulation_Click(object sender, RoutedEventArgs e)
        {
            this.Visibility = Visibility.Hidden;
            mainwindow.Visibility = Visibility.Visible;
        }

        private void Grid_Loaded(object sender, RoutedEventArgs e)
        {

        }

        private void bt_M_Click(object sender, RoutedEventArgs e)
        {
            int rows = matrix.Matrix.GetLength(0);
            int columns = matrix.Matrix.GetLength(1);

            double x = matrix.Matrix[0, 0].x;
            for(int j = 0; x<65; j++)
            {
                it.Add(x.ToString());
                x = matrix.Matrix[0, j].x;
            }

            datagrid.ItemsSource = it;
        }
    }
}
