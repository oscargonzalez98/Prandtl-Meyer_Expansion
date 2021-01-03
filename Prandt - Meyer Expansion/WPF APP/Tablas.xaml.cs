using DocumentFormat.OpenXml.Spreadsheet;
using LIBRERIA_CLASES;
using System;
using System.Collections.Generic;
using System.Data;
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
        public Matriz matrix;

        DataTable dt_M = new DataTable();
        DataTable dt_P = new DataTable();
        DataTable dt_Rho = new DataTable();
        DataTable dt_T = new DataTable();
        DataTable dt_U = new DataTable();
        DataTable dt_V = new DataTable();

        public Tablas(Matriz matrix)
        {
            InitializeComponent();
            this.matrix = matrix;
        }

        private void btn_Simulation_Click(object sender, RoutedEventArgs e)
        {
            this.Visibility = Visibility.Hidden;
            mainwindow.Visibility = Visibility.Visible;
        }

        private void bt_M_Click(object sender, RoutedEventArgs e)
        {
            int rows = matrix.Matrix.GetLength(0);
            int columns = matrix.Matrix.GetLength(1);

            // Primero añadimos las columnas
            double x = matrix.Matrix[0, 0].x;
            int j = 0;

            while (x < 65 && j < columns)
            {
                DataColumn row = new DataColumn(Convert.ToString(""), typeof(double));
                dt_M.Columns.Add(row);
                x = matrix.Matrix[0, j].x;
                j++;
            }


            for (int i = 0; i< rows; i++)
            {
                DataRow row = dt_M.NewRow();
                x = matrix.Matrix[0, 0].x;
                j = 0;

                while (x < 65 && j < columns)
                {
                    row[j] = Math.Round(this.matrix.Matrix[i, j].M, 4);
                    x = matrix.Matrix[0, j].x;
                    j++;
                }
                dt_M.Rows.Add(row);
            }

            datagrid.ItemsSource = dt_M.DefaultView;
        }

        private void bt_P_Click(object sender, RoutedEventArgs e)
        {
            int rows = matrix.Matrix.GetLength(0);
            int columns = matrix.Matrix.GetLength(1);

            // Primero añadimos las columnas
            double x = matrix.Matrix[0, 0].x;
            int j = 0;

            while (x < 65 && j < columns)
            {
                DataColumn row = new DataColumn(Convert.ToString(""), typeof(double));
                dt_P.Columns.Add(row);
                x = matrix.Matrix[0, j].x;
                j++;
            }


            for (int i = 0; i < rows; i++)
            {
                DataRow row = dt_P.NewRow();
                x = matrix.Matrix[0, 0].x;
                j = 0;

                while (x < 65 && j < columns)
                {
                    row[j] = Math.Round(this.matrix.Matrix[i, j].P, 4);
                    x = matrix.Matrix[0, j].x;
                    j++;
                }
                dt_P.Rows.Add(row);
            }

            datagrid.ItemsSource = dt_P.DefaultView;
        }

        private void bt_Rho_Click(object sender, RoutedEventArgs e)
        {
            int rows = matrix.Matrix.GetLength(0);
            int columns = matrix.Matrix.GetLength(1);

            // Primero añadimos las columnas
            double x = matrix.Matrix[0, 0].x;
            int j = 0;

            while (x < 65 && j < columns)
            {
                DataColumn row = new DataColumn(Convert.ToString(""), typeof(double));
                dt_Rho.Columns.Add(row);
                x = matrix.Matrix[0, j].x;
                j++;
            }


            for (int i = 0; i < rows; i++)
            {
                DataRow row = dt_Rho.NewRow();
                x = matrix.Matrix[0, 0].x;
                j = 0;

                while (x < 65 && j < columns)
                {
                    row[j] = Math.Round(this.matrix.Matrix[i, j].rho, 4);
                    x = matrix.Matrix[0, j].x;
                    j++;
                }
                dt_Rho.Rows.Add(row);
            }

            datagrid.ItemsSource = dt_Rho.DefaultView;
        }

        private void bt_T_Click(object sender, RoutedEventArgs e)
        {
            int rows = matrix.Matrix.GetLength(0);
            int columns = matrix.Matrix.GetLength(1);

            // Primero añadimos las columnas
            double x = matrix.Matrix[0, 0].x;
            int j = 0;

            while (x < 65 && j < columns)
            {
                DataColumn row = new DataColumn(Convert.ToString(""), typeof(double));
                dt_T.Columns.Add(row);
                x = matrix.Matrix[0, j].x;
                j++;
            }


            for (int i = 0; i < rows; i++)
            {
                DataRow row = dt_T.NewRow();
                x = matrix.Matrix[0, 0].x;
                j = 0;

                while (x < 65 && j < columns)
                {
                    row[j] = Math.Round(this.matrix.Matrix[i, j].T, 4);
                    x = matrix.Matrix[0, j].x;
                    j++;
                }
                dt_T.Rows.Add(row);
            }

            datagrid.ItemsSource = dt_T.DefaultView;
        }

        private void bt_U_Click(object sender, RoutedEventArgs e)
        {
            int rows = matrix.Matrix.GetLength(0);
            int columns = matrix.Matrix.GetLength(1);

            // Primero añadimos las columnas
            double x = matrix.Matrix[0, 0].x;
            int j = 0;

            while (x < 65 && j < columns)
            {
                DataColumn row = new DataColumn(Convert.ToString(""), typeof(double));
                dt_U.Columns.Add(row);
                x = matrix.Matrix[0, j].x;
                j++;
            }


            for (int i = 0; i < rows; i++)
            {
                DataRow row = dt_U.NewRow();
                x = matrix.Matrix[0, 0].x;
                j = 0;

                while (x < 65 && j < columns)
                {
                    row[j] = Math.Round(this.matrix.Matrix[i, j].U, 4);
                    x = matrix.Matrix[0, j].x;
                    j++;
                }
                dt_U.Rows.Add(row);
            }

            datagrid.ItemsSource = dt_U.DefaultView;
        }

        private void bt_V_Click(object sender, RoutedEventArgs e)
        {
            int rows = matrix.Matrix.GetLength(0);
            int columns = matrix.Matrix.GetLength(1);

            // Primero añadimos las columnas
            double x = matrix.Matrix[0, 0].x;
            int j = 0;

            while (x < 65 && j < columns)
            {
                DataColumn row = new DataColumn(Convert.ToString(""), typeof(double));
                dt_V.Columns.Add(row);
                x = matrix.Matrix[0, j].x;
                j++;
            }

            for (int i = 0; i < rows; i++)
            {
                DataRow row = dt_V.NewRow();
                x = matrix.Matrix[0, 0].x;
                j = 0;

                while (x < 65 && j < columns)
                {
                    row[j] = Math.Round(this.matrix.Matrix[i, j].V, 4);
                    x = matrix.Matrix[0, j].x;
                    j++;
                }
                dt_V.Rows.Add(row);
            }

            datagrid.ItemsSource = dt_V.DefaultView;
        }
    }
}
