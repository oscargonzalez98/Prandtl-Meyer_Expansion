using System;
using System.Collections.Generic;
using System.Text;

namespace LIBRERIA_CLASES
{
    public class Matriz
    {
        int rows;
        int columns;

        Celda[,] Matrix;
        

        public Matriz(int rows, int columns)
        {
            this.rows = rows;
            this.columns = columns;
            this.Matrix = new Celda[rows, columns];
        }

        public void SetInitialConditions(double M, double U, double V, double P, double rho, double T, double theta, double gamma, double R)
        {
            // Rellenamos la primera columna de la matriz de celdas, asignandoles un valor incial
            for(int i = 0; i<rows; i++)
            {
                this.Matrix[i, 0] = new Celda(M, U, V, P, rho, T, theta, gamma, R);
            }

            // Calculamos F's y G's de la primera columna
            for (int i = 0; i < rows; i++)
            {
                this.Matrix[i, 0].F1 = rho * U;
                this.Matrix[i, 0].F2 = rho * U * U + P;
                this.Matrix[i, 0].F3 = rho * U * V;
                this.Matrix[i, 0].F4 = (gamma / (gamma - 1)) * P * U + rho * U * (((U * U) + (V * V)) / 2);

                this.Matrix[i, 0].G1 = rho * V;
                this.Matrix[i, 0].G2 = rho * U * V;
                this.Matrix[i, 0].G3 = rho * V * V + P;
                this.Matrix[i, 0].G4 = (gamma / (gamma - 1)) * P * V + rho * V * (((U * U) + (V * V)) / 2);
            }
        }


    }
}
