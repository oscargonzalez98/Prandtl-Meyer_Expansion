using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace LIBRERIA_CLASES
{
    public class Matriz
    {
        int rows;
        int columns;

        Celda[,] Matrix;

        double E = 10;
        double Theta = 5.352;
        double H = 40;
        double L = 65;
        double C = 0.5;

        double deltaETA = 0.025;
        

        public Matriz(int rows, int columns)
        {
            this.rows = rows;
            this.columns = columns;
            this.Matrix = new Celda[rows, columns];
        }

        public void SetInitialConditions(double M, double U, double V, double P, double rho, double T, double Theta, double gamma, double R)
        {
            // Rellenamos la primera columna de la matriz de celdas, asignandoles un valor incial
            for(int i = 0; i<rows; i++)
            {
                this.Matrix[i, 0] = new Celda(M, U, V, P, rho, T,Theta, gamma, R);
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

            // Calculamos los parametros del cambio de variables para cada celda de la columna inicial
            this.Calculate_ParametrosCambioVariable(E, H, Theta, 0);
        }

        public void Calculate_ParametrosCambioVariable(double E, double H, double Theta, int j)
        {
            // Calculamos ETA de cada celda de la columna

            for (int i = 1; i < rows; i++)
            {
                this.Matrix[i, j].ETA = this.Matrix[i - 1, j].ETA + deltaETA;
            }

            // Calculamos h, ys, de cada posicion de cada celda de la columna

            for (int i = 0; i < rows; i++)
            {
                double x = this.Matrix[i, j].x;
                double theta = Theta * Math.PI / 180;

                if (x < E)
                {
                    this.Matrix[i, j].h = H;
                    this.Matrix[i, j].y_s = 0;
                }
                else
                {
                    this.Matrix[i, j].h = H + (x - E) * Math.Tan(theta);
                    this.Matrix[i, j].y_s = - (x - E) * Math.Tan(theta);
                }

                this.Matrix[i,j].y = this.Matrix[i, j].h * this.Matrix[i, j].ETA + this.Matrix[i, j].y_s;
            }
        }

        public void RellenarColumnadeCeldas(int j)
        {
            for (int i = 0; i < rows; i++)
            {
                this.Matrix[i, j] = new Celda();
            }
        }

        public void Compute_Assign_DeltaXI(int j)
        {
            List<double> listaTanMax = new List<double>();

            for (int i = 0; i < rows; i++)
            {
                listaTanMax.Add(this.Matrix[i, j].CalculateTanMax());
            }
            double deltaY = this.Matrix[2, j].y - this.Matrix[1, j].y;
            double deltaXI = C * deltaY / listaTanMax.Max();

            // Asignamos el valor de x correspondiente de x en cada posicion de la columna de turno
            for (int i = 0; i < rows; i++)
            {
                this.Matrix[i, j + 1].x = this.Matrix[i, j].x + deltaXI;

            }
        }

        public void CalculateDeltaETA_DeltaX(int j)
        {
            for(int i = 0; i<rows; i++)
            {
                this.Matrix[i, j].CalculateDeltaETA_DeltaX(E);
            }
        }

        public void CalculatePRedictorStep(int j)
        {
            // Calculamos deltaF / deltaETA para cada posición de la columna. Eq 8.36, pag 387
            for(int i = 0; i < rows; i++)
            {
                if(i+1 < rows) // no estmaos en la foundary asi que hacemos forward difference
                {
                    this.Matrix[i, j].CalculateDeltaF_DeltaETA(this.Matrix[i + 1, j].F1, this.Matrix[i + 1, j].F2, this.Matrix[i + 1, j].F3, this.Matrix[i + 1, j].F4, this.Matrix[i + 1, j].G1, this.Matrix[i + 1, j].G2, this.Matrix[i + 1, j].G3, this.Matrix[i + 1, j].G4, deltaETA);
                }
                else // estamos en la boundary, asi que hacemos una backward difference
                {
                    this.Matrix[i, j].CalculateDeltaF_DeltaETA_Boundary(this.Matrix[i - 1, j].F1, this.Matrix[i - 1, j].F2, this.Matrix[i - 1, j].F3, this.Matrix[i - 1, j].F4, this.Matrix[i - 1, j].G1, this.Matrix[i - 1, j].G2, this.Matrix[i - 1, j].G3, this.Matrix[i - 1, j].G4, deltaETA);
                }
            }

            // Calculamos las Fs predecidas en la siguiente columna Eq 8.37, pag 389
            for(int i = 0; i < rows; i++)
            {
                if( i == 0) // estamos en una boundary, no tenemos en cuenta la friccion artificial
                {
                    this.Matrix[i, j + 1].CalculatePredictedF_NextColumn_Boundary(this.Matrix[i, j].F1, this.Matrix[i, j].F2, this.Matrix[i, j].F3, this.Matrix[i, j].F4, this.Matrix[i, j].deltaF1_deltaETA, this.Matrix[i, j].deltaF2_deltaETA, this.Matrix[i, j].deltaF3_deltaETA, this.Matrix[i, j].deltaF4_deltaETA, deltaETA);
                }
                else if(i+1 == rows)
                {
                    this.Matrix[i, j + 1].CalculatePredictedF_NextColumn_Boundary(this.Matrix[i, j].F1, this.Matrix[i, j].F2, this.Matrix[i, j].F3, this.Matrix[i, j].F4, this.Matrix[i, j].deltaF1_deltaETA, this.Matrix[i, j].deltaF2_deltaETA, this.Matrix[i, j].deltaF3_deltaETA, this.Matrix[i, j].deltaF4_deltaETA, deltaETA);
                }
                else
                {
                    this.Matrix[i, j + 1].CalculatePredictedF_NextColumn(this.Matrix[i, j].F1, this.Matrix[i + 1, j].F1, this.Matrix[i - 1, j].F1, this.Matrix[i, j].F2, this.Matrix[i + 1, j].F2, this.Matrix[i - 1, j].F2, this.Matrix[i, j].F3, this.Matrix[i + 1, j].F3, this.Matrix[i - 1, j].F3, this.Matrix[i, j].F4, this.Matrix[i + 1, j].F4, this.Matrix[i - 1, j].F4, this.Matrix[i, j].deltaF1_deltaETA, this.Matrix[i, j].deltaF2_deltaETA, this.Matrix[i, j].deltaF3_deltaETA, this.Matrix[i, j].deltaF4_deltaETA, deltaETA, this.Matrix[i, j].P, this.Matrix[i + 1, j].P, this.Matrix[i - 1, j].P);
                }
            }

            // Calculamos las Gs predecidas en la siguiente columna
            for(int i = 0; i<rows; i++)
            {
                this.Matrix[i, j + 1].CalculatePredictedG_NextColumn();
            }
        }

        public void CalculateCorrectorStep(int j)
        {
            // Ahora, con los valores de F y G predecidos calculamos la derivada de F respecto de Xi predecida 
            for (int i = 0; i < rows; i++)
            {
                if (i == 0)
                {
                    this.Matrix[i, j + 1].CalculatePRedicted_dF_dETA_NextColumn_Boundary(this.Matrix[i + 1, j + 1].predictedF1, this.Matrix[i + 1, j + 1].predictedF2, this.Matrix[i + 1, j + 1].predictedF3, this.Matrix[i + 1, j + 1].predictedF4, this.Matrix[i + 1, j + 1].predictedG1, this.Matrix[i + 1, j + 1].predictedG2, this.Matrix[i + 1, j + 1].predictedG3, this.Matrix[i + 1, j + 1].predictedG4, deltaETA);
                }
                else
                {
                    this.Matrix[i, j + 1].CalculatePRedicted_dF_dETA_NextColumn(this.Matrix[i - 1, j + 1].predictedF1, this.Matrix[i - 1, j + 1].predictedF2, this.Matrix[i - 1, j + 1].predictedF3, this.Matrix[i - 1, j + 1].predictedF4, this.Matrix[i - 1, j + 1].predictedG1, this.Matrix[i - 1, j + 1].predictedG2, this.Matrix[i - 1, j + 1].predictedG3, this.Matrix[i - 1, j + 1].predictedG4, deltaETA);
                }
            }

            // Ahora calculamos dF / dETA average
            for(int i = 0; i < rows; i++)
            {
                this.Matrix[i, j].CalculateDeltaF_DeltaETA_Average(this.Matrix[i,j + 1].predicted_deltaF1_deltaETA, this.Matrix[i, j + 1].predicted_deltaF2_deltaETA, this.Matrix[i, j + 1].predicted_deltaF3_deltaETA, this.Matrix[i, j + 1].predicted_deltaF4_deltaETA);
            }

            // Ahora calculamos las Fs (las de verdad) de cada celda de la columna siguiente

            // Primero calculamos Ppredecida en cada celda de la columna siguiente (la necesitamos para el calculo de SF)
            for(int i = 0; i < rows; i++)
            {
                this.Matrix[i, j].CalculatePredictedP();
            }

            for(int i = 0; i< rows; i++)
            {
                if(i == 0 || (i+1 == rows))
                {
                    this.Matrix[i, j + 1].CalculateF_NextColumn_Boundary(this.Matrix[i, j].F1, this.Matrix[i, j].F2, this.Matrix[i, j].F3, this.Matrix[i, j].F4, this.Matrix[i, j].deltaF1_deltaETA_Average, this.Matrix[i, j].deltaF2_deltaETA_Average, this.Matrix[i, j].deltaF3_deltaETA_Average, this.Matrix[i, j].deltaF4_deltaETA_Average, deltaETA);
                }
                else
                {
                    this.Matrix[i, j + 1].CalculateF_NextColumn(this.Matrix[i, j].F1, this.Matrix[i + 1, j + 1].predictedF1, this.Matrix[i - 1, j + 1].predictedF1, this.Matrix[i, j].F2,
                        this.Matrix[i + 1, j + 1].predictedF2, this.Matrix[i - 1, j + 1].predictedF3, this.Matrix[i, j].F3, this.Matrix[i + 1, j + 1].predictedF3,
                        this.Matrix[i - 1, j + 1].predictedF3, this.Matrix[i, j].F4, this.Matrix[i + 1, j + 1].predictedF4, this.Matrix[i - 1, j + 1].predictedF4, 
                        this.Matrix[i, j].deltaF1_deltaETA_Average, this.Matrix[i, j].deltaF2_deltaETA_Average, this.Matrix[i, j].deltaF3_deltaETA_Average, this.Matrix[i, j].deltaF4_deltaETA_Average, 
                        this.Matrix[i + 1, j + 1].predictedP, this.Matrix[i - 1, j + 1].predictedP, deltaETA);
                }
            }
        }



        public void CalculateSteps()
        {
            // Bucle que recorra todas las columnas

            double x = this.Matrix[0, 0].x;
            int j = 0;
            while((x < L) && (j + 1) < columns)
            {
                // Cremos las celdas de la nueva columna
                this.RellenarColumnadeCeldas(j + 1);

                // Calculamos deltaX entre la columna anterior y esta
                this.Compute_Assign_DeltaXI(j);

                // Calculamos deltaETA / deltaX
                this.CalculateDeltaETA_DeltaX(j);

                // Calculamos el Predictor Step:
                this.CalculatePRedictorStep(j);

                // Calculamos el corrector Step
                this.CalculateCorrectorStep(j);



                j++;
            }
        }


    }
}
