using System;
using System.Collections.Generic;
using System.Linq;

namespace LIBRERIA_CLASES
{
    public class Celda
    {
        double Cy = 0.6;

        public double x;
        public double y;

        public double M;
        public double U;
        public double V;
        public double P;
        public double rho;
        public double T;
        public double Theta;
        public double gamma;
        public double R;

        public double ETA;
        public double XI;

        public double h;
        public double y_s;

        public double dETA_dX;

        public double F1, F2, F3, F4;
        public double G1, G2, G3, G4;

        public double predictedF1, predictedF2, predictedF3, predictedF4;
        public double predictedG1, predictedG2, predictedG3, predictedG4;

        public double deltaF1_deltaETA, deltaF2_deltaETA, deltaF3_deltaETA, deltaF4_deltaETA;
        public double predicted_deltaF1_deltaETA, predicted_deltaF2_deltaETA, predicted_deltaF3_deltaETA, predicted_deltaF4_deltaETA;
        public double deltaF1_deltaETA_Average, deltaF2_deltaETA_Average, deltaF3_deltaETA_Average, deltaF4_deltaETA_Average;

        public Celda(double M, double U, double V, double P, double rho, double T, double Theta, double gamma, double R)
        {
            this.M = M;
            this.U = U;
            this.V = V;
            this.P = P;
            this.rho = rho;
            this.T = T;
            this.Theta = Theta;
            this.gamma = gamma;
            this.R = R;
        }

        public Celda()
        {

        }

        public double CalculateTanMax()
        {
            double theta = Theta * Math.PI / 180;
            double mu = Math.Asin(1 / this.M);

            List<double> lista = new List<double>();
            lista.Add(Math.Abs(Math.Tan(theta + mu)));
            lista.Add(Math.Abs(Math.Tan(theta - mu)));

            return lista.Max();
        }

        public void CalculateDeltaETA_DeltaX(double E)
        {
            double theta = Theta * Math.PI / 180;

            if(this.x < E)
            {
                dETA_dX = 0;
            }
            else
            {
                dETA_dX = (1 - ETA) * Math.Tan(theta) / this.h;
            }
        }

        public void CalculateDeltaF_DeltaETA(double F1_arriba, double F2_arriba, double F3_arriba, double F4_arriba, double G1_arriba, double G2_arriba, double G3_arriba, double G4_arriba, double DeltaETA)
        {
            this.deltaF1_deltaETA = dETA_dX * ((this.F1 - F1_arriba) / DeltaETA) + (1 / this.h) * ((this.G1 - G1_arriba) / DeltaETA);
            this.deltaF2_deltaETA = dETA_dX * ((this.F2 - F2_arriba) / DeltaETA) + (1 / this.h) * ((this.G2 - G2_arriba) / DeltaETA);
            this.deltaF3_deltaETA = dETA_dX * ((this.F3 - F3_arriba) / DeltaETA) + (1 / this.h) * ((this.G3 - G3_arriba) / DeltaETA);
            this.deltaF4_deltaETA = dETA_dX * ((this.F4 - F4_arriba) / DeltaETA) + (1 / this.h) * ((this.G4 - G4_arriba) / DeltaETA);
        }
        public void CalculateDeltaF_DeltaETA_Boundary(double F1_abajo, double F2_abajo, double F3_abajo, double F4_abajo, double G1_abajo, double G2_abajo, double G3_abajo, double G4_abajo, double DeltaETA)
        {
            this.deltaF1_deltaETA = dETA_dX * ((F1_abajo - this.F1) / DeltaETA) + (1 / this.h) * ((G1_abajo - this.G1) / DeltaETA);
            this.deltaF2_deltaETA = dETA_dX * ((F2_abajo - this.F2) / DeltaETA) + (1 / this.h) * ((G2_abajo - this.G2) / DeltaETA);
            this.deltaF3_deltaETA = dETA_dX * ((F3_abajo - this.F3) / DeltaETA) + (1 / this.h) * ((G3_abajo - this.G3) / DeltaETA);
            this.deltaF4_deltaETA = dETA_dX * ((F4_abajo - this.F4) / DeltaETA) + (1 / this.h) * ((G4_abajo - this.G4) / DeltaETA);
        }

        public void CalculatePredictedF_NextColumn(double F1, double F1_arriba, double F1_abajo, double F2, double F2_arriba, double F2_abajo, double F3, double F3_arriba, double F3_abajo, double F4, double F4_arriba, double F4_abajo, double dF1_dETA, double dF2_dETA, double dF3_dETA, double dF4_dETA,double p, double p_arriba, double p_abajo,  double deltaETA)
        {
            double SF1 = (Cy * Math.Abs(p_arriba - 2 * p + p_abajo) / (p_arriba + 2 * p + p_abajo)) * (F1_arriba - 2 * F1 + F1_abajo);
            double SF2 = (Cy * Math.Abs(p_arriba - 2 * p + p_abajo) / (p_arriba + 2 * p + p_abajo)) * (F2_arriba - 2 * F2 + F2_abajo);
            double SF3 = (Cy * Math.Abs(p_arriba - 2 * p + p_abajo) / (p_arriba + 2 * p + p_abajo)) * (F3_arriba - 2 * F3 + F3_abajo);
            double SF4 = (Cy * Math.Abs(p_arriba - 2 * p + p_abajo) / (p_arriba + 2 * p + p_abajo)) * (F4_arriba - 2 * F4 + F4_abajo);

            this.predictedF1 = F1 + dF1_dETA * deltaETA + SF1;
            this.predictedF2 = F2 + dF2_dETA * deltaETA + SF2;
            this.predictedF3 = F3 + dF3_dETA * deltaETA + SF3;
            this.predictedF4 = F4 + dF4_dETA * deltaETA + SF4;
        }
        public void CalculatePredictedF_NextColumn_Boundary(double F1, double F2, double F3, double F4, double dF1_dETA, double dF2_dETA, double dF3_dETA, double dF4_dETA, double deltaETA)
        {
            this.predictedF1 = F1 + dF1_dETA * deltaETA;
            this.predictedF2 = F2 + dF2_dETA * deltaETA;
            this.predictedF3 = F3 + dF3_dETA * deltaETA;
            this.predictedF4 = F4 + dF4_dETA * deltaETA;
        }

        public void CalculatePredictedG_NextColumn()
        {
            double A = (Math.Pow(this.predictedF3, 2) / (2 * this.predictedF1)) - this.predictedF4;
            double B = (this.gamma / (this.gamma - 1)) * this.predictedF1 * this.predictedF2;
            double C = -1 * ((this.gamma + 1) / (2 * (this.gamma - 1))) * Math.Pow(this.predictedF3, 3);

            double predictedRho = ((-B + Math.Sqrt(B * B - 4 * A * C)) / 2 * A);

            this.predictedG1 = predictedRho * (this.predictedF3 / this.predictedF1);
            this.predictedG2 = this.predictedF3;
            this.predictedG3 = predictedRho * Math.Pow(this.predictedF3 / this.predictedF1, 2) + this.predictedF2 - (Math.Pow(this.predictedF1, 2) / predictedRho);
            this.predictedG4 = (this.gamma / (this.gamma - 1)) * (this.predictedF2 - (Math.Pow(this.predictedF1, 2) / predictedRho)) * (this.predictedF3 / this.predictedF1) + (predictedRho / 2) * (this.predictedF3 / this.predictedF1) * (Math.Pow(this.predictedF1 / predictedRho, 2) + Math.Pow(this.predictedF3 / this.predictedF1, 2));
        }

        public void CalculatePRedicted_dF_dETA_NextColumn_Boundary(double predictedF1_arriba, double predictedF2_arriba, double predictedF3_arriba, double predictedF4_arriba, double predictedG1_arriba, double predictedG2_arriba, double predictedG3_arriba, double predictedG4_arriba, double deltaETA)
        {
            this.predicted_deltaF1_deltaETA = dETA_dX * ((this.predictedF1 - predictedF1_arriba) / deltaETA) + (1 / this.h) * ((this.predictedG1 - predictedG1_arriba) / deltaETA);
            this.predicted_deltaF2_deltaETA = dETA_dX * ((this.predictedF2 - predictedF2_arriba) / deltaETA) + (1 / this.h) * ((this.predictedG2 - predictedG2_arriba) / deltaETA);
            this.predicted_deltaF3_deltaETA = dETA_dX * ((this.predictedF3 - predictedF3_arriba) / deltaETA) + (1 / this.h) * ((this.predictedG3 - predictedG3_arriba) / deltaETA);
            this.predicted_deltaF4_deltaETA = dETA_dX * ((this.predictedF4 - predictedF4_arriba) / deltaETA) + (1 / this.h) * ((this.predictedG4 - predictedG4_arriba) / deltaETA);
        }

        public void CalculatePRedicted_dF_dETA_NextColumn(double predictedF1_abajo, double predictedF2_abajo, double predictedF3_abajo, double predictedF4_abajo, double predictedG1_abajo, double predictedG2_abajo, double predictedG3_abajo, double predictedG4_abajo, double deltaETA)
        {
            this.predicted_deltaF1_deltaETA = dETA_dX * ((predictedF1_abajo - this.predictedF1) / deltaETA) + (1 / this.h) * ((predictedG1_abajo - this.predictedG1) / deltaETA);
            this.predicted_deltaF2_deltaETA = dETA_dX * ((predictedF2_abajo - this.predictedF2) / deltaETA) + (1 / this.h) * ((predictedG2_abajo - this.predictedG2) / deltaETA);
            this.predicted_deltaF3_deltaETA = dETA_dX * ((predictedF3_abajo - this.predictedF3) / deltaETA) + (1 / this.h) * ((predictedG3_abajo - this.predictedG3) / deltaETA);
            this.predicted_deltaF4_deltaETA = dETA_dX * ((predictedF4_abajo - this.predictedF4) / deltaETA) + (1 / this.h) * ((predictedG4_abajo - this.predictedG4) / deltaETA);
        }

        public void CalculateDeltaF_DeltaETA_Average(double predicted_deltaF1_deltaETA_Siguiente, double predicted_deltaF2_deltaETA_Siguiente, double predicted_deltaF3_deltaETA_Siguiente, double predicted_deltaF4_deltaETA_Siguiente)
        {
            this.deltaF1_deltaETA_Average = (1 / 2) * (this.deltaF1_deltaETA + predicted_deltaF1_deltaETA_Siguiente);
            this.deltaF2_deltaETA_Average = (1 / 2) * (this.deltaF2_deltaETA + predicted_deltaF2_deltaETA_Siguiente);
            this.deltaF3_deltaETA_Average = (1 / 2) * (this.deltaF3_deltaETA + predicted_deltaF3_deltaETA_Siguiente);
            this.deltaF4_deltaETA_Average = (1 / 2) * (this.deltaF4_deltaETA + predicted_deltaF4_deltaETA_Siguiente);
        }
    }
}
