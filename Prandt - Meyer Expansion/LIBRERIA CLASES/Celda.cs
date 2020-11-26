using System;

namespace LIBRERIA_CLASES
{
    public class Celda
    {
        double x;
        double y;

        double M;
        double U;
        double V;
        double P;
        double rho;
        double T;
        double theta;
        double gamma;
        double R;

        public double F1, F2, F3, F4;
        public double G1, G2, G3, G4;

        public Celda(double M, double U, double V, double P, double rho, double T, double theta, double gamma, double R)
        {
            this.M = M;
            this.U = U;
            this.V = V;
            this.P = P;
            this.rho = rho;
            this.T = T;
            this.theta = theta;
            this.gamma = gamma;
            this.R = R;
        }
    }
}
