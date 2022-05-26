import static java.lang.Math.*;
import static java.lang.Math.cos;

public class FindCountPursuersBetaAboveZero {
    public static double[] firstPursuer(double r, double l, double alpha, double beta) {
        double a = (pow(alpha, 2) + pow(beta, 2));
        double b = 2 * beta * (r - l);
        double c = -4 * r * l;
        double d = pow(b, 2) - 4 * a * c;
        double t1 = (-b + sqrt(d)) / (2 * a); // positive value
        double t2 = (-b - sqrt(d)) / (2 * a); // negative value
        double phi = asin((r + beta * t1 - l) / (r + l)); // phi в радианах
        double var = pow(alpha * t1, 2) + pow(r + beta * t1 - l, 2);
//        System.out.println("var = " + var);
//        System.out.println("(r+l)^2 = " + pow(r + l, 2));
        return new double[]{t1, phi};
    }


    public static int findCountPursuersBZ(double r, double l, double alpha, double beta) {
        int k;
        if (r + beta * (r + l) / alpha < l) {
            k = 1;
        } else {
            k = 2;

        }
        return k;
    }

    public static void main(String[] args) {
        double r = 80;
        double l = 50;
        double alpha = 2;
        double beta = 1;
        double t1 = firstPursuer(r, l, alpha, beta)[0];
        double phi = firstPursuer(r, l, alpha, beta)[1];
        System.out.println("t1 = " + t1);
        double phiInDeg = Math.toDegrees(phi); // Угол в градусах
        System.out.println("phi = " + phiInDeg);

    }
}
