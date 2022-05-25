package sample.others;

import static java.lang.Math.*;

public class CalcK {

    /**
     * Метод для нахождения угла по двум точкам
     *
     * @param Ax
     * @param Ay
     * @param Bx
     * @param By
     * @return
     */
    private static double findArc(double Ax, double Ay, double Bx, double By) {
        double arc = 0;
        if ((By > Ay) && (Bx > Ax)) arc = atan((By - Ay) / (Bx - Ax));
        if ((By < Ay) && (Bx > Ax)) arc = 2 * PI - atan((Ay - By) / (Bx - Ax));
        if ((By > Ay) && (Bx < Ax)) arc = PI - atan((By - Ay) / (Ax - Bx));
        if ((By < Ay) && (Bx < Ax)) arc = PI + atan((Ay - By) / (Ax - Bx));
        if ((Ax == Bx) && (By > Ay)) arc = PI / 2;
        if ((Ax == Bx) && (By < Ay)) arc = -PI / 2;
        if ((Ax < Bx) && (By == Ay)) arc = 0;
        if ((Ax > Bx) && (By == Ay)) arc = PI;
        return arc;
    }

    /**
     * Метод нахождения угла прямой по расстоянию d до точки
     * точка B - точка на прямой
     * y = k*(x-Ax) + Ay
     * точка A - точка удаленная от прямой на расстояние d
     *
     * @return
     */
    private static double[] calcArc(double Ax, double Ay, double Bx, double By, double d) {
        double a = pow(Bx - Ax, 2) - pow(d, 2);
        double b = 2 * (Bx - Ax) * (By - Ay);
        double c = pow(By - Ay, 2) - pow(d, 2);

        double D = pow(b, 2) - 4 * a * c;
        double arc1 = atan((-b + sqrt(D)) / (2 * a));
        double arc2 = atan((-b - sqrt(D)) / (2 * a));
        return new double[]{arc1, arc2};
    }

    private static double distanceBetweenTwoPoints(double Ax, double Ay, double Bx, double By) {
        return (sqrt(pow(Bx - Ax, 2) + pow(By - Ay, 2)));
    }

    private static double[] moveToNextPoint(double angle, double alpha, double t, double x, double y) {
        double k = tan(toRadians(angle));
        double Nx = x + alpha * cos(k) * t;
        double Ny = y + alpha * sin(k) * t;
        return new double[]{Nx, Ny};
    }

    public static void main(String[] args) {
        double Ax = 1;
        double Ay = 2;
        double Bx = 5;
        double By = -2;

        System.out.println("Finding Angle between two points");
        double psiOne = findArc(Ax, Ay, Bx, By);
        System.out.println("Angle in radians " + psiOne);
        System.out.println("Angle in degrees " + toDegrees(psiOne));
        System.out.println("-----------------------------------------");

        System.out.println("Finding Angle by distance from point to line");
        double Cx = 3;
        double Cy = 4;
        double d = distanceBetweenTwoPoints(Ax, Ay, Cx, Cy);
        double psiTwo = calcArc(Cx, Cy, 3, 0, d)[0];
        System.out.println("Distance Between A and C " + d);
        System.out.println("Angle in radians " + psiTwo);
        System.out.println("Angle in degrees " + toDegrees(psiTwo));
        System.out.println("-----------------------------------------");
        System.out.println("psiOne + psiTwo = " + (psiOne + abs(psiTwo)));
        System.out.println("2*PI = " + 2 * PI);
        System.out.println(-107.01241856962162 + 107.09403501783302);
        System.out.println(-28.830720133443833+28.85576324282347);

    }
}
