import static java.lang.Math.*;

public class FindCountPursuersBetaZero {

    // Beta = 0

    public static double[] firstPursuer(double r, double l) {
        double phi = asin((r - l) / (r + l)); // phi в радианах
        double phiInDeg = Math.toDegrees(phi); // Угол в градусах
        double xs = -(r + l) + sin(phi) * l;
        double ys = -cos(phi) * l;
        double k = tan(phi);
        double a = (1 + pow(k, 2));
        double b = 2 * (pow(k, 2) * (r + l + sin(phi) * l) - k * cos(phi) * l);
        double c = pow(k, 2) * pow((r + l + sin(phi) * l), 2) + pow(cos(phi) * l, 2) - 2 * k * cos(phi) * l * (r + l + sin(phi) * l) - pow(r, 2);
        double Bx = (-b + sqrt(pow(b, 2) - 4 * a * c)) / (2 * a);
        double By = k * (Bx + r + l + sin(phi) * l) - cos(phi) * l;
        double Ax = (-b - sqrt(pow(b, 2) - 4 * a * c)) / (2 * a);
        double Ay = k * (Ax + r + l + sin(phi) * l) - cos(phi) * l;
        return new double[]{phi, Ax, Ay, Bx, By};
    }

    private static double[] findPointOnCircleLow(double r, double l, double angleInRad) {
        double Sx = 0;
        double Sy = 0;
        Sx = -(r + l) + sin(angleInRad) * l;
        Sy = -cos(angleInRad) * l;
        return new double[]{Sx, Sy};
    }

    private static double[] findPointOnCircleUp(double r, double l, double angleInRad) {
        double Sx = 0;
        double Sy = 0;
        Sx = -(r + l) - sin(angleInRad) * l;
        Sy = cos(angleInRad) * l;
        return new double[]{Sx, Sy};
    }

    private static double[][] findPointsOnCircleAndLine(double k, double Ax, double Ay, double r) {
        double a = 1 + pow(k, 2);
        double b = 2 * k * (Ay - k * Ax);
        double c = pow(k * Ax - Ay, 2) - pow(r, 2);
        double d = pow(b, 2) - 4 * a * c;
        if (d < 0) {
            return new double[][]{{Double.NaN, Double.NaN}, {Double.NaN, Double.NaN}};
        }
        double x1 = (-b + sqrt(d)) / (2 * a);
        double y1 = k * (x1 - Ax) + Ay;
        double x2 = (-b - sqrt(d)) / (2 * a);
        double y2 = k * (x2 - Ax) + Ay;
        return new double[][]{{x1, y1}, {x2, y2}};
    }

    private static double[] calcArc(double Bx, double By, double Ox, double Oy, double d) {
        double a = pow(Ox - Bx, 2) - pow(d, 2);
        double b = 2 * (Ox - Bx) * (By - Oy);
        double c = pow(By - Oy, 2) - pow(d, 2);
        double D = pow(b, 2) - 4 * a * c;
        double arc1 = atan((-b + sqrt(D)) / (2 * a));
        double arc2 = atan((-b - sqrt(D)) / (2 * a));
        return new double[]{arc1, arc2};
    }

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

    public static int findCountPursuersBZ(double r, double l) {
        double Cx = -(r + l);
        double Cy = 0;
        if (r <= l) {
            System.out.println("phiOne = 0");
            return 1;
        } else {
            int k = 1;

            // Информация о первом игроке
            System.out.println("--------------------------------------");
            System.out.println("Информация о первом игроке");
            double phiOne = firstPursuer(r, l)[0];
            System.out.println("phi" + k + "  = " + toDegrees(phiOne)); // InRad
            // Нахождение точек на окружности
            double Sx = findPointOnCircleLow(r, l, phiOne)[0];
            double Sy = findPointOnCircleLow(r, l, phiOne)[1];
            System.out.println("S" + k + "(" + Sx + "; " + Sy + ")");
            double Ux = findPointOnCircleUp(r, l, phiOne)[0];
            double Uy = findPointOnCircleUp(r, l, phiOne)[1];
            System.out.println("U" + k + "(" + Ux + "; " + Uy + ")");
            // Нахождение двух точке пересечения нижней прямой с окружностью
            double Ax = findPointsOnCircleAndLine(tan(phiOne), Sx, Sy, r)[1][0];
            double Ay = findPointsOnCircleAndLine(tan(phiOne), Sx, Sy, r)[1][1];
            double Bx = findPointsOnCircleAndLine(tan(phiOne), Sx, Sy, r)[0][0];
            double By = findPointsOnCircleAndLine(tan(phiOne), Sx, Sy, r)[0][1];
            System.out.println("A" + k + "(" + Ax + "; " + Ay + ")");
            System.out.println("B" + k + "(" + Bx + "; " + By + ")");

            // Информация о втором игроке
            System.out.println("--------------------------------------");
            System.out.println("Информация о втором игроке");
            k++;
            double phiNext = calcArc(Bx, By, Cx, Cy, l)[1];
            System.out.println("phi" + k + "  = " + toDegrees(phiNext));
            double SxNext = findPointOnCircleLow(r, l, phiNext)[0];
            double SyNext = findPointOnCircleLow(r, l, phiNext)[1];
            double AxNext = findPointsOnCircleAndLine(tan(phiNext), SxNext, SyNext, r)[1][0];
            double AyNext = findPointsOnCircleAndLine(tan(phiNext), SxNext, SyNext, r)[1][1];
            double BxNext = findPointsOnCircleAndLine(tan(phiNext), SxNext, SyNext, r)[0][0];
            double ByNext = findPointsOnCircleAndLine(tan(phiNext), SxNext, SyNext, r)[0][1];
            double UxNext = findPointOnCircleUp(r, l, phiNext)[0];
            double UyNext = findPointOnCircleUp(r, l, phiNext)[1];
            System.out.println("U" + k + "(" + UxNext + "; " + UyNext + ")");
            System.out.println("S" + k + "(" + SxNext + "; " + SyNext + ")");
            System.out.println("A" + k + "(" + AxNext + "; " + AyNext + ")");
            System.out.println("B" + k + "(" + BxNext + "; " + ByNext + ")");

            while (true) {
                k++;
                System.out.println("--------------------------------------");
                System.out.println("Информация о "+ k + " игроке");
                phiNext = calcArc(BxNext, ByNext, Cx, Cy, l)[1];
                System.out.println("phi" + k + " = " + toDegrees(phiNext));
                SxNext = findPointOnCircleLow(r, l, phiNext)[0];
                SyNext = findPointOnCircleLow(r, l, phiNext)[1];
                AxNext = findPointsOnCircleAndLine(tan(phiNext), SxNext, SyNext, r)[1][0];
                AyNext = findPointsOnCircleAndLine(tan(phiNext), SxNext, SyNext, r)[1][1];
                BxNext = findPointsOnCircleAndLine(tan(phiNext), SxNext, SyNext, r)[0][0];
                ByNext = findPointsOnCircleAndLine(tan(phiNext), SxNext, SyNext, r)[0][1];
                UxNext = findPointOnCircleUp(r, l, phiNext)[0];
                UyNext = findPointOnCircleUp(r, l, phiNext)[1];
                System.out.println("U" + k + "(" + UxNext + "; " + UyNext + ")");
                System.out.println("S" + k + "(" + SxNext + "; " + SyNext + ")");
                System.out.println("A" + k + "(" + AxNext + "; " + AyNext + ")");
                System.out.println("B" + k + "(" + BxNext + "; " + ByNext + ")");
                if (Double.isNaN(AxNext) || Double.isNaN(AyNext)) {
//                    Ux = findPointOnCircleUp(r, l, phiNext)[0];
//                    Uy = findPointOnCircleUp(r, l, phiNext)[1];
//                    phiNext = findArc(Ux, Uy, BxNext, ByNext);
//                    System.out.println("phi" + k + " = " + toDegrees(phiNext));

                    break;
                }
            }
            return k;
        }

    }

    public static String numberToString(int n){
        String [] strings = {"первом", "втором", "третьем","четвертом","пятом","шестом","седьмом","восьмом","девятом", "десятом"};
        return strings[n-1];
    }
    public static void main(String[] args) {
        System.out.println("k = " + findCountPursuersBZ(350, 40));
    }
}
