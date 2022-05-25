import static java.lang.Math.*;


public class FindNBzero {
    private static double Cx, Cy; // initial position of pursuers
    private static int k;

    /**
     * Метод должен вернуть угол phi и координаты точек D N
     */
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

    /* Method which finds coordinates of point on circle*/
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

    public static double[] secondPursuer(double r, double l, double BxOne, double ByOne) {
        double phiTwo = calcArc(BxOne, ByOne, Cx, Cy, l)[1]; // InRad
        double kTwo = tan(phiTwo);
        double SxTwo = findPointOnCircleLow(r, l, phiTwo)[0];
        double SyTwo = findPointOnCircleLow(r, l, phiTwo)[1];
        System.out.println("STwo: " + "(" + SxTwo + "; " + SyTwo + ")");

        double AxTwo = findPointsOnCircleAndLine(kTwo, SxTwo, SyTwo, r)[0][0];
        double AyTwo = findPointsOnCircleAndLine(kTwo, SxTwo, SyTwo, r)[0][1];
        double BxTwo = findPointsOnCircleAndLine(kTwo, SxTwo, SyTwo, r)[1][0];
        double ByTwo = findPointsOnCircleAndLine(kTwo, SxTwo, SyTwo, r)[1][1];
        System.out.println("phiOne: " + phiTwo); // InRad
        System.out.println("AxTwo: " + "(" + AxTwo + "; " + AyTwo + ")");
        System.out.println("BxTwo: " + "(" + BxTwo + "; " + ByTwo + ")");
        return new double[]{phiTwo, AxTwo, AyTwo, BxTwo, ByTwo};
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


    public static double findMN(double phi, double r) {
        double k = tan(phi);
        double Mx = findPointsOnCircleAndLine(k, 0, 0, r)[0][0];
        double My = findPointsOnCircleAndLine(k, 0, 0, r)[0][1];
        double Nx = findPointsOnCircleAndLine(k, 0, 0, r)[1][0];
        double Ny = findPointsOnCircleAndLine(k, 0, 0, r)[1][1];
        double MN = distanceBetweenTwoPoints(Mx, My, Nx, Ny);
        return MN;
    }

    public static double[] pointOnIntersectionTwoLines(double k1, double Bx, double By, double k2, double Mx, double My) {
        double x = (My - k2 * Mx + k1 * Bx - By) / (k1 - k2);
        double y = k1 * (x - Bx) + By;
        return new double[]{x, y};
    }

    /*Calculates distance between two points*/
    private static double distanceBetweenTwoPoints(double Ax, double Ay, double Bx, double By) {
        return (sqrt(pow(Bx - Ax, 2) + pow(By - Ay, 2)));
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

    public static void main(String[] args) {
        double r = 6;
        double l = 3;
        Cx = -(r + l);
        Cy = 0;

        k = 1;
        // Один игрок ловит
        if (r <= l) {
            System.out.println("k = " + k);
            double phiOne = 0;
            System.out.println("phi:" + k + " = " + phiOne); // InRad
        } else {
            // Информация о первом игроке
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
            double phiNext = calcArc(Bx, By, Cx, Cy, l)[1];
            k++;
            System.out.println("phi" + k + "  = " + toDegrees(phiNext));
            // Информация о втором игроке
            double SxNext = findPointOnCircleLow(r, l, phiNext)[0];
            double SyNext = findPointOnCircleLow(r, l, phiNext)[1];
            double AxNext = findPointsOnCircleAndLine(tan(phiNext), SxNext, SyNext, r)[1][0];
            double AyNext = findPointsOnCircleAndLine(tan(phiNext), SxNext, SyNext, r)[1][1];
            double BxNext = findPointsOnCircleAndLine(tan(phiNext), SxNext, SyNext, r)[0][0];
            double ByNext = findPointsOnCircleAndLine(tan(phiNext), SxNext, SyNext, r)[0][1];
            double UxNext = findPointOnCircleUp(r, l, phiNext)[0];
            double UyNext = findPointOnCircleUp(r, l, phiNext)[1];

            System.out.println("S" + k + "(" + SxNext + "; " + SyNext + ")");
            System.out.println("A" + k + "(" + AxNext + "; " + AyNext + ")");
            System.out.println("B" + k + "(" + BxNext + "; " + ByNext + ")");
            phiNext = calcArc(BxNext, ByNext, Cx, Cy, l)[1];
            System.out.println("phi" + k + " = " + toDegrees(phiNext));
            System.out.println("A" + k + "(" + AxNext + "; " + AyNext + ")");
            System.out.println("B" + k + "(" + BxNext + "; " + ByNext + ")");
            phiNext = calcArc(BxNext, ByNext, Cx, Cy, l)[1];
            System.out.println("phi" + k + " = " + toDegrees(phiNext));
            while (true) {
                // Информация о игроках следующий после второго
                SxNext = findPointOnCircleLow(r, l, phiNext)[0];
                SyNext = findPointOnCircleLow(r, l, phiNext)[1];
                AxNext = findPointsOnCircleAndLine(tan(phiNext), SxNext, SyNext, r)[1][0];
                AyNext = findPointsOnCircleAndLine(tan(phiNext), SxNext, SyNext, r)[1][1];
                BxNext = findPointsOnCircleAndLine(tan(phiNext), SxNext, SyNext, r)[0][0];
                ByNext = findPointsOnCircleAndLine(tan(phiNext), SxNext, SyNext, r)[0][1];
                UxNext = findPointOnCircleUp(r, l, phiNext)[0];
                UyNext = findPointOnCircleUp(r, l, phiNext)[1];
                k++;

                if (Double.isNaN(AxNext) || Double.isNaN(AyNext)) {
                    phiNext = findArc(UxNext,UyNext,BxNext,ByNext);
                    System.out.println("phi" + k + " = " + toDegrees(phiNext));
                    k++;
                    break;
                } else {

                }
                System.out.println("A" + k + "(" + AxNext + "; " + AyNext + ")");
                System.out.println("B" + k + "(" + BxNext + "; " + ByNext + ")");
                phiNext = calcArc(BxNext, ByNext, Cx, Cy, l)[1];
                System.out.println("phi" + k + " = " + toDegrees(phiNext));
            }
        }
        System.out.println(k);
    }
}

//            double Hx = pointOnIntersectionTwoLines(tan(phiTwo), Bx, By, -1 / tan(phiTwo), 0, 0)[0];
//            double Hy = pointOnIntersectionTwoLines(tan(phiTwo), Bx, By, -1 / tan(phiTwo), 0, 0)[1];
//            System.out.println("H" + "(" + Hx + "; " + Hy + ")");
//            double Mx =  findPointsOnCircleAndLine(tan(phiTwo), 0, 0, r)[0][0];
//            double My =  findPointsOnCircleAndLine(tan(phiTwo), 0, 0, r)[0][1];
//            double ZB = distanceBetweenTwoPoints(Mx, My, Bx, By);

//     if (Double.isNaN(AxOne)||Double.isNaN(AyOne) ) {
//            System.out.println("k ");
//        }
//        else{
//            System.out.println("AOne: " + "(" + AxOne + "; " + AyOne + ")");
//            System.out.println("BOne: " + "(" + BxOne + "; " + ByOne + ")");
//            double phiPrev = phiOne;
//            double phiNext;
//            double Ax = AxOne;
//            double Ay = AyOne;
//            double Bx = BxOne;
//            double By = ByOne;
//            do{
//                Ax = secondPursuer(r, l, BxOne, ByOne)[1];
//                Ay = secondPursuer(r, l, BxOne, ByOne)[2];
//                Bx = secondPursuer(r, l, BxOne, ByOne)[3];
//                By = secondPursuer(r, l, BxOne, ByOne)[4];
//            }
//            while(Double.isNaN(Ax) || Double.isNaN(Ay));
//        }
