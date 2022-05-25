package sample.others;

public class MyMath {

    public static double lineToOne(double a, double b, double x) {
        return (a * x + b);
    }

    public static double lineToTwo(double a1, double b1, double x) {
        return (a1 * x + b1);
    }

    public static double[] calcPointZ(double a, double b, double a1, double b1) {
        double Zx = -(b - b1) / (a - a1);
        double Zy = a * Zx + b;
        return new double[]{Zx, Zy};
    }

    public static void main(String[] args) {
        double Ax = -2;
        double Ay = -1;
        double Bx = 2;
        double By = 3;

        double a1 = Ay - (Ax * (By - Ay)) / (Bx - Ax);
        double b1 = (By - Ay) / (Bx - Ax);

        double a = -1;
        double b = -2;

        System.out.println("Zx = " + calcPointZ(a, b, a1, b1)[0]);
        System.out.println("Zy = " + calcPointZ(a, b, a1, b1)[1]);

        double Zx = calcPointZ(a, b, a1, b1)[0];
        double Zy = calcPointZ(a, b, a1, b1)[1];

        System.out.println(lineToOne(a, b, Zx));
     /*   System.out.println("a1 " + a1);
        System.out.println("b1 " + b1);*/
        System.out.println(lineToTwo(a1, b1, Zx));
    }
}
