package sample;

import com.sun.nio.zipfs.ZipFileAttributes;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.fxml.Initializable;
import javafx.scene.Node;
import javafx.scene.control.Label;
import javafx.scene.control.TextField;
import javafx.stage.Stage;
import sample.pursuer.Pursuer;
import sample.pursuer.Segment;

import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.ResourceBundle;

import static java.lang.Math.*;
import static sample.Controller.*;

public class StrategyDialogController implements Initializable {
    @FXML
    TextField numberOfStrategy;
    @FXML
    Label labelForT1;
    @FXML
    Label labelForT2;
    @FXML
    Label labelForT3;
    @FXML
    TextField valueT1;
    @FXML
    TextField valueT2;
    @FXML
    TextField valueT3;
    @FXML
    Label labelOfStrategy;
    @FXML
    Label labelForMessage;

    double r, l, alpha, beta;
    int s, k;
    static final double eps = 0.02;

    private ObservableList<Pursuer> appMaintvPursuers = FXCollections.observableArrayList();
    private ArrayList<Pursuer> appMainListPursuers = new ArrayList<>();

    private Pursuer pursuerOne, pursuerTwo, pursuerThree;
    private double AxOne, AyOne, BxOne, ByOne; // points after first pursuer
    private double AxTwo, AyTwo, BxTwo, ByTwo; // points after second pursuer
    private double Cx, Cy; // initial position of pursuers
    private double kphiOne, kpsiOne, kphiTwo, kpsiTwo, kpsiThree;

    public void setparameters(int vs, int vk, double vr, double vl, double valpha, double vbeta) {
        s = vs;
        k = vk;
        r = vr;
        l = vl;
        alpha = valpha;
        beta = vbeta;
        switch (s) {
            case 1: {
                labelOfStrategy.setText("First Strategy");
                break;
            }
            case 2: {
                labelOfStrategy.setText("Second Strategy");
                break;
            }
        }
        switch (k) {
            case 2: {
                labelForT2.setVisible(true);
                valueT2.setVisible(true);
                break;
            }
            case 3: {
                labelForT2.setVisible(true);
                valueT2.setVisible(true);
                labelForT3.setVisible(true);
                valueT3.setVisible(true);
                break;
            }
        }
    }

    /* Method which finds coordinates of point on circle*/
    private static double[][] findPointsOnCircle(double Sx, double Sy, double radius, double angleInRad) {
        double Hx = Sx + sin(angleInRad) * radius;
        double Hy = Sy - cos(angleInRad) * radius;
        double Vx = Sx - cos(angleInRad) * radius;
        double Vy = Sy + sin(angleInRad) * radius;
        return new double[][]{{Hx, Hy}, {Vx, Vy}};
    }

    public static double[] calcPointZ(double a, double b, double a1, double b1) {
        double Zx = (b - b1) / (a1 - a);
        double Zy = a * Zx + b;
        return new double[]{Zx, Zy};
    }

    /* Method for finding coordinates of point on circle and line*/
    // point center of circle in this method (0,0)
    private static double[][] findPointsOnCircleAndLine(double k, double Ax, double Ay, double r) {
        double a = 1 + pow(k, 2);
        double b = 2 * k * (Ay - k * Ax);
        double c = pow(k * Ax - Ay, 2) - pow(r, 2);
        double d = pow(b, 2) - 4 * a * c;
        if (d < 0) {
            /*      System.out.println("Прямая и точка не пересекаются");*/
            return new double[][]{{0, 0}, {0, 0}};
        }
        double x1 = (-b + sqrt(d)) / (2 * a);
        double y1 = k * (x1 - Ax) + Ay;
        double x2 = (-b - sqrt(d)) / (2 * a);
        double y2 = k * (x2 - Ax) + Ay;
        return new double[][]{{x1, y1}, {x2, y2}};
    }

    private static double[][] findPointsOnCircleAndLineTwo(double k, double Sx, double Sy, double r) {
        double a = 1 + pow(k, 2);
        double b = 2 * k * Sy - 2 * pow(k, 2) * Sx;
        double c = pow(k * Sx, 2) - 2 * k * Sy * Sx + pow(Sy, 2) - pow(r, 2);
        double d = pow(b, 2) - 4 * a * c;
        if (d < 0) {
            System.out.println("Прямая и точка не пересекаются");
            return new double[][]{{0, 0}, {0, 0}};
        }
        double x1 = (-b + sqrt(d)) / (2 * a);
        double y1 = k * (x1 - Sx) + Sy;
        double x2 = (-b - sqrt(d)) / (2 * a);
        double y2 = k * (x2 - Sx) + Sy;
        return new double[][]{{x1, y1}, {x2, y2}};
    }

    /**
     * Method that calculates k
     *
     * @param Ax - coordinates of the point on the axis ox from line equals d
     * @param Ay - coordinates of the point on the axis oy from line equals d
     * @param Bx - coordinates of the point on the axis ox on the line
     * @param By - coordinates of the point on the axis oy on the line
     * @param d  - distance from point to line
     * @return k - line slope
     */
    private static double[] distanceFromPointToLine(double Ax, double Ay, double Bx, double By, double d) {
        double a = pow(Ax - Bx, 2) - pow(d, 2);
        double b = 2 * (Ax - Bx) * (By - Ay);
        double c = pow(Ay - By, 2) - pow(d, 2);
        double D = pow(b, 2) - 4 * a * c;
        double k1 = (-b + sqrt(D)) / (2 * a);
        double k2 = (-b - sqrt(D)) / (2 * a);
        return new double[]{k1, k2};
    }

    private static double[] distanceFromPointToLineTwo(double Bx, double By, double Cx, double Cy, double l) {
        double a = pow(l, 2) - pow((Bx - Cx), 2);
        double b = 2 * (Bx - Cx) * (By - Cy);
        double c = pow(l, 2) - pow((Cy - By), 2);
        double D = pow(b, 2) - 4 * a * c;
        double k1 = (-b + sqrt(D)) / (2 * a);
        double k2 = (-b - sqrt(D)) / (2 * a);
        return new double[]{k1, k2};
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

    /**
     * Метод нахождения углового коэффицента прямой по расстоянию d до точки O
     * точка B - точка на прямой
     * y = k*(x-Bx) + By
     *
     * @return угол в радианах
     */

    private static double[] calcArc(double Bx, double By, double Ox, double Oy, double d) {
        double a = pow(Ox - Bx, 2) - pow(d, 2);
        double b = 2 * (Ox - Bx) * (By - Oy);
        double c = pow(By - Oy, 2) - pow(d, 2);
        double D = pow(b, 2) - 4 * a * c;
        double arc1 = atan((-b + sqrt(D)) / (2 * a));
        double arc2 = atan((-b - sqrt(D)) / (2 * a));
        return new double[]{arc1, arc2};
    }

    private static double checkArc(double arc, double Bx, double By, double Ox, double Oy) {
        double k = tan(arc);
        return (abs(k * Ox - Oy - k * Bx + By) / (sqrt(pow(k, 2) + 1)));
    }

    /* Rounding double number*/
    private double roundingNumber(double x) {
        x = x * 100;
        int i = (int) round(x);
        return ((double) i / 100);
    }

    @FXML
    public void calcAndCreatePursuers(ActionEvent event) {
        switch (s) {
            case 1: {
                System.out.println("First Strategy");
                calcForFirstStrategy();
                break;
            }
            case 2: {
                System.out.println("Second Strategy");
                calcForSecondStrategy();
                break;
            }
        }
        drawCircle(0, 0, r);
        closeStage(event);
    }

    public double calcAndAddFirstPursuer() {
        double t1 = Double.parseDouble(valueT1.getText());
        Cx = -2 * sqrt(r * l);
        Cy = r - l;
        double phiOne = 0;
        kphiOne = tan(0);
        /* Calculation For First Pursuer*/
        double G1x = Cx + alpha * t1 * cos(phiOne);
        double G1y = Cy + alpha * t1 * sin(phiOne);
        /* H - точка на окружности радиуса l с центром в точке C*/
        double Hx = -2 * sqrt(r * l);
        double Hy = r - 2 * l;
        /* D - точка пересечения прямой и  окружности*/
        double Dx = -2 * sqrt(l * (r - l));
        double Dy = r - 2 * l;

        System.out.println("-------------Информация о первом игроке---------------");
        System.out.println("C (" + Cx + " ; " + Cy + ")");
        System.out.println("G1 (" + G1x + " ; " + G1y + ")");
        System.out.println("H (" + Hx + " ; " + Hy + ")");
        System.out.println("D (" + Dx + " ; " + Dy + ")");

        /*kpsiOne - угловой коэффицент для обратной траектории P1*/
        double psiOne = calcArc(G1x, G1y, Dx, Dy, l)[0];
        kpsiOne = tan(psiOne);

        System.out.println("psiOne in radians = " + psiOne);
        System.out.println("psiOne in degrees = " + toDegrees(psiOne));
        System.out.println("kpsiOne = " + kpsiOne);
        System.out.println(checkArc(psiOne, G1x, G1y, Dx, Dy));

        double U1x = -2 * sqrt(l * (r - l));
        double U1y = r - l;

        double U2x = Dx + sin(psiOne) * l;
        double U2y = Dy - cos(psiOne) * l;

        System.out.println("U1 (" + U1x + " ; " + U1y + ")");
        System.out.println("U2 (" + U2x + " ; " + U2y + ")");

/*
        System.out.println("Distance between U1 and G1 " + distanceBetweenTwoPoints(U1x, U1y, G1x, G1y));
        System.out.println("Distance between U2 and G1 " + distanceBetweenTwoPoints(U2x, U2y, G1x, G1y));*/
        /* расстояния U1G1 и U2G1 должны быть равны*/

        /* Нахождение угла по точкам U2 и G1  */
     /*   double tetaOne = findArc(U2x, U2y, G1x, G1y);
        System.out.println("teta " + tetaOne);
        *//*Должна быть равна нулю*//*
        System.out.println("Погрешность " + (tetaOne - psiOne)); // = 0*/

        /*Q - точка на окружности радиуса l с центром в точке G1*/
        double Qx = findPointsOnCircle(G1x, G1y, l, psiOne)[0][0];
        double Qy = findPointsOnCircle(G1x, G1y, l, psiOne)[0][1];
        System.out.println("Q(" + Qx + "; " + Qy + ")");

        double IOnex = G1x - sin(psiOne) * l;
        double IOney = G1y + cos(psiOne) * l;
        double ITwox = G1x + sin(psiOne) * l;
        double ITwoy = G1y - cos(psiOne) * l;

        System.out.println(ITwox + " " + ITwoy);

        double Mx = findPointsOnCircleAndLine(kpsiOne, Qx, Qy, r)[0][0];
        double My = findPointsOnCircleAndLine(kpsiOne, Qx, Qy, r)[0][1];
        double Nx = findPointsOnCircleAndLine(kpsiOne, Qx, Qy, r)[1][0];
        double Ny = findPointsOnCircleAndLine(kpsiOne, Qx, Qy, r)[1][1];

        /*G1F - конечная точка первого игрока*/
        double G1fx = Nx - sin(psiOne) * l;
        double G1fy = Ny + cos(psiOne) * l;
        System.out.println("G1x = " + G1x);
        System.out.println("G1y = " + G1y);
        System.out.println("G1fx = " + G1fx);
        System.out.println("G1Fy = " + G1fy);
        System.out.println("-----------------------------------------");

        double CG = distanceBetweenTwoPoints(Cx, Cy, G1x, G1y);
        double NQ = distanceBetweenTwoPoints(Nx, Ny, Qx, Qy);

/*        drawPoint(U1x, U1y);
        drawPoint(U2x, U2y);
        drawLine(U1x, U1y, Dx, Dy);
        drawLine(U2x, U2y, Dx, Dy);
        drawPoint(G1x, G1y);
        drawPoint(Hx, Hy);
        drawPoint(Dx, Dy);
        // drawPoint(Ux, Uy);
        drawPoint(IOnex, IOney);
        drawPoint(ITwox, ITwoy);
        drawPoint(Mx, My);
        drawPoint(Nx, Ny);
        drawPoint(G1fx, G1fy);
        drawLine(Cx, Cy, G1x, G1y);
        drawLine(G1x, G1y, G1fx, G1fy);
        drawLine(Hx, Hy, G1x, Hy);
        drawLine(Hx, Hy + 2 * l, G1x, Hy + 2 * l);
        drawLine(IOnex, IOney, Dx, Dy);
        drawLine(ITwox, ITwoy, Nx, Ny);
        drawCircle(Cx, Cy, l);
        drawCircle(G1x, G1y, l);
        drawCircle(G1fx, G1fy, l);*/

  /*      double G1G1F = distanceBetweenTwoPoints(G1x, G1y, G1fx, G1fy);
        if (G1G1F == NQ) System.out.println("distance is equals");*/
        double t1r = (NQ / alpha);
        double T1 = (NQ + CG) / alpha;

        /* Add First Pursuer in ArrayList*/
        ArrayList<Segment> sForFirst = new ArrayList<>();
        sForFirst.add(new Segment(1, t1, 0));
        sForFirst.add(new Segment(2, t1r, 180 + toDegrees(psiOne)));

        pursuerOne = new Pursuer(1, Cx, Cy, l, alpha, sForFirst);
        appMainListPursuers.add(pursuerOne);
        appMaintvPursuers.add(pursuerOne);
        System.out.println(pursuerOne);
        AxOne = Nx;
        AyOne = Ny;
        BxOne = Mx;
        ByOne = My;
        return T1;
    }

    public double calcAndAddSecondPursuer() {
        double t2 = Double.parseDouble(valueT2.getText());
        double phiTwo = calcArc(BxOne, ByOne, Cx, Cy, l)[1];
        kphiTwo = tan(phiTwo);

        System.out.println("-----------------Информация о втором игроке-------------");
        System.out.println("phiTwo in radians = " + phiTwo);
        System.out.println("psiTwo in degrees = " + toDegrees(phiTwo));
        System.out.println("kphiTwo = " + kphiTwo);
        System.out.println(checkArc(phiTwo, BxOne, ByOne, Cx, Cy));

        double Vx = BxOne - sin(abs(phiTwo)) * l;
        double Vy = ByOne - cos(abs(phiTwo)) * l;

        /*  Точка поворота G2*/
        double G2x = Cx + t2 * cos(phiTwo);
        double G2y = Cy + t2 * sin(phiTwo);
       /* System.out.println("V (" + Vx + " ; " + Vy + ")");
        System.out.println("G2 (" + G2x + " ; " + G2y + ")");*/


    /*    double tetaTwo = 2 * PI - findArc(Cx, Cy, Vx, Vy);
        System.out.println("tetaTwo " + tetaTwo);
        *//*  Должна быть равна нулю*//*
        System.out.println("Погрешность " + (tetaTwo + phiTwo));
*/
        double SxOne = Cx - sin(phiTwo) * l;
        double SyOne = Cy + cos(phiTwo) * l;
        double SxTwo = Cx + sin(phiTwo) * l;
        double SyTwo = Cy - cos(phiTwo) * l;

        double JxOne = G2x - sin(phiTwo) * l;
        double JyOne = G2y + cos(phiTwo) * l;
        double JxTwo = G2x + sin(phiTwo) * l;
        double JyTwo = G2y - cos(phiTwo) * l;

        drawCircle(Cx, Cy, l);

        drawPoint(Vx, Vy);
        /*    drawLine(Vx, Vy, BxOne, ByOne);*/
        drawCircle(G2x, G2y, l);
        drawPoint(SxOne, SyOne);
        drawPoint(SxTwo, SyTwo);
        drawPoint(JxOne, JyOne);
        drawPoint(JxTwo, JyTwo);
        drawPoint(Cx, Cy);
        drawPoint(AxOne, AyOne);
        drawPoint(BxOne, ByOne);

        drawLine(JxOne, JyOne, JxTwo, JyTwo);
        drawLine(SxOne, SyOne, SxTwo, SyTwo);
        drawLine(AxOne, AyOne, BxOne, ByOne);
        drawLine(Cx, Cy, G2x, G2y);
        drawLine(SxOne, SyOne, JxOne, JyOne);
        drawLine(SxTwo, SyTwo, JxTwo, JyTwo);
     /*   System.out.println("kpsiOne = " + kpsiOne);
        double x = -180;
        drawLine(BxOne, ByOne, x, kpsiOne * (x - BxOne) + ByOne);

        System.out.println("kphiTwo = " + kphiTwo);
        double x1 = 200;
        drawLine(SxTwo, SyTwo, x1, kphiTwo * (x1 - SxTwo) + SyTwo);*/

        /*  Z - точка пересечения двух прямых
         * y = kphiTwo(x - Bx) + By
         * y = kpsiOne(x - Ax) + Ay
         * */
        double a = kpsiOne;
        double b = -BxOne * kpsiOne + ByOne;
        double a1 = kphiTwo;
        double b1 = -SxTwo * kphiTwo + SyTwo;
        double Zx = calcPointZ(a, b, a1, b1)[0];
        double Zy = calcPointZ(a, b, a1, b1)[1];
       /* System.out.println("Zx = " + Zx);
        System.out.println("Zy = " + Zy);*/
        drawPoint(Zx, Zy);

        /* Точка поворота W */
        double k1 = kpsiOne;
        double k2 = kphiTwo;

        double Wx = (k1 * Zx - k2 * BxOne - Zy + ByOne + l * (sqrt(pow(k1, 2) + 1) - sqrt(pow(k2, 2) + 1))) / (k1 - k2);
        double Wy = k2 * (Wx - Vx) + Vy;

        /*drawPoint(Wx, Wy);
        drawCircle(Wx, Wy, l);*/

        double CW = distanceBetweenTwoPoints(Cx, Cy, Wx, Wy) / alpha;
        System.out.println("t* = " + CW);

        /*  psiTwo - угол для обратной траектории P2*/
        double psiTwo, psiTwo1, psiTwo2;
        psiTwo1 = calcArc(Zx, Zy, G2x, G2y, l)[0];
        psiTwo2 = calcArc(AxOne, AyOne, G2x, G2y, l)[0];
        double t2r, T2;
        if (t2 <= CW) {
            System.out.println("По точке A");
            psiTwo = psiTwo2;
            kpsiTwo = tan(psiTwo);
            double CG2 = distanceBetweenTwoPoints(Cx, Cy, G2x, G2y);
            double AG = distanceBetweenTwoPoints(AxOne, AyOne, G2x, G2y);
            t2r = sqrt(pow(AG, 2) - pow(l, 2));
        } else {
            System.out.println("По точке Z");
            psiTwo = psiTwo1;
            kpsiTwo = tan(psiTwo);
            if (k == 2) {
                double A = -(1 / kpsiTwo);
                double B = (Zx / kpsiTwo) + Zy;
                double A1 = kpsiTwo;
                double B1 = -G2x * kpsiTwo + G2y;

                double Ix = calcPointZ(A, B, A1, B1)[0];
                double Iy = calcPointZ(A, B, A1, B1)[1];
                drawPoint(Ix, Iy);
                /*System.out.println("Проверка точки I = " + distanceBetweenTwoPoints(Ix, Iy, Zx, Zy));*/

                double F = -(1 / kpsiTwo);
                double G = (AxOne / kpsiTwo) + AyOne;
                double F1 = kpsiTwo;
                double G1 = -G2x * kpsiTwo + G2y;
                double Jx = calcPointZ(F, G, F1, G1)[0];
                double Jy = calcPointZ(F, G, F1, G1)[1];

                double AJ = distanceBetweenTwoPoints(AxOne, AyOne, Jx, Jy);
                double JV = sqrt(pow(l, 2) - pow(AJ, 2));
                drawPoint(Jx, Jy);
                System.out.println("J (" + Jx + "; " + Jy);
                System.out.println("AJ = " + AJ);
                System.out.println("JV =" + JV);

                double CG2 = distanceBetweenTwoPoints(Cx, Cy, G2x, G2y);
                double AG = distanceBetweenTwoPoints(AxOne, AyOne, G2x, G2y);
                t2r = sqrt(pow(AG, 2) - pow(l, 2)) - JV;
            } else {
                double CG2 = distanceBetweenTwoPoints(Cx, Cy, G2x, G2y);
                double AG = distanceBetweenTwoPoints(AxOne, AyOne, G2x, G2y);
                t2r = sqrt(pow(AG, 2) - pow(l, 2));
            }
        }
        T2 = (t2 + t2r) / alpha;

        System.out.println("psiTwo in radians = " + psiTwo);
        System.out.println("psiTwo in degrees = " + toDegrees(psiTwo));
        System.out.println("kpsiTwo = " + kpsiTwo);
        System.out.println(checkArc(psiTwo, Zx, Zy, G2x, G2y));

        double G2fx = G2x - t2r * cos(psiTwo);
        double G2fy = G2y - t2r * sin(psiTwo);
        double QxOne = G2x - sin(psiTwo) * l;
        double QyOne = G2y + cos(psiTwo) * l;
        double QxTwo = G2x + sin(psiTwo) * l;
        double QyTwo = G2y - cos(psiTwo) * l;

        double Mx = findPointsOnCircleAndLine(kpsiTwo, QxTwo, QyTwo, r)[0][0];
        double My = findPointsOnCircleAndLine(kpsiTwo, QxTwo, QyTwo, r)[0][1];
        double Nx = findPointsOnCircleAndLine(kpsiTwo, QxTwo, QyTwo, r)[1][0];
        double Ny = findPointsOnCircleAndLine(kpsiTwo, QxTwo, QyTwo, r)[1][1];

        AxTwo = Nx;
        AyTwo = Ny;
        BxTwo = Mx;
        ByTwo = My;

        drawLine(G2x, G2y, G2fx, G2fy);
        drawLine(QxOne, QyOne, QxOne - t2r * cos(psiTwo), QyOne - t2r * sin(psiTwo));
        drawLine(QxTwo, QyTwo, QxTwo - t2r * cos(psiTwo), QyTwo - t2r * sin(psiTwo));
        drawPoint(QxOne, QyOne);
        drawPoint(QxTwo, QyTwo);
        drawCircle(G2fx, G2fy, l);


        ArrayList<Segment> sForSecond = new ArrayList<>();
        sForSecond.add(new Segment(1, t2, toDegrees(phiTwo)));
        sForSecond.add(new Segment(2, t2r, toDegrees(PI + psiTwo)));
        pursuerTwo = new Pursuer(2, Cx, Cy, l, alpha, sForSecond);
        System.out.println(pursuerTwo);
        appMainListPursuers.add(pursuerTwo);
        appMaintvPursuers.add(pursuerTwo);
        return T2;
    }

    public double calcAndAddThirdPursuer() {
        double t3 = Double.parseDouble(valueT3.getText());
        double phiThree = calcArc(BxTwo, ByTwo, Cx, Cy, l)[1];
        double kphiThree = tan(phiThree);

        System.out.println("-----------------Информация о третьем игроке-------------");
        System.out.println("phiThree in radians = " + phiThree);
        System.out.println("psiThree in degrees = " + toDegrees(phiThree));
        System.out.println("kphiThree = " + kphiThree);
        System.out.println(checkArc(phiThree, BxTwo, ByTwo, Cx, Cy));

        double Vx = BxTwo - sin(abs(phiThree)) * l;
        double Vy = ByTwo - cos(abs(phiThree)) * l;
        drawPoint(Vx, Vy);

        double G3x = Cx + t3 * cos(phiThree);
        double G3y = Cy + t3 * sin(phiThree);
        System.out.println("G3(" + G3x + "; " + G3y + ")");

        double SxOne = Cx - sin(phiThree) * l;
        double SyOne = Cy + cos(phiThree) * l;
        double SxTwo = Cx + sin(phiThree) * l;
        double SyTwo = Cy - cos(phiThree) * l;

        double JxOne = G3x - sin(phiThree) * l;
        double JyOne = G3y + cos(phiThree) * l;
        double JxTwo = G3x + sin(phiThree) * l;
        double JyTwo = G3y - cos(phiThree) * l;

        drawCircle(Cx, Cy, l);
        drawCircle(G3x, G3y, l);
        drawPoint(SxOne, SyOne);
        drawPoint(SxTwo, SyTwo);
      /*  drawPoint(JxOne, JyOne);
        drawPoint(JxTwo, JyTwo);*/
        drawPoint(AxTwo, AyTwo);
        drawPoint(BxTwo, ByTwo);
//        drawLine(AxOne, AyOne, BxOne, ByOne);
        drawLine(AxTwo, AyTwo, BxTwo, ByTwo);
        drawLine(Cx, Cy, G3x, G3y);
        drawLine(SxOne, SyOne, JxOne, JyOne);
        drawLine(SxTwo, SyTwo, JxTwo, JyTwo);


        double a = kpsiTwo;
        double b = -BxTwo * kpsiTwo + ByTwo;
        double a1 = kphiThree;
        double b1 = -SxTwo * kphiThree + SyTwo;
        double Zx = calcPointZ(a, b, a1, b1)[0];
        double Zy = calcPointZ(a, b, a1, b1)[1];
        System.out.println("Zx = " + Zx);
        System.out.println("Zy = " + Zy);
        drawPoint(Zx, Zy);

        double k1 = kpsiTwo;
        double k2 = kphiThree;

        double Wx = (k1 * Zx - k2 * BxTwo - Zy + ByTwo + l * (sqrt(pow(k1, 2) + 1) - sqrt(pow(k2, 2) + 1))) / (k1 - k2);
        double Wy = k2 * (Wx - Vx) + Vy;
        drawPoint(Wx, Wy);
        /*drawCircle(Wx, Wy, l);*/
        double CW = distanceBetweenTwoPoints(Cx, Cy, Wx, Wy) / alpha;
        System.out.println("t* = " + CW);

        double psiThree, psiThree1, psiThree2;
        psiThree1 = calcArc(Zx, Zy, G3x, G3y, l)[0];
        psiThree2 = calcArc(AxTwo, AyTwo, G3x, G3y, l)[0];

        double t3r, T3;
        if (t3 <= CW) {
            System.out.println("По точке A");
            psiThree = psiThree2;
            kpsiThree = tan(psiThree);
            double CG2 = distanceBetweenTwoPoints(Cx, Cy, G3x, G3y);
            double AG = distanceBetweenTwoPoints(AxTwo, AyTwo, G3x, G3y);
            t3r = sqrt(pow(AG, 2) - pow(l, 2));

        } else {
            System.out.println("По точке Z");
            psiThree = psiThree1;
            kpsiThree = tan(psiThree);

            double A = -(1 / kpsiThree);
            double B = (Zx / kpsiThree) + Zy;
            double A1 = kpsiThree;
            double B1 = -G3x * kpsiThree + G3y;
            double Ix = calcPointZ(A, B, A1, B1)[0];
            double Iy = calcPointZ(A, B, A1, B1)[1];
            drawPoint(Ix, Iy);

            System.out.println("Проверка точки I = " + distanceBetweenTwoPoints(Zx, Zy, Ix, Iy));

            double F = -(1 / kpsiThree);
            double G = (AxTwo / kpsiThree) + AyTwo;
            double F1 = kpsiThree;
            double G1 = -G3x * kpsiThree + G3y;
            double Jx = calcPointZ(F, G, F1, G1)[0];
            double Jy = calcPointZ(F, G, F1, G1)[1];


            double AJ = distanceBetweenTwoPoints(AxTwo, AyTwo, Jx, Jy);
            double JV = sqrt(pow(l, 2) - pow(AJ, 2));
            drawPoint(Jx, Jy);

            System.out.println("J (" + Jx + "; " + Jy);
            System.out.println("AJ = " + AJ);
            System.out.println("JV =" + JV);

            double CG2 = distanceBetweenTwoPoints(Cx, Cy, G3x, G3y);
            double AG = distanceBetweenTwoPoints(AxTwo, AyTwo, G3x, G3y);
            t3r = sqrt(pow(AG, 2) - pow(l, 2)) - JV + 0.2;
        }
        drawCircle(G3x - t3r * cos(psiThree), G3y - t3r * sin(psiThree), l);

        System.out.println("psiThree1 = " + psiThree1);
        System.out.println("psiThree2 = " + psiThree2);

        kpsiThree = tan(psiThree);


        double QxOne = G3x - sin(psiThree) * l;
        double QyOne = G3y + cos(psiThree) * l;
        double QxTwo = G3x + sin(psiThree) * l;
        double QyTwo = G3y - cos(psiThree) * l;


       /* drawPoint(QxOne, QyOne);
        drawPoint(QxTwo, QyTwo);
*/

        T3 = (t3 + t3r) / alpha;
        ArrayList<Segment> sForThird = new ArrayList<>();
        sForThird.add(new Segment(1, t3, toDegrees(phiThree)));
        sForThird.add(new Segment(1, t3r, toDegrees(PI + psiThree)));
        pursuerThree = new Pursuer(3, Cx, Cy, l, alpha, sForThird);
        appMainListPursuers.add(pursuerThree);
        appMaintvPursuers.add(pursuerThree);
        return T3;
    }

    public void calcForFirstStrategy() {
        double T1 = calcAndAddFirstPursuer();
        System.out.println("AOne: " + AxOne + " " + AyOne);
        System.out.println("BOne: " + BxOne + " " + ByOne);
        System.out.println("T1 = " + T1);
        if (k == 2) {
            double T2 = calcAndAddSecondPursuer();
            System.out.println("ATwo: " + AxTwo + " " + AyTwo);
            System.out.println("BTwo: " + BxTwo + " " + ByTwo);
            System.out.println("T2 = " + T2);
            System.out.println("Max Time of two Pursuers = " + max(T1, T2));
        }
        if (k == 3) {
            double T2 = calcAndAddSecondPursuer();
            System.out.println("T2 = " + T2);
            double T3 = calcAndAddThirdPursuer();
            System.out.println("T3 = " + T3);
            System.out.println("Max Time of two Pursuers = " + max(T1, max(T2, T3)));
        }

    }

    public void calcForSecondStrategy() {

        switch (k) {
            case 2: {
                System.out.println("Two Pursuers");
                break;
            }
            case 3: {
                System.out.println("Three Pursuers");
                break;
            }
        }
    }

    private void closeStage(ActionEvent event) {
        Node source = (Node) event.getSource();
        Stage stage = (Stage) source.getScene().getWindow();
        stage.close();
    }

    public void setObservableList(ObservableList<Pursuer> tvPursuers) {
        this.appMaintvPursuers = tvPursuers;
    }

    public void setPursuerList(ArrayList<Pursuer> pursuersList) {
        this.appMainListPursuers = pursuersList;
    }

    @Override
    public void initialize(URL location, ResourceBundle resources) {
        valueT1.setText("450");
        valueT2.setText("420");
        valueT3.setText("460");
    }
}
