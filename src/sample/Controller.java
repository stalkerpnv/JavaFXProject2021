package sample;

import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.fxml.FXMLLoader;
import javafx.fxml.Initializable;
import javafx.scene.Parent;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.control.*;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.paint.Color;
import javafx.stage.Modality;
import javafx.stage.Stage;
import sample.pursuer.Pursuer;
import sample.pursuer.Segment;

import java.io.IOException;
import java.net.URL;
import java.util.*;

import static java.lang.Math.*;

public class Controller implements Initializable {
    @FXML
    public Canvas canvas;
    @FXML
    public TextField valuel;
    @FXML
    public TextField valuer;
    @FXML
    public TextField valuen;
    @FXML
    public TextField valuealpha;
    @FXML
    public TextField valuebeta;
    @FXML
    public TextField valuedt;
    @FXML
    private Button startButton;
    @FXML
    private Button pauseButton;
    @FXML
    private TableView<Pursuer> tableP;
    @FXML
    private TableColumn colId;
    @FXML
    private TableColumn colX;
    @FXML
    private TableColumn colY;
    @FXML
    private TableColumn colSegments;
    @FXML
    private ComboBox numberOfStrategy;
    @FXML
    private ComboBox numberOfPursuers;
    @FXML
    private CheckBox kCoverageFlag;

    Pursuer pursuer;
    Segment singleSegment;
    int step_i;
    int mark;
    private ObservableList<Pursuer> tvPursuers = FXCollections.observableArrayList();
    private ArrayList<Pursuer> pursuersList = new ArrayList<>();

//    private ObservableList<String> listOfStrategies = FXCollections.observableArrayList("1", "2");
//    private ObservableList<String> listOfPursuers = FXCollections.observableArrayList("1", "2", "3");

    public static GraphicsContext gc;
    private static double r, l, alpha, beta, dt, Cx, Cy;
    double tStep, tStepOne, tStepTwo;
    private static double t1, t2, t3;
    private static double T1, T2, T3;
    private static int n;
    private static double m;
    static int w, h;
    public static int Ox, Oy;
    private static int[][] infSet;
    private static int[][] infSetOne;
    private static int[][] infSetTwo;
    private static int[][] infSetIncOne;
    private static int[][] infSetIncTwo;
    static int[] stepArr;

    static Timer timer;
    static boolean timerState;
    static TimerTask timerTask;
    private static final double eps = 0.001;

    /**
     * Fill the infSet - array n*n
     *
     * @param r
     * @param n
     */
    private static void fillInfSet(double r, int n, double m) {
        int halfN = n / 2;
        infSet = new int[n][n];
        infSetOne = new int[n][n];
        infSetTwo = new int[n][n];

        // Основная область
        for (int i = 0; i < infSet.length; i++) {
            for (int j = 0; j < infSet.length; j++) {
                if (pow(r / m, 2) >= (pow((i - halfN), 2) + pow((j - halfN), 2))) {
                    infSet[i][j] = 1;
                }
            }
        }
        // Первая вспомогательная область
        infSetOne = new int[n][n];
        for (int i = 0; i < infSetOne.length; i++) {
            for (int j = 0; j < infSetOne.length; j++) {
                if (pow(r / m, 2) >= (pow((i - halfN), 2) + pow((j - halfN), 2))) {
                    infSetOne[i][j] = 1;
                }
            }
        }
        // Вторая вспомогательная область
        infSetTwo = new int[n][n];
        for (int i = 0; i < infSetTwo.length; i++) {
            for (int j = 0; j < infSetTwo.length; j++) {
                if (pow(r / m, 2) >= (pow((i - halfN), 2) + pow((j - halfN), 2))) {
                    infSetTwo[i][j] = 1;
                }
            }
        }
    }

    private static void increaseInfSetFirstCondition(double rm) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (infSetOne[i][j] == 1) {
                    for (int i2 = 0; i2 < n; i2++) {
                        for (int j2 = 0; j2 < n; j2++) {
                            if (pow(rm / m, 2) >= (pow((i2 - i), 2) + pow((j2 - j), 2))) {
                                infSetIncOne[i2][j2] = 1;
                            }
                        }
                    }
                }
            }
        }
    }

    private static void increaseInfSetSecondCondition(double rm) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (infSetTwo[i][j] == 1) {
                    //  System.out.print("[ " + i + "]" + "[" + j + " ]" );
                    for (int i2 = 0; i2 < n; i2++) {
                        for (int j2 = 0; j2 < n; j2++) {
                            if (pow(rm / m, 2) >= (pow((i2 - i), 2) + pow((j2 - j), 2))) {
                                infSetIncTwo[i2][j2] = 1;
                            }
                        }
                    }
                }
            }
        }
    }

    private static void copyInfSetFC() {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (infSetIncOne[i][j] == 1) {
                    infSetOne[i][j] = 1;
                }
            }
        }
    }

    private static void copyInfSetSC() {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (infSetIncTwo[i][j] == 1) {
                    infSetTwo[i][j] = 1;
                }
            }
        }
    }

    private static void copyInfSetFCToMain() {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (infSetOne[i][j] == 1) {
                    infSet[i][j] = 1;
                }
            }
        }
    }

    private static void copyInfSetSCToMain() {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (infSetTwo[i][j] == 1) {
                    infSet[i][j] = 1;
                }
            }
        }
    }

    /*Draw the point on canvas*/
    public static void drawPoint(double x, double y) {
        gc.setFill(Color.ORANGE);
        double d = 2;
        gc.fillOval((Ox + (x - d)), (Oy - (y + d)), 2 * d, 2 * d);
    }

    public static void drawLine(double Ax, double Ay, double Bx, double By) {
        gc.setStroke(Color.RED);
        gc.strokeLine(Ox + Ax, Oy - Ay, Ox + Bx, Oy - By);
    }

    /*Method for draw infSet on canvas*/
    public static void drawInfSet() {
        for (int i = 0; i < infSet.length; i++) {
            for (int j = 0; j < infSet.length; j++) {
                if (infSet[i][j] > 0) {
                    if (infSet[i][j] == 1) gc.setFill(Color.NAVY);
                    if (infSet[i][j] == 3) gc.setFill(Color.RED);
                    if (infSet[i][j] == 5) gc.setFill(Color.GREEN);
                    if (infSet[i][j] == 7) gc.setFill(Color.BLUE);
                    if (infSet[i][j] == 8) gc.setFill(Color.YELLOW);
                    if (infSet[i][j] == 10) gc.setFill(Color.MAGENTA);
                    if (infSet[i][j] == 12) gc.setFill(Color.CYAN);
                    if (infSet[i][j] == 15) gc.setFill(Color.WHITE);
                    gc.fillOval(i * m - 0.5, j * m - 0.5, 1, 1);
                }
            }
        }
    }

    private void pursuerKnowSet(double xc, double yc) {
        for (int i = 0; i < infSet.length; i++) {
            for (int j = 0; j < infSet.length; j++) {
                if (pow(l / m + eps, 2) >= (pow((i - xc / m), 2) + pow((j - yc / m), 2))) {
                    infSet[i][j] = 0;
                    infSetOne[i][j] = 0;
                    infSetTwo[i][j] = 0;
                    infSetIncOne[i][j] = 0;
                    infSetIncTwo[i][j] = 0;
                }
            }
        }
    }

    private void pursuerKnowSetKCoverage(double xc, double yc, int id) {
        // xc, yc нужно поделить на m
        if (id == 1) mark = 3;
        if (id == 2) mark = 5;
        if (id == 3) mark = 7;
        for (int i = 0; i < infSet.length; i++) {
            for (int j = 0; j < infSet.length; j++) {
                if (pow(l + eps, 2) >= (pow(((i - 1) - xc), 2) + pow(((j + 1) - yc), 2))) {
                    if ((infSet[i][j] == 1) || (infSet[i][j] == mark)) {
                        infSet[i][j] = mark;
                    } else {
                        if (mark == 3) {
                            if ((infSet[i][j] == 5) || (infSet[i][j] == 7) || (infSet[i][j] == 12)) {
                                infSet[i][j] = infSet[i][j] + mark;
                            }
                        }
                        if (mark == 5) {
                            if ((infSet[i][j] == 3) || (infSet[i][j] == 7) || (infSet[i][j] == 10)) {
                                infSet[i][j] = infSet[i][j] + mark;
                            }
                        }
                        if (mark == 7) {
                            if ((infSet[i][j] == 3) || (infSet[i][j] == 5) || (infSet[i][j] == 8)) {
                                infSet[i][j] = infSet[i][j] + mark;
                            }
                        }
                    }
                }
            }
        }
    }

    /*Draw the circle with center(x,y)*/
    public static void drawCircle(double x, double y, double radius) {
        gc.setStroke(Color.BLUE);
        gc.strokeOval((Ox + x) - radius, (Oy - y) - radius, 2 * radius, 2 * radius);
    }

    /*Draw the circle with center(x,y)*/
    public static void drawCircleForPursuer(double x, double y, double radius, int mark) {
        switch (mark) {
            case 1: {
                gc.setStroke(Color.RED);
                break;
            }
            case 2: {
                gc.setStroke(Color.GREEN);
                break;
            }
            case 3: {
                gc.setStroke(Color.BLUE);
                break;
            }
        }

        gc.strokeOval((Ox + x) - radius, (Oy - y) - radius, 2 * radius, 2 * radius);
    }

    /* Method which finds coordinates of point on circle*/
    private static double[][] findPointsOnCircle(double Sx, double Sy, double radius, double kangle) {
        double Hx = Sx + sin(kangle) * radius;
        double Hy = Sy - cos(kangle) * radius;
        double Vx = Sx - cos(kangle) * radius;
        double Vy = Sy + sin(kangle) * radius;
        return new double[][]{{Hx, Hy}, {Vx, Vy}};
    }

    /* Method for finding coordinates of point on circle and line*/
    // point center of circle in this method (0,0)
    private static double[][] findPointsOnCircleAndLine(double k, double Ax, double Ay, double r) {
        double a = 1 + pow(k, 2);
        double b = 2 * (k * Ay - pow(k, 2) * Ax);
        /*   double c = pow((k * Ax - Ay), 2) - pow(r, 2);*/
        double c = pow(Ax * k, 2) - 2 * Ay * k * Ax + pow(Ay, 2) - pow(r, 2);
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

    /*Calculates distance between two points*/
    private static double distanceBetweenTwoPoints(double Ax, double Ay, double Bx, double By) {
        return (sqrt(pow(Bx - Ax, 2) + pow(By - Ay, 2)));
    }

    @FXML
    public void draw() {
        l = Double.parseDouble(valuel.getText());
        r = Double.parseDouble(valuer.getText());
        n = Integer.parseInt(valuen.getText());
        alpha = Double.parseDouble(valuealpha.getText());
        beta = Double.parseDouble(valuebeta.getText());
        dt = Double.parseDouble(valuedt.getText());

        /* Field of play w x h*/
        w = (int) canvas.getWidth();
        h = (int) canvas.getHeight();
        m = (double) w / n;

        Ox = (int) canvas.getWidth() / 2;
        Oy = (int) canvas.getHeight() / 2;

        infSetIncOne = new int[n][n];
        infSetIncTwo = new int[n][n];

        fillInfSet(r, n, m);
        drawInfSet();
        drawCircle(0, 0, r);


        for (int i = 0; i < pursuersList.size(); i++) {
            double x = pursuersList.get(i).getXc();
            double y = pursuersList.get(i).getYc();
            drawCircle(x, y, l);
        }
        /*
        tStep - фактическое количество шагов
        tStepOne - шаг по горизонтали и вертикали(приращение +1), обнуляется тогда когда tStepOne*beta*dt>m
        tStepTwo - шаг по диагонали(приращение +1), обнуляется тогда когда tStepTwo*beta*dt>sqrt(2)*m
         */
        tStep = 0;
        tStepOne = 0;
        tStepTwo = 0;
        startButton.setVisible(true);
        startButton.requestFocus();
    }

    @FXML
    void startGame() {
        l = Double.parseDouble(valuel.getText());
        r = Double.parseDouble(valuer.getText());
        n = Integer.parseInt(valuen.getText());
        alpha = Double.parseDouble(valuealpha.getText());
        beta = Double.parseDouble(valuebeta.getText());
        dt = Double.parseDouble(valuedt.getText());
        int w = (int) canvas.getWidth();
        int h = (int) canvas.getHeight();
        m = (double) w / n;

        Cx = -2 * sqrt(r * l);
        Cy = r - l;

        gc = canvas.getGraphicsContext2D();

        stepArr = new int[pursuersList.size()]; // шаг для каждого ищущего
        Arrays.fill(stepArr, 0);


        startButton.setDisable(false);
        pauseButton.setVisible(true);
        pauseButton.requestFocus();

        timerState = true;
        timer = new Timer();
        timerTask = new TimerTask() {
            @Override
            public void run() {
                gc.clearRect(0, 0, canvas.getWidth(), canvas.getHeight());
           /*     int Xamin = (int) (n / 2 - (r + tStep * dt * beta) / m);
                int Xamax = (int) (n / 2 + (r + tStep * dt * beta) / m);
                int Ybmin = (int) (n / 2 - (r + tStep * dt * beta) / m);
                int Ybmax = (int) (n / 2 + (r + tStep * dt * beta) / m);*/


                if (beta != 0) {
                    if (dt * tStepOne * beta >= m) {
                        //System.out.println("dt*tStepOne*beta =" + dt * tStepOne * beta);
                        increaseInfSetFirstCondition(dt * tStepOne * beta);
                        copyInfSetFC();
                        copyInfSetFCToMain();
                        tStepOne = 0;
                    }
                    if ((dt * tStepTwo * beta) >= (sqrt(2) * m)) {
                        //System.out.println("dt * tStepTwo * beta) =" + dt * tStepTwo * beta));
                        increaseInfSetSecondCondition(dt * tStepTwo * beta);
                        copyInfSetSC();
                        copyInfSetSCToMain();
                        tStepTwo = 0;
                    }
                }


                for (int i = 0; i < pursuersList.size(); i++) {
                    pursuer = pursuersList.get(i);
                    try {
                        if (pursuer.getIter() != -1) {
                            singleSegment = pursuer.getSegments().get(pursuer.getIter());
                            if (stepArr[i] * dt < singleSegment.getT()) {
                                pursuer.setXc(pursuer.getXc() + cos(toRadians(singleSegment.getAngle())) * alpha * dt);
                                pursuer.setYc(pursuer.getYc() + sin(toRadians(singleSegment.getAngle())) * alpha * dt);

//                                System.out.println("Фактические координаты");
//                                System.out.println("Xc " + pursuer.getXc() + " Yc " + pursuer.getYc());
                                stepArr[i]++;
                            } else {
                                pursuer.incIter();
                                stepArr[i] = 0;
                            }
                        }
                    } catch (IndexOutOfBoundsException ie) {
                        timerTask.cancel();
                        timer.cancel();
                        timer.purge();
                        System.out.println(ie.getMessage());
                    }
                    drawCircleForPursuer(pursuer.getXc(), pursuer.getYc(), l, pursuer.getId());
                    if (kCoverageFlag.isSelected()) {
                        pursuerKnowSetKCoverage(Ox + pursuer.getXc(), Oy - pursuer.getYc(), pursuer.getId());
                    } else {
                        // 0903
                        pursuerKnowSet(Ox + pursuer.getXc(), Oy - pursuer.getYc());

                       /* System.out.println("Координаты хранящиеся в инфо об объекте");
                        System.out.println("Xc " + pursuer.getXc() + " Yc " + pursuer.getYc());
                        System.out.println("Координаты передаваемые в метод pursuerKnowSet()");
                        // поделить значения на m
                        System.out.println("Ox " + Ox + " Oy " + Oy);
                        System.out.println("Xc" + (Ox + pursuer.getXc())/m + " Yc " + (Oy - pursuer.getYc())/m);*/
                    }
                }

                drawCircle(0, 0, (r + dt * beta * tStep));
                drawInfSet();
                //   System.out.println("tStep = " + tStep + " tStepOne = " + tStepOne + " tStepTwo " + tStepTwo);
                //  System.out.println("r + dt*beta*tStep =" + (r + dt * beta * tStep));
                tStep++;
                tStepOne++;
                tStepTwo++;
            }
        };
        timer.scheduleAtFixedRate(timerTask, 0, 500);
    }

    @FXML
    void openDialogPursuerProp(ActionEvent event) throws IOException {
        FXMLLoader fxmlLoader = new FXMLLoader(getClass().getResource("PursuerDialog.fxml"));
        Parent parent = fxmlLoader.load();
        r = Double.parseDouble(valuer.getText());
        l = Double.parseDouble(valuel.getText());
        alpha = Double.parseDouble(valuealpha.getText());

//        Cx = -2 * sqrt(r * l);
//        Cy = r - l;
        Cx = -(r + l);
        Cy = 0;
        int i = tvPursuers.size() + 1;

        PursuerDialogController pdc = fxmlLoader.<PursuerDialogController>getController();
        pdc.setObservableList(tvPursuers);
        pdc.setPursuerList(pursuersList);
        pdc.setParameters(i, Cx, Cy, l, alpha);

        Scene scene = new Scene(parent);
        Stage stage = new Stage();
        stage.setTitle("Add Pursuer");
        stage.initModality(Modality.APPLICATION_MODAL);
        stage.setScene(scene);
        stage.setResizable(false);
        stage.showAndWait();
    }

    @FXML
    void findKAndSet(ActionEvent event) throws IOException {
        r = Double.parseDouble(valuer.getText());
        l = Double.parseDouble(valuel.getText());
        alpha = Double.parseDouble(valuealpha.getText());
        beta = Double.parseDouble(valuebeta.getText());
//        if(beta==0){
//
//        }
//        else{
//
//        }
        Cx = -(r + l);
        Cy = 0;
        if (r <= l) {
            /* Add First Pursuer in ArrayList*/
            ArrayList<Segment> sForFirst = new ArrayList<>();
            sForFirst.add(new Segment(1, 1000, 0));
            Pursuer pursuer = new Pursuer(1, Cx, Cy, l, alpha, sForFirst);
            pursuersList.add(pursuer);
            tvPursuers.add(pursuer);

        } else {
            int k = 1;
            // Информация о первом игроке
            double phiOne = firstPursuer(r, l)[0];
            double Sx = findPointOnCircleLow(r, l, phiOne)[0];
            double Sy = findPointOnCircleLow(r, l, phiOne)[1];

            double Ux = findPointOnCircleUp(r, l, phiOne)[0];
            double Uy = findPointOnCircleUp(r, l, phiOne)[1];


            // Нахождение двух точке пересечения нижней прямой с окружностью
            double Ax = findPointsOnCircleAndLine(tan(phiOne), Sx, Sy, r)[1][0];
            double Ay = findPointsOnCircleAndLine(tan(phiOne), Sx, Sy, r)[1][1];

            double Bx = findPointsOnCircleAndLine(tan(phiOne), Sx, Sy, r)[0][0];
            double By = findPointsOnCircleAndLine(tan(phiOne), Sx, Sy, r)[0][1];
            ArrayList<Segment> sForFirst = new ArrayList<>();
            sForFirst.add(new Segment(1, 1000, toDegrees(phiOne)));
            Pursuer pursuer = new Pursuer(k, Cx, Cy, l, alpha, sForFirst);
            pursuersList.add(pursuer);
            tvPursuers.add(pursuer);
            k++;
            double phiNext = calcArc(Bx, By, Cx, Cy, l)[1];
            double SxNext = findPointOnCircleLow(r, l, phiNext)[0];
            double SyNext = findPointOnCircleLow(r, l, phiNext)[1];
            double AxNext = findPointsOnCircleAndLine(tan(phiNext), SxNext, SyNext, r)[1][0];
            double AyNext = findPointsOnCircleAndLine(tan(phiNext), SxNext, SyNext, r)[1][1];
            double BxNext = findPointsOnCircleAndLine(tan(phiNext), SxNext, SyNext, r)[0][0];
            double ByNext = findPointsOnCircleAndLine(tan(phiNext), SxNext, SyNext, r)[0][1];
            double UxNext = findPointOnCircleUp(r, l, phiNext)[0];
            double UyNext = findPointOnCircleUp(r, l, phiNext)[1];
            ArrayList<Segment> sPursuerNext = new ArrayList<>();
            sPursuerNext.add(new Segment(1, 1000, toDegrees(phiNext)));
            Pursuer pursuerNext = new Pursuer(k, Cx, Cy, l, alpha, sPursuerNext);
            pursuersList.add(pursuerNext);
            tvPursuers.add(pursuerNext);
            while(true){
                k++;
                phiNext = calcArc(BxNext, ByNext, Cx, Cy, l)[1];
                SxNext = findPointOnCircleLow(r, l, phiNext)[0];
                SyNext = findPointOnCircleLow(r, l, phiNext)[1];
                AxNext = findPointsOnCircleAndLine(tan(phiNext), SxNext, SyNext, r)[1][0];
                AyNext = findPointsOnCircleAndLine(tan(phiNext), SxNext, SyNext, r)[1][1];
                BxNext = findPointsOnCircleAndLine(tan(phiNext), SxNext, SyNext, r)[0][0];
                ByNext = findPointsOnCircleAndLine(tan(phiNext), SxNext, SyNext, r)[0][1];
                UxNext = findPointOnCircleUp(r, l, phiNext)[0];
                UyNext = findPointOnCircleUp(r, l, phiNext)[1];
                sPursuerNext = new ArrayList<>();
                sPursuerNext.add(new Segment(1, 1000, toDegrees(phiNext)));
                pursuerNext = new Pursuer(k, Cx, Cy, l, alpha, sPursuerNext);
                pursuersList.add(pursuerNext);
                tvPursuers.add(pursuerNext);

                if (Double.isNaN(AxNext) || Double.isNaN(AyNext)) {
                    break;
                }
            }

        }

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

    @FXML
    void openDialogStrategyProp(ActionEvent event) throws IOException {
        FXMLLoader fxmlLoader = new FXMLLoader(getClass().getResource("StrategyDialog.fxml"));
        Parent parent = fxmlLoader.load();
        r = Double.parseDouble(valuer.getText());
        l = Double.parseDouble(valuel.getText());
        alpha = Double.parseDouble(valuealpha.getText());

        int s = Integer.parseInt(String.valueOf(numberOfStrategy.getValue()));
        int k = Integer.parseInt(String.valueOf(numberOfPursuers.getValue()));

        StrategyDialogController sdc = fxmlLoader.<StrategyDialogController>getController();
        sdc.setObservableList(tvPursuers);
        sdc.setPursuerList(pursuersList);
        sdc.setparameters(s, k, r, l, alpha, beta);

        Scene scene = new Scene(parent);
        Stage stage = new Stage();
        stage.setTitle("Strategy Settings");
        stage.initModality(Modality.APPLICATION_MODAL);
        stage.setScene(scene);
        stage.setResizable(false);
        stage.showAndWait();
    }


    @Override
    public void initialize(URL location, ResourceBundle resources) {

        gc = canvas.getGraphicsContext2D();
        Ox = (int) canvas.getWidth() / 2;
        Oy = (int) canvas.getHeight() / 2;
        colId.setCellValueFactory(new PropertyValueFactory<>("id"));
        colX.setCellValueFactory(new PropertyValueFactory<>("xc"));
        colY.setCellValueFactory(new PropertyValueFactory<>("yc"));
        colSegments.setCellValueFactory(new PropertyValueFactory<>("segments"));
        tableP.setItems(tvPursuers);
        tableP.setPlaceholder(new Label("List of pursuers is empty"));
        startButton.setVisible(false);
    }

    @FXML
    public void stopTimer(ActionEvent event) {
        timerState = false;
        timerTask.cancel();
        timer.cancel();
    }
}
