package sample.pursuer;

import java.util.ArrayList;
import java.util.Iterator;

public class Pursuer {
    private int id;
    private double xc;
    private double yc;
    private double l;
    private double alpha;
    private ArrayList<Segment> segments;
    private int iter;

    public Pursuer(int id, double xc, double yc, double l, double alpha, ArrayList<Segment> segment) {
        this.id = id;
        this.xc = xc;
        this.yc = yc;
        this.l = l;
        this.alpha = alpha;
        this.segments = segment;
        this.iter = 0;
    }

    public int getId() { return id; }

    public double getXc() { return xc; }

    public double getYc() { return yc; }


    public double getAlpha() {
        return alpha;
    }

    public double getL() {
        return l;
    }

    public void setXc(double xc) {
        this.xc = xc;
    }

    public void setYc(double yc) {
        this.yc = yc;
    }

    public int getIter() {
        return iter;
    }

    public void incIter() {
        if(this.iter < this.segments.size()-1){
             this.iter++;
        }
        else{
            this.iter = -1;
        }
    }

    public ArrayList<Segment> getSegments() {
        return segments;
    }

    @Override
    public String toString() {
        return "Pursuer [" +
                "id = " + id +
                ", xc = " + xc +
                ", yc = " + yc +
                ", l = " + l +
                ", alpha = " + alpha +
                ", segments = " + segments +
                ']';
    }
}
