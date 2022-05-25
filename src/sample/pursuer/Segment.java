package sample.pursuer;

public class Segment {
    private int id;
    private double t;
    private double angle;
    public Segment(int id, double t, double angle) {
        this.id = id;
        this.t = t;
        this.angle = angle;
    }

    public int getId() {
        return id;
    }

    public double getT() {
        return t;
    }

    public double getAngle() {
        return angle;
    }

    @Override
    public String toString() {
        return "Segment{" +
                "id=" + id +
                ", t=" + t +
                ", angle=" + angle +
                '}';
    }

}
