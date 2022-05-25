
public class TestClass {
    int a;
    static int b;

    void f1() {
        ++a;
    }

    static void f2() {
        ++b;
    }

    public static double []compareTwoNumbers(double a, double b){
        if(a<b){
            return new double[]{Double.NaN, Double.NaN};
        }
        else{
            return new double[]{a,b};
        }
    }
    public static void main(String[] args) {
        System.out.println(compareTwoNumbers(5,12)[0] + " " + compareTwoNumbers(5,12)[1]);
    }
}
