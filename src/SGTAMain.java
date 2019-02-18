import java.util.Arrays;

public class SGTAMain {


    //for the static method that does everything, have it take the arguments of the extrema of the antibiotic concn
    //and work out alpha in that

    public static void main(String[] args){
     /*   int[] a = {1, 2, 3};
        int[] b = a;

        b[0] = 99;

        int j = 4;
        j+= -5;

        System.out.println(Arrays.toString(a));
        System.out.println(Arrays.toString(b));
        System.out.println(String.valueOf(j));*/

        double dt = 0.01;
        double dx = 0.1;
        double D = 0.0108;
        double c_max = 10.5;
        //System.out.println(c_max);

        BioSystem.exponentialGradient_spatialAndGRateDistributions(c_max, D, dt, dx);
    }


}




