import javax.tools.Tool;
import java.util.Arrays;
import java.util.Random;

public class BioSystem {

    private int L, K;
    private double s, s_max;
    private double c, alpha, timeElapsed, dt, dx;
    private double D; //diffusion constant

    private Microhabitat[] microhabitats;
    Random rand = new Random();
    private int initialPop = 100;

    public BioSystem(int L, double S, double D, double alpha, double dt, double dx){

        this.L = L;
        this.s = S;
        this.s_max = S;
        this.D = D;
        this.alpha = alpha;
        this.dt = dt;
        this.dx = dx;
        this.microhabitats = new Microhabitat[L];
        this.timeElapsed = 0.;

        for(int i = 0; i < L; i++){
            double c_i = Math.exp(alpha*(double)i) - 1.;
            microhabitats[i] = new Microhabitat(S, c_i);
        }

        for(int i = 0; i < 20; i++){
            microhabitats[i].fillWithWildType(20-i);
        }
    }

    public int getL(){
        return L;
    }
    public double getS(){return s;}
    public double getTimeElapsed(){
        return timeElapsed;
    }

    public void updateMicrohabitats(Microhabitat[] updatedSystem){
        //this updates the system with the new values for S and N
        this.microhabitats = updatedSystem;
    }

    public Microhabitat getMicrohabitats(int i){
        return microhabitats[i];
    }

    public double getTotalN(){
        double runningTotal = 0;
        for(Microhabitat m : microhabitats) {
            runningTotal += m.getN_alive()+m.getN_dead();
        }
        return runningTotal;
    }

    public double getCurrentLivePopulation(){
        double runningTotal = 0;
        for(Microhabitat m : microhabitats) {
            runningTotal += m.getN_alive();
        }
        return runningTotal;
    }

    public double[] getLiveSpatialDistributionArray(){
        double[] mh_pops = new double[L];
        for(int i = 0; i < L; i++){
            mh_pops[i] = microhabitats[i].getN_alive();
        }
        return mh_pops;
    }

    public double[] getDeadSpatialDistributionArray(){
        double[] mh_pops = new double[L];
        for(int i = 0; i < L; i++){
            mh_pops[i] = microhabitats[i].getN_dead();
        }
        return mh_pops;
    }


    public double[] getGrowthRatesArray(){
        double[] mh_gRates = new double[L];
        for(int i = 0; i < L; i++){
            mh_gRates[i] = microhabitats[i].replication_or_death_rate();
        }
        return mh_gRates;
    }

    public double[] getNutrientsArray(){
        double[] mh_S = new double[L];
        for(int i = 0; i < L; i++){
            mh_S[i] = microhabitats[i].getS();
        }
        return mh_S;
    }

    public void performAction(){

        Microhabitat[] newSystem = microhabitats.clone();

        for(int i = 0; i < L; i++){

            double s_i = microhabitats[i].getS();
            double nalive_i = microhabitats[i].getN_alive();
            double ndead_i = microhabitats[i].getN_dead();

            if(microhabitats[i].itIsGrowing()){

                double growthRate_i = microhabitats[i].replication_rate_noisy();
                s_i += -dt*nalive_i*growthRate_i;


                if(i==0){
                    nalive_i += dt*(D*(microhabitats[i+1].getN_alive() + nalive_i - 2*nalive_i)/(dx*dx) + nalive_i*growthRate_i);
                }else if(i==L-1){
                    nalive_i += dt*(D*(nalive_i + microhabitats[i-1].getN_alive() - 2*nalive_i)/(dx*dx)+ nalive_i*growthRate_i);
                }else{
                    nalive_i += dt*(D*(microhabitats[i+1].getN_alive() + microhabitats[i-1].getN_alive() - 2*nalive_i)/(dx*dx)
                            + nalive_i*growthRate_i);
                }

            }else{
                double deathRate_i = microhabitats[i].death_rate(); //this is a -ve value
                //System.out.println("drate: "+deathRate_i);
                //growth rate is -ve so no change in nutrients
                if(i==0){
                    double deltaN = dt*(D*(microhabitats[i+1].getN_alive() + nalive_i - 2*nalive_i)/(dx*dx) + nalive_i*deathRate_i);
                    //if(deltaN < 0.) System.out.println("DEATH: "+deltaN);
                    nalive_i += deltaN;
                    ndead_i += Math.abs(deltaN);
                }else if(i==L-1){
                    double deltaN = dt*(D*(nalive_i + microhabitats[i-1].getN_alive() - 2*nalive_i)/(dx*dx) + nalive_i*deathRate_i);
                    //if(deltaN < 0.) System.out.println("DEATH: "+deltaN);
                    nalive_i += deltaN;
                    ndead_i += Math.abs(deltaN);
                }else{
                    double deltaN = dt*(D*(microhabitats[i+1].getN_alive() + microhabitats[i-1].getN_alive() - 2*nalive_i)/(dx*dx)
                            + nalive_i*deathRate_i);
                    //if(deltaN < 0.) System.out.println("DEATH: "+deltaN);
                    nalive_i += deltaN;
                    ndead_i += Math.abs(deltaN);
                }
            }
            newSystem[i].setS(s_i);
            newSystem[i].setN_alive(nalive_i);
            newSystem[i].setN_dead(ndead_i);
        }


        updateMicrohabitats(newSystem);

        timeElapsed += dt;
    }



    public static void exponentialGradient_spatialAndGRateDistributions(double biggestC, double D, double dt, double dx){

        int L = 5000, nReps = 1;
        int nTimeMeasurements = 20;

        double duration = 2000.;
        double interval = duration/(double)nTimeMeasurements;

        double alpha = BioSystem.calculateAlpha(L, biggestC);
        int S = 500;

        String filename_alive = "SGTA_death-alpha="+String.valueOf(alpha)+"-aliveDistribution-continuum-NOISY";
        String filename_dead = "SGTA_death-alpha="+String.valueOf(alpha)+"-deadDistribution-continuum-NOISY";
        String filename_gRate = "SGTA_death-alpha="+String.valueOf(alpha)+"-gRateDistribution-continuum-NOISY";
        String filename_nutrients = "SGTA_death-alpha="+String.valueOf(alpha)+"-nutrientDistribution-continuum-NOISY";


        double[][][] allN_alive = new double[nReps][][];
        double[][][] allN_dead = new double[nReps][][];
        double[][][] allGRates = new double[nReps][][];
        double[][][] allNutrients = new double[nReps][][];


        for(int r = 0; r < nReps; r++){

            boolean alreadyRecorded = false;

            double[][] alivePopsOverTime = new double[nTimeMeasurements+1][];
            double[][] deadPopsOverTime = new double[nTimeMeasurements+1][];
            double[][] gRatesOverTime = new double[nTimeMeasurements+1][];
            double[][] nutrientsOverTime = new double[nTimeMeasurements+1][];
            int timerCounter = 0;


            BioSystem bs = new BioSystem(L, S, D, alpha, dt, dx);

            while(bs.timeElapsed <= duration+0.2*interval){

                if((bs.getTimeElapsed()%interval >= 0. && bs.getTimeElapsed()%interval <= 0.1*interval) && !alreadyRecorded){

                    System.out.println("rep: "+r+"\ttime elapsed: "+String.valueOf(bs.getTimeElapsed()) );
                            //"\ttotal N: "+String.valueOf(bs.getTotalN()));

                    //System.out.println(Arrays.toString(bs.getNutrientsArray()));

                    alivePopsOverTime[timerCounter] = bs.getLiveSpatialDistributionArray();
                    deadPopsOverTime[timerCounter] = bs.getDeadSpatialDistributionArray();
                    gRatesOverTime[timerCounter] = bs.getGrowthRatesArray();
                    nutrientsOverTime[timerCounter] = bs.getNutrientsArray();

                    alreadyRecorded = true;
                    timerCounter++;
                    System.out.println(timerCounter);
                }
                if(bs.getTimeElapsed()%interval >= 0.1*interval) alreadyRecorded = false;


                bs.performAction();
            }

            //System.out.println(Arrays.deepToString(alivePopsOverTime));
            allN_alive[r] = alivePopsOverTime;
            allN_dead[r] = deadPopsOverTime;
            allGRates[r] = gRatesOverTime;
            allNutrients[r] = nutrientsOverTime;



        }
        //System.out.println(Arrays.toString(allN_alive));
        double[][] averagedAlivePopDistributions = Toolbox.averagedResults(allN_alive);
        double[][] averagedDeadPopDistributions = Toolbox.averagedResults(allN_dead);
        double[][] averagedGRateDistributions = Toolbox.averagedResults(allGRates);
        double[][] averagedNutrientDistributions = Toolbox.averagedResults(allNutrients);

        //print the live and dead results to two separate files, then just join them together in
        //gnuplot or something
        Toolbox.printAveragedResultsToFile(filename_alive, averagedAlivePopDistributions);
        Toolbox.printAveragedResultsToFile(filename_dead, averagedDeadPopDistributions);
        Toolbox.printAveragedResultsToFile(filename_gRate, averagedGRateDistributions);
        Toolbox.printAveragedResultsToFile(filename_nutrients, averagedNutrientDistributions);

    }


    public static double calculateAlpha(int L, double c_max){

        return Math.log(c_max+1)/L;
    }

}
