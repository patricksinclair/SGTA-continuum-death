

public class Microhabitat {

    //no. of bacteria in the microhabitat
    private double N_alive, N_dead;
    private double c;
    private double S, S_max;

    //death and migration rates
    private double d = 0., b = 0.1;
    private double K_prime = 33.;

    public Microhabitat(double S, double c){
        this.S = S;
        this.S_max = S;
        this.c = c;
        this.N_alive = 0;
        this.N_dead = 0;
    }



    public double getS(){return S;}
    public double getC(){return c;}
    public double getN_alive(){return N_alive;}
    public double getN_dead(){return N_dead;}
    public double getN_tot(){return N_alive+N_dead;}
    public double getD(){return d;}
    public double getB(){return b;}

    //these voids are used to update the system via the diff eqns
    public void setS(double S){
        this.S = S;
    }
    public void setN_alive(double N_alive){
        this.N_alive = N_alive;
    }
    public void setN_dead(double N_dead){
        this.N_dead = N_dead;
    }


    public double beta(){

        double mu = S/(K_prime+S);
        double mu_max = S_max/(K_prime+S_max);
        return 1. + 9.*mu/mu_max;
    }

    public double phi_c(){
        return 1. - (6.*(c/beta())*(c/beta()))/(5. + (c/beta())*(c/beta()));
    }


    public boolean itIsGrowing(){
        //returns true if phi(c) > 0, i.e. growth rate, false if death
        return phi_c() > 0;
    }

    public double replication_or_death_rate(){

        double phi_c = phi_c();
        return (phi_c >= 0) ? phi_c*(S/(K_prime + S)) : phi_c;
    }

    public double replication_rate(){
        return phi_c()*(S/(K_prime + S));
    }

    public double replication_rate_noisy(){
        double noise = (Math.random()*(0.1 - (-0.1))) + (-0.1); //random noise between -0.1 and 0.1
        return phi_c()*(S/(K_prime + S))*(1. + noise);
    }

    public double death_rate(){
        //use this if phi(c) is negative
        //it will therefore return a negative value
        return phi_c();
    }


    public void fillWithWildType(int initalPop){
        N_alive+=initalPop;
    }

    public void consumeANutrient(){
        if(S > 0) S--;
    }

    public void addABacterium(){N_alive++;}
    public void removeABacterium(){N_alive--;}

    public void replicateABacterium(){N_alive++;}

    public void killABacterium(){
        if(N_alive > 0){
            N_alive--;
            N_dead++;
        }
    }






}


