#include <cstdlib>
#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>
#include "MersenneTwister.h"
#include <stdio.h>

using namespace std;
using std::atoi;

// Random number generator:
MTRand rnd;


int main()	{


    double b, d, c, s, h, alpha, X, N, X0, N0, g, dt, dBt1, dBt2, Nt, Nhar, Nmin, typesim;
    int id, n, rep;
    int T, Tfix;
    bool step, mode;


    int no = 1;
    int end  = 0;


    FILE * fp;
    fp = fopen("param.txt","r");
    if (!fp){
        cout << "Parameter file doesn't exist!\n\n";
    }
    int x;
    do{
        do {x = fgetc(fp);} while (!((x == '*') || (x == EOF)));
        // Lines with parameter sets must begin with *
        if (x == EOF)
        {
            cout << "\nEnd of input file\n";
            end = 1;
        }
        else
        {

            fscanf(fp,"%d ",&typesim);// 0 Our model with genetic feedback, 1 Simulations without feed-back (population size remains stochastic however), and 2 Simulations with constant population size
            fscanf(fp,"%d ",&mode);// 0 writes only means for each of the n repetitions at the end of simulation, the last line representing the summary statistics over all n repetitions, 1 writes entire trajectory for 1 simulation only (n = 1)
            fscanf(fp,"%lf ",&N0);// initial population size
            fscanf(fp,"%lf ",&X0);// frequency of allele a
            fscanf(fp,"%lf ",&b);//birth rate
            fscanf(fp,"%lf ",&d);//death rate
            fscanf(fp,"%lf ",&c);//competition
            fscanf(fp,"%lf ",&alpha);//self-fertilisation rate
            fscanf(fp,"%lf ",&g);//parameter gamma for scaling the speed of birth and death rates. In order to obtain a model as close as possible to the Wright-Fisher diffusion, this parameter is set to 0.5
            fscanf(fp,"%lf ",&s);//selection coefficient
            fscanf(fp,"%lf ",&h);//dominance
            fscanf(fp,"%lf ",&dt);//Time step, set to 10e-04 for large populaiton and growth rates, and 10e-05 otherwise
            fscanf(fp,"%d ",&n);// number of repeats for each parameter set
            fscanf(fp,"%d ",&id);// identity number defined by the used and added at the end of the output file so as to avoid over-writing simulation outputs with the same parameter values
            
           
            
            // creating name of output file (which indicates the parameter values):
            
            char nomFichier[256], nomF2[256], nomF3[256];
            stringstream nomF;
            
            nomF <<"Type"<<typesim<<"_Mode" << mode<<"_N" << N0 << "_X" << X0 <<"_b" << b<<"_d" << d << "_c"<< c << "_alpha"<< alpha << "_g" << g << "_s" << s << "_h"<< h <<"_dt"<< dt<<"_"<<id<<".txt";
            nomF >> nomFichier;
            ofstream fout(nomFichier);


            // Simulation:

            double Xfmean = 0; // probaility of fixation over n repeats
            double Xmean = 0; // mean frequency of X over n repeats
            double Tmean = 0; // Mean number of time steps
           double Tfmean = 0; // Mean time to fixation
           double Textmean = 0; // Mean time to extinction (conditional to extinction)
           double Tlmean = 0; // Mean time to loss
            double Tallmean = 0; // Mean time to homogenisation of pop (fixation or loss)
            double Nharmean =0; // mean harminic mean of all repeats
            double Nmean =0; // mean population size of all repeats
            double meanNhar =0; // mean 1 / population size of all repeats
            

            if (mode == 1) n = 1;
            int nf =0;
            int nl =0;
            int next = 0;
            int nall = 0;

            for(rep = 0; rep < n; rep++)
            {
                N = N0;
                Nmin = N0;
                X = X0;
                Nhar = 0;
                int T = 0; //Initialising time counter
                Tfix = 0;
                step =0;
                double Nm =0;
                do {
                    T++;

                    // Population size follows a Brownian motion. At step dt+1 N depends on population size at step dt and a random variable sampled from a Normal distribution centered on 0 and with variance dt. It is not defined for Type 2 because population size remains constant
                    if (typesim==0){
                    double dBt1= rnd.randNorm(0,sqrt(dt));
                    N+=sqrt(2*g*N)*dBt1+(b-d-c*N+(s/(2-alpha))*(h*(1-alpha)*4*X*(1-X)+alpha*X+(1-alpha)*2*X*X))*N*dt;
                    }
                    else if (typesim==1){
                        double dBt1= rnd.randNorm(0,sqrt(dt));
                        N+=sqrt(2*g*N)*dBt1+(b-d-c*N)*N*dt;
                    }

                    if(N <= 0){
                        N=0;
                        Textmean += T;
                        next ++;

                    }

                    if (N !=0){

                    Nhar += 1/N;
                    Nm += N;

                    // frequency X of allele A at step dt+1. Like for poopulation size, X follows a Brownian motion.
                    double dBt2 = rnd.randNorm(0,sqrt(dt));
                    X+=sqrt(2*g*X*(1-X)/((2-alpha)*N))*dBt2+(s/(2-alpha))*(1-X)*X*(2*h*(1-alpha)*(1-2*X)+alpha+2*(1-alpha)*X)*dt;
                    }
                    
                    if(X >= 0.9999) X=1;
                    if(X <= 0.0001) X=0;

                    if ((step==0)&&((X==0)||(X==1)))
                    {
                        step=1;
                        Tfix = T;

                    }

                    if (mode == 1)
                    fout<<X <<" "<< N << " "<< (T/Nhar)/((4*g)/(2-alpha))/*Harmonic effective population size*/<<endl;

                    if (N<Nmin)
                    {Nmin = N;}

                } while ((N>0)&&(step==0));


               
                if (mode == 0)
                fout<<X <<" "<< Nm/T <<" "<< (T/Nhar)/((4*g)/(2-alpha))//Harmonic effective population size
                    << " "<< dt*T <<" "<< dt*Tfix <<" "<<N<<" "<<Nmin<<endl;

                Xmean +=X;

                if (X==1){
                    Tfmean += Tfix;
                    nf++;
                    Tallmean += Tfix;
                    nall++;
                    Xfmean +=X;
                }
                if (X==0){
                    Tlmean += Tfix;
                    nl++;
                    Tallmean += Tfix;
                    nall++;
                    Xfmean +=X;
                }

                Tmean +=T;
                Nharmean += T/Nhar;
                Nmean += Nm/T;
                meanNhar += Nhar;


            }
            
            
            if (mode == 0)
            fout<<Xfmean/nall <<" "<< Nmean/n<<" "<< Nharmean/n //Mean demographic harmonic population size over all simulations
            <<" "<< (Nharmean/n)/((4*g)/(2-alpha))//Mean effective harmonic population size over all simulations
             <<" "<<
            
            (Tmean)/(meanNhar*((4*g)/(2-alpha)))  // Our proposed effective population size
                << " "<< dt*(Tfmean/nf)<<" "<< dt*(Tlmean/nl)<< " "<< dt*(Tallmean/nall)<<" "<< dt*(Textmean/next) <<endl;

        }

    }while(end!=1);

    // Closes files:
    fclose(fp);
}




