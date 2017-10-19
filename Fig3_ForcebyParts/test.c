/*THE PURPOSE OF THIS PROGRAM IS TO STUDY STRESS RELAXATION OVER A GENERALIZED VISCOELASTIC SURFACE USING PRONY SERIES*/

/*HEADERS*/
#include<stdio.h>
#include<math.h>
/*HEADERS*/

/* CONSTANTS */
#define 	pi	3.1415926535
#define     n  12               //Number of Maxwell arms in the generalized model
#define     Nprint  50000      //Number of data points printed
#define     NN  25      //number of fundamental periods used to calculate phase and amplitude
/* CONSTANTS */

/*DECLARATION OF VARIABLES*/
//variables used to write txt outputfiles 
FILE		*Out1, *Out2, *Out3;						// Pointer for output files
FILE        *Inp;                          //Pointer for input file
char        input_file[100];
char		output_file_1[100], output_file_2[100], output_file_3[100];	
char        output_compu[100];						// Output file root name
char        output_analytical[100];                 // Output file root name

double      a, H;   //inter atomic distance and Hamaker constant

double		eta[n];  //array storing the values of viscosity of the dash-pots of each of the Maxwell arms
double      Xc[n], Xc_rate[n];
double      G_xb_xc[n], G_xc[n];
double 		Gg, Ge, G_loss, G_storage, G_absolute;    //Glassy modulus, equilibrium, loss, storage, and absolute moduli, respectively
double      alfa, R, nu;                  //cell constant, radius of tip, Poisson's ratio of sample
double      Z, A, omega, f1, period1, dt, simultime;  
double      sum_transient_E, E_analytical, E_storage, E_loss, E_equilibrium, E_transient, E_analytical_vdw;
double      sum_G_xc, sum_G_xb_xc;
double      sum_Gstorage;  //variable used to calculate storage modulus
double      z_base_counter, z_base, z_base_max;
double      time, TipPos, TipV, Fts, Xb, NP1, Energy;
double      t_in, t_fin, t_fin_2;
double      TipPrescribed, CosRef;
double      sum_transient, F_analytical, F_Ge, F_DMA, F_transient;
double      MaxTipPos, MinTipPos;   //variables to calculate maximum and minimum tip positions, respectively
double		PeakForce;  

//Print Counters
long int	print_counter1;		// Variable that controls the printing step of SumA and SumB
double      startprintF, stopprintF, startCalcAmp, startvirial;            // Variables that control the print length of the ForceSLS file
int         i, j,k, print_counter;

//Verlet variables
double		v1, z1, z1_old, z1_new, v2, z2, z2_old, z2_new, v3, z3, z3_old, z3_new;     	// Initial velocity and velocity at time t 
double		a1, a2, a3;						            // acceleration at time t for the various eigenmodes
double		TipPos, TipPos_old;							        	// Instantaneous tip position
double		TipV;								        // Instantaneous tip velocity

//Variables related to dynamics of the cantilever
double		A1, A2, A3, A2_free;						                            // Target free amplitude of the first eigenmode
double      F_drive1, F_drive2, F_drive3;                               // drive amplitude of the force applied to the tip to reach the target amplitude 
double      Asp;                                                        // setpoint amplitude  
double		f1, f1_for, w1, w1_free, f2, w2, f3, w3;				            // Eigenfrequency of the first eigenmode and excitation frequency
double		k_m1, Q1, mass, period1, k_m2, Q2, period2, k_m3, Q3, period3, period_timestep;		        // Stiffness of the first mode of the cantilever

//variables to calculate phase, amplitudes, virials an power dissipated

double		SumA, SumB, SumC, SumD;				// Variables used for calculating the phase and amplitude
double		Output_Phase1, Output_Phase2;			// Phases of individual eigenmodes
double		Output_Amp1, Output_Amp2;				// Amplitudes of individual eigenmodes 
double      Ets1, Vts1, Ets2, Vts2;                 // Power dissipated and virial for the first two modes
double      Virial_compu, Av_Force_an, Av_Force_compu;
double      Zs;
long int	counterPh1, counterPh2;			      	// Counter variables used in the calculation of the phase
double      V_storage, V_loss, V_transient, V_eq, V_analytical, sum_transient_V, V_analytical_vdw;



int         NOL;   //number of lines of input file with operator coefficients
double     ch;  //variable storing current number that is being read in the input file
double     d;  //temporal variable containing the coefficient
int        x, y;
double mintau;
/*DECLARATION OF VARIABLES*/


/*BEGINING OF MAIN PROGRAM*/
main()
{
	mintau = 1.0e-9;
		//Reading the Prony coefficients from input files// Opening input file
    sprintf(input_file,"PIB.txt");
    Inp = fopen(input_file,"r");
	
	NOL = 0;
	do 
	{
		ch = fgetc(Inp);
		if(ch == '\n')
    	NOL++;
	} while (ch != EOF);
	
    printf("%d\n", NOL);
	fclose(Inp);
	
	double  Q[NOL][2];    
			
	sprintf(input_file,"PIB.txt");
    Inp = fopen(input_file,"r");
	
	/*STORING MATRIX CONTAINING MODULI AND RELAXATION TIMES*/
	for(x= 0; x < NOL; x++ ) 
	{
        for(y = 0; y < 2; y++ ) 
		
		{
			fscanf(Inp,"%lf",&d);
			Q[x][y]= d;
		}
		
    }	
	/*STORING MATRIX CONTAINING MODULI AND RELAXATION TIMES*/
	fclose(Inp);
	
	j = 0;
	k = 0;
	do 
	{
		if (Q[j][0] > mintau)
		{
			printf("%d\n", j);
			printf("%20.20f\n", Q[j][0]);
			k++;
		}
		j++;
	} while (j<NOL);
	
	double  tau[k];
	double  G[k];
	
	for (i=0; i< k; i++)
	{
		tau[i] = Q[i+(j-k)][0];
		G[i] = Q[i+(j-k)][1];
	}
	//End of reading coefficients from input file
		
	
	/*MODEL PARAMETERS*/  //Taken from Tobolski data for polyisobutylene, see Brinson, Polymer Engineering Science and Viscoelasticity p249.
//opening output file 3
sprintf(output_file_3,"Model_params.txt");
Out3 = fopen(output_file_3,"w");
  
// Printing the output file header of output file 3
fprintf(Out3,"Relaxation_time(s)\tG_modulus(Pa)\tGe(Pa)\n");  

Ge = 0.0;  //Equilibrium Modulus in Pascals
Gg = Ge; //initializing glassy modulus to value of equilibrium modulus


for (i=0; i < (k) ; i++)
{
	fprintf(Out3, "%20.20f\t%20.20f\t%5.5f\n", tau[i], G[i], Ge);	
}
fclose(Out3);
/*MODEL PARAMETERS*/	
	
printf("%d\n", j);
printf("%d\n", k);
	getch();
}//closing main program
            