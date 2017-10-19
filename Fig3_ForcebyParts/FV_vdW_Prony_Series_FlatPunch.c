/*THE PURPOSE OF THIS PROGRAM IS TO STUDY STRESS RELAXATION OVER A GENERALIZED VISCOELASTIC SURFACE USING PRONY SERIES*/

/*HEADERS*/
#include<stdio.h>
#include<math.h>
/*HEADERS*/

/* CONSTANTS */
#define 	pi	3.1415926535
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
int         i, j,n, print_counter;

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
	n = 0;
	do 
	{
		if (Q[j][0] > mintau)
		{
			printf("%d\n", j);
			printf("%20.20f\n", Q[j][0]);
			n++;
		}
		j++;
	} while (j<NOL);
	
	double  tau[n];
	double  G[n];
	double		eta[n];  //array storing the values of viscosity of the dash-pots of each of the Maxwell arms
	double      Xc[n], Xc_rate[n];
	double      G_xb_xc[n], G_xc[n];
	
	for (i=0; i< n; i++)
	{
		tau[i] = Q[i+(NOL-n)][0];
		G[i] = Q[i+(NOL-n)][1];
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


for (i=0; i < (n) ; i++)
{
	fprintf(Out3, "%20.20f\t%20.20f\t%5.5f\n", tau[i], G[i], Ge);	
}
fclose(Out3);
/*MODEL PARAMETERS*/
	
	
	

/*DEFINING VALUES OF VISCOSITY FOR EACH DASHPOT, CALCULATING Gg*/
//Initializing values
sum_Gstorage = 0.0;
G_loss = 0.0;

for (i=0; i < (n) ; i = i+1)
{
	eta[i] = tau[i] * G[i];
	Gg = G[i]+Gg;
}
/*DEFINING VALUES OF VISCOSITY FOR EACH DASHPOT, CALCULATING Gg*/
	
	
	
	
	
	
	
	
	
/*INSERTING VALUE OF VARIABLES RELATED TO DYNAMICS OF CANTILEVER*/
k_m1 = 10.0;   
Q1 = 150;
Q2 = 450;
Q3 = 750;

A1 = 100;                    //target free amplitude of the first eigenmode
A1 = A1*1e-9;	     	    // Converting to meters
A2 = 0;
A2 = A2*1e-9;
A3 = 0;                     //target free amplitude of the third eigenmode
A3 = A3*1e-9;               // Converting to meters 
   
z_base_max = A1*1e9;        //cantilever base maximum position in nanometers
/*INSERTING VALUE OF VARIABLES RELATED TO DYNAMICS OF CANTILEVER*/

/*PARAMETERS CALCULATED FOR THE SECOND AND THIRD MODES FROM THE DATA INSERTED BY USER*/
//calculation of frequency of the three eigenmodes

f1 = 100.0e3;  //excitation frequency
omega = 2*pi*f1;   //first mode frequency in rad/s
f2 = 6.27*f1; // excitation frequency of the second eigenmode 
f3 = 17.6*f1; // excitation frequency of the second eigenmode

period1 = 1.0/f1;  
period2 = period1/6.27;
period3 = period1/17.6;

// Calculating the equivalent mass from the resonance frequency and the harmonic spring constants of second and third eigenmodes
mass = k_m1 / pow((f1*2*pi), 2);
k_m2 = k_m1*(f2/f1)*(f2/f1);
k_m3 = k_m1*(f3/f1)*(f3/f1);

//Calculation of the drive amplitude of the three Eigenmodes to achieve the desired free amplitude (When driven at Resonance!)
F_drive1 = A1*k_m1/Q1;
F_drive2 = A2*k_m2/Q2;
F_drive3 = A3*k_m3/Q3;

//Start of printing
startCalcAmp = period1*10.0*Q1;
startprintF = period1*10.0*Q1 + NN*period1;
simultime = startprintF + NN*period1;

//Defining timestep
dt = period1/2.0e4;
/*CALCULATION FOR THE SECOND AND THIRD MODES ENDS HERE*/  

/*PARAMETERS RELATED TO AFM DYNAMICS when given prescribed trajectory*/
NP1 = ((simultime-startprintF)/dt)*(1.0/Nprint);  //variable used for counting the print steps
/*PARAMETERS RELATED TO AFM DYNAMICS when given prescribed trajectory*/

/*CELL CONSTANT (APARATUS) PARAMETERS*/
a = 0.2e-9;  //Minimum intermolecular distance
H = 2.0e-19;  //Hamaker constant of sample
R = 10.0e-9; //cylindrical flat punch end diameter
nu = 0.5;  //material's Poisson ratio
alfa = 4.0*R/(1.0-nu);  //Cell constant converting stress/strain to force/displacement. This is derived form application of correspondence principle to elastic solution obtained by Sneddon 1965 (Int Journal of Engng Sci. vol3, pp.47-57), see eq 6.1
/*CELL CONSTANT (APARATUS) PARAMETERS*/






//opening output file 2
sprintf(output_file_2,"Summary.txt");
Out2 = fopen(output_file_2,"w");
  
// Printing the output file header of output file 2
fprintf(Out2,"Amp1_free(nm)\tk1(N/m)\tQ1\tRadiusPunch(nm)\tPoisonRation\tHammaker1e19\tZeq(nm)\tZr(nm)\tAmplitude1(nm)\tPhase1(deg)\tFreq1(kHz)\tEdiss_compu(aJ)\tVts_compu(aJ)\tAv_Fts_compu(nN)\tEdiss_Tamayo(aJ)\tVts_Lozano(aJ)\tAv_Fts_Lozano(aJ)\n");  

/*START OF FOR LOOP WHERE THERE IS GOING TO BE A CHANGE OF THE AMPLITUDE OF THE SECOND EIGENMODE*/ 
//for (z_base_counter = (z_base_max-10.0); z_base_counter > 0.1; z_base_counter = z_base_counter - 10.0)
for (z_base_counter = (z_base_max-60.0); z_base_counter > -160.0; z_base_counter = z_base_counter - 1000.0)
{
	z_base = z_base_counter*1e-9;
	
	/*WRITING NAME, OPENING, AND WRITING HEADER OF OUTPUT FILE WITH NUMERICAL RESULTS*/
	// Changing the name of the output file1
    sprintf(output_file_1,"%scompu_%.2f.txt",output_compu,z_base*1e9);
  
    // Opening the output file1
    Out1 = fopen(output_file_1,"w"); 

    // Printing the output file header of output file 1
    fprintf(Out1,"Time(us)\tTipPos(nm)\tForce(nN)\tSamplePos(nm)\tTipVel(nm/s)\n");
    /*END OF WRITING NAME, OPENING, AND WRITING HEADER OF OUTPUT FILE WITH NUMERICAL RESULTS*/
	
					
    /*INITIALIZING VERLET VARIABLES*/
    a1 = 0.0;
    a2 = 0.0;
    a3 = 0.0;
  
    v1 = 0.0;
    v2 = 0.0;
    v3 = 0.0;
  
    z1 = 0.0;
    z2 = 0.0;
    z3 = 0.0;
  
    z1_old = 0.0;
    z2_old = 0.0;
    z3_old = 0.0;
  
    z1_new = 0.0;
    z2_new = 0.0;
    z3_new = 0.0;
	TipPos = 0.0;
    TipPos_old = 0.0;
	TipV = 0.0;
	/*INITIALIZING VERLET VARIABLES*/
	
	/*INITIALIZING VARIABLES RELATED TO VIRIAL, DISS POWER AND CALCULATION OF AMPLITUDES, PHASE*/
	//Initialization of virial and dissipated power variables
    Ets1 = 0.0;
    Vts1 = 0.0;
    Ets2 = 0.0;
    Vts2 = 0.0;
   	
    //Initalization of phase and amplitude calculation variables
    counterPh1 = 0;
    counterPh2 = 0;
    print_counter1 = 1;
    SumA = 0.0;
    SumB = 0.0;
    SumC = 0.0;
    SumD = 0.0;
    Output_Amp1 = 0.0;
    Output_Amp2 = 0.0;
    Output_Phase1 = 0.0;
    Output_Phase2 = 0.0;
	/*END OF INITIALIZING VARIABLES RELATED TO VIRIAL, DISS POWER AND CALCULATION OF AMPLITUDES, PHASE*/
	
	/*INITIALIZATION OF VARIABLES*/
	time = 0.0;
	print_counter = 0;
	Energy = 0.0;
	
	PeakForce = 0.0;
	MaxTipPos = 0.0;
	MinTipPos = 0.0;
	
	Virial_compu = 0.0;
	Av_Force_compu = 0.0;
	Av_Force_an = 0.0;
	/*INITIALIZATION OF VARIABLES*/
	
	/*INITIALIZING VARIABLES TO GET NUMERICAL FORCE*/
	Fts = 0.0;
	Xb = 0.0;
	Zs = 0.0;
	sum_G_xb_xc = 0.0;
	sum_G_xc = 0.0;
	for (i=0; i < (n) ; i = i+1)
	{
		Xc[i] = 0.0;
		Xc_rate[i] = 0.0;
		G_xb_xc[i] = 0.0;
		G_xc[i] = 0.0;
	}
	/*INITIALIZING VARIABLES TO GET NUMERICAL FORCE*/
	
	/*CALCULATION OF THE TIP-SAMPLE FORCES*/  
	do  //advancing in time
	{
		time = time + dt;
		
		/*INTEGRATION OF THE EOM THROUGH THE VERLET ALGORITHM  */
        // Calculating the acceleration for the three eigenmodes using the equation of motion of a tip excited cantilever
        a1 = ( -k_m1*z1 - (mass*(f1*2*pi)*v1/Q1) + F_drive1*cos((f1*2*pi)*time) + F_drive2*cos((f2*2*pi)*time) + F_drive3*cos((f3*2*pi)*time) + Fts ) /mass;
        
        a2 = ( -k_m2*z2 - (mass*(f2*2*pi)*v2/Q2) + F_drive1*cos((f1*2*pi)*time) + F_drive2*cos((f2*2*pi)*time) + F_drive3*cos((f3*2*pi)*time) + Fts ) /mass;
        
        a3 = ( -k_m3*z3 - (mass*(f3*2*pi)*v3/Q3) + F_drive1*cos((f1*2*pi)*time) + F_drive2*cos((f2*2*pi)*time) + F_drive3*cos((f3*2*pi)*time) + Fts ) /mass;
         
        // Verlet algorithm to calculate velocity and position of the tip
        z1_new = 2*z1 - z1_old + a1*pow(dt, 2);
        z2_new = 2*z2 - z2_old + a2*pow(dt, 2);
        z3_new = 2*z3 - z3_old + a3*pow(dt, 2);
        
        v1 = (z1_new - z1_old)/(2*dt);
        v2 = (z2_new - z2_old)/(2*dt);
        v3 = (z3_new - z3_old)/(2*dt);
             
        // Updating z1_old and z1 for the next run
        z1_old = z1;
        z1 = z1_new;
        
        z2_old = z2;
        z2 = z2_new;
        
        z3_old = z3;
        z3 = z3_new;
    
        // Calculating tip position and velocity
        TipPos_old = TipPos;         //storing the previous value of tipPos to make the integration of the dissipation loop
        TipPos = z1_new + z2_new + z3_new + z_base;
		TipV = v1 + v2 + v3;
        /*END OF INTEGRATION OF EOM WITH VERLET ALGORITHM*/
		
		
		
        /*CALCULATION OF THE TIP-SAMPLE FORCES*/  		
		if (TipPos > Xb)
		{   
			for (i=0; i < (n) ; i=i+1) 
				{
					G_xc[i] = G[i]*Xc[i];
					sum_G_xc = sum_G_xc + G_xc[i];
				}
			if (  (sum_G_xc)/(Gg) > TipPos  )
			{ //contact, the surface surpassed the sample surface
				Xb = TipPos;
				for (i=0; i < (n) ; i=i+1) 
				{
					G_xb_xc[i] = G[i]*(Xb-Xc[i]);
					sum_G_xb_xc = sum_G_xb_xc + G_xb_xc[i];
					Xc_rate[i] = G[i]*(Xb-Xc[i])/eta[i];
					Xc[i] = Xc[i] + Xc_rate[i]*dt;
				}
				Fts = alfa*(-Ge*Xb - sum_G_xb_xc)- H*R/(6*pow(a,2));							// - H*R/(6*pow(a,2)) of force is added from the vdW interaction
			}
			else
			{ //true non contact
				Xb = (sum_G_xc)/(Gg);
				Fts = -H*R/(6*pow((  (TipPos - Xb) + a  ),2) );	// only vdw interaction
				for (i=0; i < (n) ; i=i+1) 
				{
					Xc_rate[i] = G[i]*(Xb-Xc[i])/eta[i];
					Xc[i] = Xc[i] + Xc_rate[i]*dt;
				}
			}
			sum_G_xc = 0.0;
			sum_G_xb_xc = 0.0;
		}
		
		else
		{ //Contact region TipPos is lower than the base
			if (TipPos_old > Xb)
				Zs = Xb;
			Xb = TipPos;
			for (i=0; i < (n) ; i++) 
			{
				G_xb_xc[i] = G[i]*(Xb-Xc[i]);
				sum_G_xb_xc = sum_G_xb_xc + G_xb_xc[i];
				Xc_rate[i] = G[i]*(Xb-Xc[i])/eta[i];
				Xc[i] = Xc[i] + Xc_rate[i]*dt;
			}
			Fts = alfa*(-Ge*Xb - sum_G_xb_xc)- H*R/(6*pow(a,2));							// - H*R/(6*pow(a,2)) of force is added from the vdW interaction
			sum_G_xb_xc = 0.0;
		}
		/*END OF CALCULATION OF THE TIP-SAMPLE FORCES*/ 
		
		/*NUMERICAL CALCULATION OF PEAK FORCE*/
		if (time>startprintF)
		{
			if (Fts > PeakForce)
				PeakForce = Fts;
		}
		/*NUMERICAL CALCULATION OF PEAK FORCE*/
		
		/*NUMERICAL CALCULATION OF TOTAL DISSIPATED ENERGY AND VIRIAL*/
		if (time > startprintF)
		{
		    Energy = -Fts*TipV*dt + Energy;    //calculating dissipated energy
			
			Virial_compu = Fts*((z1_new + z2_new + z3_new))*dt + Virial_compu;
			Av_Force_compu = Fts*dt + Av_Force_compu;
		}
		
		
	   	/*NUMERICAL CALCULATION OF TOTAL DISSIPATED ENERGY AND VIRIAL*/ 
					
		/*CALCULATION OF THE PHASE AND AMPLITUDE*/
		if (time > 0.9*startCalcAmp)
		{
			SumA = SumA+cos(f1*2*pi*time)*z1_new*dt;
			SumB = SumB+sin(f1*2*pi*time)*z1_new*dt;
					
			if (time > counterPh1*NN*period1)
			{
			  Output_Phase1 = atan(SumB/SumA)*180/pi;
			  if (Output_Phase1 < 0)
				Output_Phase1 = Output_Phase1 + 180;
			  Output_Amp1 = sqrt(SumA*SumA+SumB*SumB)*f1*2/NN;
				  
			  counterPh1 = counterPh1 + 1;
				 
			  SumA = 0;
			  SumB = 0;
			}
		
		   
		}  
		/*END OF CALCULATION OF PHASE AND AMPLITUDE*/
		
		/*CALCULATION OF POWER DISSIPATED AND VIRIAL*/
		Ets1 = (pi*k_m1*Output_Amp1*Output_Amp1/Q1)*( (A1/Output_Amp1)*sin(Output_Phase1*pi/180)-1 );  //Energy dissipated per fundamental period. From Tamayo and Garcia 1998 APL, Relationship between phase shift and energy dissipation in tapping-mode scanning force microscopy
		Vts1 = (k_m1*Output_Amp1/2)*(-A1*cos(Output_Phase1*pi/180)/Q1);
		Av_Force_an = Vts1/-Output_Amp1;
		/*END OF CALCULATION OF THE DISSIPATED POWER AND VIRIAL*/
		
			
		/*PRINTING TIME NUMERICAL RESULTS*/
		if (time>startprintF)
		{
			print_counter = print_counter+1;
			if (print_counter > NP1)
			{
				fprintf(Out1, "%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\n", time*1e6, TipPos*1e9, Fts*1e9, Xb*1e9, TipV*1e9);
				print_counter = 0;
			}
		}
        /*END OF PRINTING TIME NUMERICAL RESULTS*/

	
	} while (time < (simultime)); // closing do loop

    
	
	/*CALCULATION OF NUMERICAL VIRIAL, DISSIPATED ENERGY, AND AVERAGE TIP-SAMPLE FORCE*/
	Virial_compu = Virial_compu/(simultime-startprintF);
	Av_Force_compu = Av_Force_compu/(simultime-startprintF);
	Energy = Energy/(simultime-startprintF)*period1;  //dissipated energy per fundamental period, using simulation results
	/*CALCULATION OF NUMERICAL VIRIAL AND AVERAGE TIP-SAMPLE FORCE*/

	fprintf(Out2, "%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\n", A1*1.0e9, k_m1, Q1, R*1.0e9, nu, H*1.0e19, z_base*1.0e9, Zs*1.0e9, Output_Amp1*1.0e9, Output_Phase1, f1/1.0e3, Energy*1.0e18, Virial_compu*1.0e18, Av_Force_compu*1.0e9, Ets1*1.0e18, Vts1*1.0e18, Av_Force_an*1.0e9);  //printing file with summarized results
	
    
	
	
	fclose(Out1);  //Closing output file with numerical results for each cantilever base position
    

} //Closing For loop that takes care of approaching the cantilever towards the sample	

fclose(Out2);  //Closing file with average quantities summarized

  
}//closing main program
            