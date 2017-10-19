/*THE PURPOSE OF THIS PROGRAM IS TO STUDY STRESS RELAXATION OVER A GENERALIZED VISCOELASTIC SURFACE USING PRONY SERIES*/

/*HEADERS*/
#include<stdio.h>
#include<math.h>
/*HEADERS*/

/* CONSTANTS */
#define 	pi	3.1415926535
#define     n  1               //Number of arms in the Generalized Maxwell model
#define     Nprint  50000      //Number of data points printed
#define     NN  50      //number of fundamental periods used to calculate phase and amplitude
/* CONSTANTS */

/*DECLARATION OF VARIABLES*/
//variables used to write txt outputfiles 
FILE		*Out1, *Out2, *Out3;						// Pointer for output files
char		output_file_1[100], output_file_2[100], output_file_3[100];	
char        output_compu[100];						// Output file root name
char        output_analytical[100];                 // Output file root name

double		G[n];  //array storing the values of Modulus of each Maxwell arm
double      tau[n];  //array storing the values of Relaxation times of each Maxwell arm
double		eta[n];  //array storing the values of viscosity of the dashpots of each of the Maxwell arms
double      Xc[n], Xc_rate[n];
double      G_xb_xc[n], G_xc[n];
double 		Gg, Ge, G_loss, G_storage, G_absolute;    //Glassy modulus, equilibrium, loss, storage, and absolute moduli, respectively
double      alfa, R, nu;                  //cell constant, radius of tip, Poisson's ratio of sample
double      Z, A, omega, f1, period1, dt, simultime;  
double      sum_transient_E, E_analytical, E_storage, E_loss, E_equilibrium, E_transient;
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
double      startprintF, stopprintF, startCalcAmp;            // Variables that control the print length of the ForceSLS file
int         i, print_counter;

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
double      V_storage, V_loss, V_transient, V_eq, V_analytical, sum_transient_V;
/*DECLARATION OF VARIABLES*/


/*BEGINING OF MAIN PROGRAM*/
main()
{
/*INSERTING VALUE OF VARIABLES RELATED TO DYNAMICS OF CANTILEVER*/
k_m1 = 3000.0;                 //stiffness of the first eigenmode
Q1 = 150;
Q2 = 450;
Q3 = 750;

A1 = 10;                    //target free amplitude of the first eigenmode
A1 = A1*1e-9;	     	    // Converting to meters
A2 = 0;
A2 = A2*1e-9;
A3 = 0;                     //target free amplitude of the third eigenmode
A3 = A3*1e-9;               // Converting to meters 
   
z_base_max = A1*1e9;        //cantilever base maximum position in nanometers
/*INSERTING VALUE OF VARIABLES RELATED TO DYNAMICS OF CANTILEVER*/

/*PARAMETERS CALCULATED FOR THE SECOND AND THIRD MODES FROM THE DATA INSERTED BY USER*/
//calculation of frequency of the three eigenmodes
f1 = 6.0e5;
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
startCalcAmp = period1*3.0*Q1;
startprintF = period1*3.0*Q1 + NN*period1;
stopprintF = startprintF + NN*period1;


//Defining run length in seconds,
simultime = stopprintF;

//Defining timestep
dt = period1/1.0e4;
/*CALCULATION FOR THE SECOND AND THIRD MODES ENDS HERE*/  

/*PARAMETERS RELATED TO AFM DYNAMICS when given prescribed trajectory*/
omega = 2.0*pi*f1;
NP1 = ((stopprintF-startprintF)/dt)*(1.0/Nprint);  //variable used for counting the print steps
/*PARAMETERS RELATED TO AFM DYNAMICS when given prescribed trajectory*/

/*MODEL PARAMETERS*/
//Ge = 1.0e4; //Equilibrium modulus in Pa
R = 10.0e-9; //cylindrical flat punch end diameter
nu = 0.5;  //material's Poisson ratio
alfa = 4.0*R/(1.0-nu);  //Cell constant converting stress/strain to force/displacement. This is derived form application of correspondence principle to elastic solution obtained by Sneddon 1965 (Int Journal of Engng Sci. vol3, pp.47-57), see eq 6.1
//alfa = (2.0*R)/(1.0-pow(nu,2));//temporal test
Ge = 0.2e9;  //temporal test
Gg = Ge; //initializing it to value of equilibrium modulus

G[0] = 2.0e9;
 
tau[0] = 2.0e-8; 

/*MODEL PARAMETERS*/


/*DEFINING VALUES OF VISCOSITY FOR EACH DASHPOT, CALCULATING Gg, Storage and Loss moduli */
for (i=0; i < (n) ; i = i+1)
{
	eta[i] = tau[i] * G[i];
	Gg = G[i]+Gg;
	sum_Gstorage = sum_Gstorage + (G[i]*pow(omega,2)*pow(tau[i],2))/( 1.0+pow(omega,2)*pow(tau[i],2) );  //see eq. (3.5-9)1 Nicholas Tschoegl, Phenomenological theory of linear viscoelastic behavior
	G_loss = G_loss + (G[i]*omega*tau[i])/( 1.0+pow(omega,2)*pow(tau[i],2) ); //see eq. (3.5-10) Nicholas Tschoegl, Phenomenological theory of linear viscoelastic behavior
}

G_storage = Ge + sum_Gstorage;   //see eq. (3.5-9)1 Nicholas Tschoegl, Phenomenological theory of linear viscoelastic behavior
G_absolute = sqrt(pow(G_loss, 2)+pow(G_storage,2));
/*DEFINING VALUES OF VISCOSITY FOR EACH DASHPOT, CALCULATING Gg*/

//opening output file 2
sprintf(output_file_2,"Summary.txt");
Out2 = fopen(output_file_2,"w");
  
// Printing the output file header of output file 2
fprintf(Out2,"Diss_Energy_compu(aJ)\tDiss_Energy_Cleveland(aJ)\tz_base\tAmp_Afree\tPhase1\tDiss_Analytical(aJ)\tLoss_Diss_analytic(aJ)\tStore_Diss_analytic(aJ)\tGe_Diss_analytic(aJ)\tTransient_diss_analytic(aJ)\tPeakForce(nN)\tZ_sample(nm)\tIndentation(nm)\tVirial_compu(aJ)\tVirial_SPaulo(aJ)\tAvForce_compu(nN)\tAvForce_an(nN)\tt_biprime_an(us)\tt_biprime_semian(us)\tVirial_an(aJ)\tVirial_store(aJ)\tVirial_loss(aJ)\tVirial_Ge(aJ)\tVirial_transient(aJ)\n");  


/*START OF FOR LOOP WHERE THERE IS GOING TO BE A CHANGE OF THE AMPLITUDE OF THE SECOND EIGENMODE*/ 
//for (z_base_counter = (z_base_max-10.0); z_base_counter > 0.1; z_base_counter = z_base_counter - 10.0)
for (z_base_counter = (z_base_max-6.0); z_base_counter > 0.1; z_base_counter = z_base_counter - 100.0)
{
	z_base = z_base_counter*1e-9;
	
	/*WRITING NAME, OPENING, AND WRITING HEADER OF OUTPUT FILE WITH NUMERICAL RESULTS*/
	// Changing the name of the output file1
    sprintf(output_file_1,"%scompu_%.2f.txt",output_compu,z_base*1e9);
  
    // Opening the output file1
    Out1 = fopen(output_file_1,"w"); 

    // Printing the output file header of output file 1
    fprintf(Out1,"Time(us)\tTipPos(nm)\tCosRef(nm)\tForce(nN)\tSamplePos(nm)\n");
    /*END OF WRITING NAME, OPENING, AND WRITING HEADER OF OUTPUT FILE WITH NUMERICAL RESULTS*/
	
	/*WRITING NAME, OPENING, AND WRITING HEADER OF OUTPUT FILE WITH ANALYTICAL RESULTS*/
	//opening output file 3
	sprintf(output_file_3,"%sanalytical_%.2f.txt",output_analytical,z_base*1e9);
	Out3 = fopen(output_file_3,"w");
  
	// Printing the output file header of output file 3
	fprintf(Out3,"Time_an(us)\tTip_Prescribed(nm)\tF_Analytical(nN)\tF_DMA_F_Ge(nN)\tF_Transient(nN)\n");  
	/*WRITING NAME, OPENING, AND WRITING HEADER OF OUTPUT FILE WITH ANALYTICAL RESULTS*/
				
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
		
		/*CALCULATING MAX AND MIN TIP POS*/
		if (time>startprintF)
		{
			if (TipPos>MaxTipPos)
				MaxTipPos = TipPos;
			if (TipPos<MinTipPos)
				MinTipPos = TipPos;
		}
		/*CALCULATING MAX AND MIN TIP POS*/
		
        /*CALCULATION OF THE TIP-SAMPLE FORCES*/  		
		if (TipPos > Xb)   //apparent non contact
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
				Fts = alfa*(-Ge*Xb - sum_G_xb_xc); //- 7.5e-9;							// -7.5e-9 of force is added from the vdW interaction
			}
			else
			{ //true non contact
				Xb = (sum_G_xc)/(Gg);
				Fts = 0.0; //-7.5e-9*0.2e-9*0.2e-9/((TipPos - Xb)+0.2e-9)/((TipPos - Xb)+0.2e-9);		// only vdw interaction 
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
			Fts = alfa*(-Ge*Xb - sum_G_xb_xc); // - 7.5e-9;							// -7.5e-9 of force is added from the vdW interaction
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
			Virial_compu = Fts*(z1_new + z2_new + z3_new)*dt + Virial_compu;
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
		
		if (time>startprintF)
			CosRef = Output_Amp1*cos(omega*time - Output_Phase1*pi/180) - Z;
		
		/*PRINTING TIME NUMERICAL RESULTS*/
		if (time>startprintF)
		{
			print_counter = print_counter+1;
			if (print_counter > NP1)
			{
				fprintf(Out1, "%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\n", time*1e6, TipPos*1e9, CosRef*1e9, Fts*1e9, Xb*1e9);
				print_counter = 0;
			}
		}
        /*END OF PRINTING TIME NUMERICAL RESULTS*/

	
	} while (time < (simultime)); // closing do loop

    
	/*CALCULATION OF VARIABLES FOR ANALYTICAL EXPRESSIONS*/
	Z = -z_base + Zs;     //Sample position, this value will be used for the analytical expressions
	A = Output_Amp1;
	sum_transient_E =0.0;
	
	t_in = (asin(-Z/A)+pi)*(1.0/omega); 
	t_fin = 1.0/omega*(  asin(Ge*(Z/A)/G_absolute) + atan(G_storage/G_loss) + 3.0/2*pi );
	t_fin_2 = t_fin;
	/*CALCULATION OF VARIABLES FOR ANALYTICAL EXPRESSIONS*/
	
	/*ANALYTICAL CALCULATION OF TIP SAMPLE FORCES*/
	time = 0.0; //initializing time for the calculation of analytical Tip-Sample Force
	print_counter = 0;
	/*Initializing analytical forces*/
	F_analytical = 0.0;
	E_analytical = 0.0;
	sum_transient = 0.0;
	sum_transient_E = 0.0;
	V_analytical = 0.0;
	sum_transient_V = 0.0;
	
	E_storage = 0.0;
	E_equilibrium = 0.0;
	E_loss = 0.0;
	E_transient = 0.0;
	
	V_storage = 0.0;
	V_eq = 0.0;
	V_loss = 0.0;
	V_transient = 0.0;
	/*Initializing analytical forces*/
		
	do
	{
		time = time +dt;
		TipPrescribed = A*sin(omega*time) + z_base;
		sum_transient = 0.0; //initializing time for the calculation of analytical Tip-Sample Force
		
			for (i=0; i<(n); i++)
			{
			sum_transient = sum_transient + G[i]*tau[i]/(1.0+pow(omega,2)*pow(tau[i],2))*exp(-(time-t_in)/tau[i])*( pow(omega,2)*Z*tau[i]-omega*sqrt(pow(A,2)-pow(Z,2)) );  
			}
			
			if (time > t_in)
			{
				if(time > t_fin_2)
				{
					F_analytical = 0.0;
				}
				else
				{
					F_DMA = -alfa*(  G_storage*A*sin(omega*time) + G_loss*A*cos(omega*time) );
					F_Ge = -alfa*(-Ge*Z);
					F_transient = -alfa*(-sum_transient);
					F_analytical = F_DMA + F_Ge + F_transient;
					if (F_analytical<0.0)
						t_fin_2 = time;
				}
			}
	
		print_counter = print_counter+1;
		
		fprintf(Out3, "%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\n", time*1e6, TipPrescribed*1e9, F_analytical*1e9, F_DMA+F_Ge*1e9, F_transient*1e9);
		print_counter = 0;
				
	 
	} while (time<period1);
	/*END OF CALCULATION OF ANALYTICAL FORM OF TIP SAMPLE FORCES*/

	/*ANALYTICAL CALCULATION OF DISSIPATED ENERGY*/
	for (i=0; i<(n); i++)
		{
		sum_transient_E = sum_transient_E + G[i]*pow(tau[i],2)/pow((1.0+pow(omega,2)*pow(tau[i],2)),2)*(tau[i]*A*omega*Z-A*sqrt(pow(A,2)-pow(Z,2)))
							*(  exp(-(t_fin_2-t_in)/tau[i])*( tau[i]*omega*sin(omega*t_fin_2)-cos(omega*t_fin_2) ) - tau[i]*omega*Z/A- 1.0/A*sqrt(pow(A,2)-pow(Z,2))  );
		}
	E_storage = alfa*(  G_storage/2.0*pow(Z,2)* (pow((A/Z*sin(omega*t_fin_2)),2)-1)   );
	E_loss = alfa*( G_loss/4.0*( 2*pow(A,2)*omega*(t_fin_2-t_in) + 2*Z*sqrt(pow(A,2)-pow(Z,2)) + pow(A,2)*sin(2*omega*t_fin_2) )  );
	E_equilibrium = - alfa*( Ge*pow(Z,2)*(A/Z*sin(omega*t_fin_2)-1));
	E_transient = - alfa*( pow(omega,2)*sum_transient_E );
	E_analytical = E_storage + E_loss + E_equilibrium + E_transient;
	/*ANALYTICAL CALCULATION OF DISSIPATED ENERGY*/
	
	/*ANALYTICAL CALCULATION OF VIRIAL*/
	for (i=0; i<(n); i++)
		{
		sum_transient_V = sum_transient_V + G[i]*pow(tau[i],2)/(1.0+pow(omega,2)*pow(tau[i],2))*(exp(-(t_fin_2-t_in)/tau[i])-1)*( pow(omega,2)*Z*tau[i]-omega*sqrt(pow(A,2)-pow(Z,2)) );  
		}
	V_storage = - alfa*A/(2*pi)*( A*G_storage*(cos(omega*t_fin_2)+1.0/A*sqrt(pow(A,2)-pow(Z,2))) );
	V_loss = alfa*A/(2*pi)*(A*G_loss*(sin(omega*t_fin_2)-Z/A) );
	V_eq = - alfa*omega*A/(2*pi)*( Ge*Z*(t_fin_2-t_in) );
	V_transient = alfa*omega*A/(2*pi)*( sum_transient_V );
	V_analytical = V_storage + V_loss + V_eq + V_transient;
	/*ANALYTICAL CALCULATION OF VIRIAL*/
	
	/*CALCULATION OF NUMERICAL VIRIAL AND AVERAGE TIP-SAMPLE FORCE*/
	Virial_compu = Virial_compu/(stopprintF-startprintF);
	Av_Force_compu = Av_Force_compu/(stopprintF-startprintF);
	/*CALCULATION OF NUMERICAL VIRIAL AND AVERAGE TIP-SAMPLE FORCE*/
	

    fprintf(Out2, "%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\n", Energy*1e18/(NN), Ets1*1e18, z_base*1e9, Output_Amp1/A1, Output_Phase1, E_analytical*1e18, E_loss*1e18, E_storage*1e18, E_equilibrium*1e18, E_transient*1e18, PeakForce*1e9, Zs*1e9, (MinTipPos-Zs)*1e9, Virial_compu*1e18, Vts1*1e18, Av_Force_compu*1e9, Av_Force_an*1e9, t_fin*1e6, t_fin_2*1e6, V_analytical*1e18, V_storage*1e18, V_loss*1e18, V_eq*1e18, V_transient*1e18);
	
    
	
	
	fclose(Out1);  //Closing output file with numerical results for each cantilever base position
	fclose(Out3);  //Closing output file with analytical results for each cantilever base position
    

} //Closing For loop that takes care of approaching the cantilever towards the sample	

fclose(Out2);  //Closing file with average quantities summarized

printf("%8.8f\n",  Energy*1e18/(NN));
printf("%8.8f\n",  Ets1*1e18);
printf("%8.8f\n",  E_analytical*1e18);
printf("%8.8f\n",  Z*1e9);
printf("%8.8f\n",  A*1e9 );
printf("%8.8f\n",  MaxTipPos*1e9 );
printf("%8.8f\n",  MinTipPos*1e9 );


getch();
  
}//closing main program
            