// C++ code for LTS cells


#define  ICOUNT       10   // number of i-cells in the network 
#define  SVCOUNT2     200  // number of e-cells in the network 

  
// parameters for i-cells
	double Ena = 50;  // all reversal potentials in mV
	double Ek  = -100;
	double El  = -67;  
	double Ei = -80; 
	double igna = 100; // all conductances in mS/cm^2
	double igk = 80;
	double gmi = 2;
	double igl = 0.1;
	double Q10 = 3.209;  // for the m-current
	double vhalf = -30;  // units of mV

// heterogeneity added to the applied current to the i-cells
	double Iappi[ICOUNT];
	Iappi[0] = 2.0;
	for (z=1; z<ICOUNT; z++)
		Iappi[z] = Iappi[0] + 0.005; 

// parameters for GABAa synapse
	 gi_i = 0.25; // mS/cm^2
	 taui = 10.0;   // msec

// equations for inhibitory GABAa synaptic activity (i to i kinetics)
	double sgi[SVCOUNT2];
	sgi[z] = 2*(1+tanh(y[V][z]/4)); // activation variable for i to i (or e to i kinetics)
   	k[ 97 ][z] = sgi[z]*(1-y[97][z])-y[97][z]/taui;  

	double Sum97 = 0;	
	for (z=0; z<ICOUNT; z++)
		Sum97 = Sum97 + y[97][z];


// derivative equation for i-cells
	k[ V ][z] = -igna*y[H][z]*pow(y[M][z],3.0)*(y[V][z]-Ena)
		    -igk*pow(y[N][z],4.0)*(y[V][z]-Ek)-igl*(y[V][z]-El)
		    -gmi*y[N2][z]*(y[V][z]-Ek) 
		    -gi_i/(ICOUNT-1.0)*(Sum97-y[97][z])*(y[V][z]-Ei)
		     + Iappi[z];     

// derivative equations for the gating variables
   	k[ M ][z] = 0.32*(y[V][z]+54)/(1-exp(-(y[V][z]+54)/4))*(1-y[M][z])-0.28*(y[V][z]+27)/(exp((y[V][z]+27)/5)-1)*y[M][z];     // derivative equation for sodium activation variable
   	k[ H ][z] = 0.128*exp(-(y[V][z]+50)/18)*(1-y[H][z])-4/(1+exp(-(y[V][z]+27)/5))*y[H][z];   // derivative equation for sodium inactivation variable
   	k[ N ][z] = 0.032*(y[V][z]+52)/(1-exp(-(y[V][z]+52)/5))*(1-y[N][z])-0.5*exp(-(y[V][z]+57)/40)*y[N][z];    // derivative equation for fast potassium activation variable
   	k[ N2 ][z] = (Q10*0.0001*(y[V][z]-vhalf)/(1-exp(-(y[V][z]-vhalf)/9)))*(1-y[N2][z]) + Q10*0.0001*(y[V][z]-vhalf)/(1-exp((y[V][z]-vhalf)/9))*y[N2][z];  // derivative equation for for m-current activation variable -- Mainen and Sejowski



