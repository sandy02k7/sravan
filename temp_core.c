#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double min(double x, double y){
	return ((x < y) ? x : y);
}
//double K(double T){
//	double C = 0.22;
//	double Mn = 0.6;
//	double Ni = 0.11;
//	double Cu = 0.12;
//	double Cr = 0.453;
//	double Mo = 0.066;
//	double Ar3 = 910 - 273*C - 74*Mn - 56*Ni - 16*Cr - 9*Mo - 5*Cu;  
//	if	(T > Ar3) return 20.45 - 4.125*C + 0.0089*T;
//	else if (T < 769) return 41.971 + 0.018183*T;
//	else		  return 55.95 + 0.00248*(T - 769);
//}
int main()
{
	FILE *fptr;
        fptr = fopen("20mm_data(Q = 120).csv", "w"); 
	const double Cp = 600,
			 k  = 45,
			 dr = 0.0004,
		     R  = 0.010,
		     Rc = 0.0235,
		   Qtot = 120.0,
	      v_bar = 3.6,
		  no    = 1;
	double Qno   = Qtot/no;
	double p     = 2*3.14*((Rc + R)/2);
	double A	 = 3.14*((Rc*Rc) - (R*R));
	int    n     = (int)(R/dr);
	double time_taken = 0;
	double tot_time = 10;
	double T[n];
	for(int i = 0 ;i < n; i++) T[i] = 1000.0;
	double Tf[n];
	for(int i = 0 ;i < n; i++) Tf[i] = 1.0;
	double Dh = 4*A/p;
	printf("%.10f\n",A);
	double v_w = Qno/(3600*A);
	double v_rel = (v_w - v_bar) < 0 ? -1*(v_w - v_bar) : v_w - v_bar;
	double ro_w = 996;
	double mu = 0.000993;
	double Re = (v_rel*ro_w*Dh)/mu;
	double kw = 0.597;
	double pr = 7.0;
	double h = 0.0;
	double length_cl = 9.828;
	printf("%.10f\n",Re);
//	printf("Enter the value of the flow rate and  mill speed:");
//	scanf("%lf %lf",&Qno, &v_bar);

	if(Re < 5000) h = (1.86*pow((Re*pr),(0.33))*kw)/Dh;
	else          h = (0.023*pow(pr,(0.3))*pow(Re,(0.8)))/Dh;

	double ro = 7850;
	double abR = 2*11/(dr);
	double wbR = 2*h/(dr);
	double x = (dr*dr*ro)/2.0;
	double y = (dr*dr);
	double z = (R*dr);
	double T_water_air = 25.0;
	double time_thr = length_cl/v_bar;
	double dt = 9e9;
	printf("%.10f\n",h);
//	for(int i = 20; i <= 1200; i++)
//	{
	double dt0 = x*Cp/(k);
	double dtR = dt0/2.0;
	double dtr = 2*x*Cp/((2*k/y) - (k/z) + (wbR));
	dt = min(dt, min(dt0, min(dtr, dtR)));
//	}

	long itr = (long)(tot_time/dt);
	printf("%ld\n",itr);
	int samples = 100;
	double arr[samples][n];
	double time_taken_arr[samples];

	double vec[n];
	double a = k/y;
	double e = ro*Cp;
	double d = dt/e;
	double f = (1/d);
	double c0 = (f) - (4*a);

	double b = k/(2*z);
	double cr = (f) - 2*a;

	double bR = wbR;
    	double aR = 2*a - k/(z);
	double cR = 0;
	int is_thr = 1;
	for(long i = 0; i < itr; i++)
	{
	//------------T[0]-------------------------------------------
//		a = K(T[0])/y;
//		c0 = (f)-(4*a);
		Tf[0] = (4*a*T[1] + c0*T[0])*d;
		if(i%(itr/samples) == 0) vec[0] = Tf[0];
		
	//------------T[j]-------------------------------------------
		for(int j = 1; j < n; j++)
		{
//			a = K(T[j])/y;
//			b = K(T[j])/(2*z);
//			cr = (f) - (2*(K(T[j])/y));
			Tf[j] = (a*(T[j-1] + T[j+1]) + b*(T[j-1] - T[j+1]) + cr*(T[j]))*d;
			if(i%(itr/samples) == 0) vec[j] = Tf[j];
		}

	//------------T[n-1]-------------------------------------------
		
		if(time_taken > time_thr) bR = abR;
//		aR = 2*(K(T[n-1])/y) - K(T[n-1])/z;
		cR = (f) - (aR) - (bR);
		Tf[n-1] = (aR*T[n-2] + bR*T_water_air + cR*T[n-1])*d;
		if(i%(itr/samples) == 0) vec[n-1] = Tf[n-1];

		for(int j = 0; j < n; j++) T[j] = Tf[j];
		time_taken += dt;

		if(time_taken > time_thr && is_thr)
		{ 
			is_thr = 0;
			for(int j =0; j < n; j++)
			{
				fprintf(fptr, "%.13f,\n", Tf[j]);
			}
		}

		if(i%(itr/samples) == 0){
			for(int j = 0; j < n; j++) arr[i/(itr/samples)][j] = vec[j];
			for(int j = 0; j < n; j++) time_taken_arr[i/(itr/samples)] = time_taken;
			printf("%.13f %ld %.20f %ld %ld\n", Tf[n-1], i/(itr/samples), time_taken, i, itr/samples);
		}
	}

	for(int i = 0; i < samples; i++){
		for(int j = 0; j < n; j++){
			fprintf(fptr, "%.13f,", arr[i][j]);
		}
		fprintf(fptr, "%.20f\n",time_taken_arr[i]);
	} 
	fclose(fptr);

	return 0;

}
