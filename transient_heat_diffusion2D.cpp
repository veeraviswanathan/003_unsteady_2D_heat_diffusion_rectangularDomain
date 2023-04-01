#include<stdio.h>
#include<stdlib.h>
#include<array>
#include<cmath>


double Flux_calc(double flux_in, double flux_out,double dx);

int main()
{
	//declare domain variables
	double end_time=300 ; // seconds
	double time_now=0.0;
	double dt=1e-5; // time step
printf(" DEbug 0: \n");
	
	int n_nodes_x = 20; // no of grid points
	int n_nodes_y = 20; // no of grid points
	double length_x = 0.05;// length of domain (m)
	double length_y = 0.05;// length of domain (m)
	double delta_x = (double)length_x / n_nodes_x;
	double delta_y = (double)length_y / n_nodes_y;
	
	double Tl = 600, Tr = 273; //left and right boundary condition (K)
	double Tinit = 300; // initial temperature for domain
	double k = 2e-3; // Diffusivity (m2/s)

	double** T = nullptr; // Temperature array
	double** Told = nullptr; // Temperature array old
	double** S = nullptr; // Source array		
	double rhs_coefficient=0.0;
	double Fr, Fl;
	double FE,FW,FN,FS;
	FILE *fp;

	// declare index variables
	int i,j,indP;
	int iter = 0;
	char fname[80];
	// fp=fopen("user_heat_diffusion_transient_explicit.out","w");
	// fprintf(fp,"%s \t ","Time");
	// for (int kk=0;kk<n_nodes;kk++)
	// {
	// 	fprintf(fp,"%s[%d]\t","T",kk);
	// }
	// fprintf(fp,"\n");

	// fclose(fp);

	printf(" DEbug 1: \n");
	//memory allocation
	T 		= (double**)malloc(n_nodes_x*sizeof(double*));
	Told 	= (double**)malloc(n_nodes_x*sizeof(double*));
	S 		= (double**)malloc(n_nodes_x*sizeof(double*));
	for(j=0;j<n_nodes_x;j++)
	{
		T[j] 		= (double*)malloc(n_nodes_y*sizeof(double));
		Told[j] 	= (double*)malloc(n_nodes_y*sizeof(double));
		S[j] 		= (double*)malloc(n_nodes_y*sizeof(double));

	}

	printf(" DEbug 2: \n");


	// initialisation block

	rhs_coefficient=k*dt;
	for (i = 0; i< n_nodes_x;i++)
	{
		for (j = 0; j<n_nodes_y;j++)
		{
			T[i][j] = Tinit;
			Told[i][j] = T[i][j];
			S[i][j] = 0 * (dt);
		}

	}
	
	//S[int(n_nodes/2)]= 5000 * (delta_x * delta_x / k);	
	
	//Gauss Siedel iterative method
	// fp=fopen("user_heat_diffusion_transient_explicit.out","a+");
	while(time_now<end_time)
	{
		
		printf("Iteration: %d\n",iter);
		// Assign old array
		for (i = 0;i < n_nodes_x;i++)
		{
			for (j = 0;j < n_nodes_x;j++)
			{		
				Told[i][j] = T[i][j];
			}			
		}

				
		// Boundary conditions for the rectangular domain; LEFT WALL= Tl, RIGHT wall=TR; TOP and bottom wall , Flux=0;
		//TOP and bottom
		for(i=0;i<n_nodes_x;i++)
		{
			FN=Flux_calc(Told[i][0],Told[i][1],delta_y);
			FS=0.0;
			if(i==0)
			{
				FW=Flux_calc(2*Tl,2*Told[i][0],delta_x);
				// FW=0.0;
			}
			else
			{
				FW=Flux_calc(Told[i-1][0],Told[i][0],delta_x);
			}
			if(i==n_nodes_x-1)
			{
				// FE=Flux_calc(2*Told[i][0],2*Tr,delta_x);
				FE=0.0;
			}
			else
			{
				FE=Flux_calc(Told[i][0],Told[i+1][0],delta_x);
			}			
			
			// BOTTOM WALL TEMPERATURE			
			T[i][0] = Told[i][0] + (rhs_coefficient/delta_x)*(FW-FE)+(rhs_coefficient/delta_y)*(FS-FN)+dt*S[i][0];
			
			// #################################//
			FN=0.0;
			FS=Flux_calc(Told[i][n_nodes_y-2],Told[i][n_nodes_y-1],delta_y);				;
			if(i==0)
			{
				FW=Flux_calc(2*Tl,2*Told[i][n_nodes_y-1],delta_x);
			}
			else
			{
				FW=Flux_calc(Told[i-1][n_nodes_y-1],Told[i][n_nodes_y-1],delta_x);
			}
			if(i==n_nodes_x-1)
			{
				// FE=Flux_calc(2*Told[i][n_nodes_y-1],2*Tr,delta_x);
				FE=0.0;
			}
			else
			{
				FE=Flux_calc(Told[i][n_nodes_y-1],Told[i+1][n_nodes_y-1],delta_x);
			}	
			// TOP WALL TEMPERATURE		
			T[i][n_nodes_y-1] = Told[i][n_nodes_y-1] + (rhs_coefficient/delta_x)*(FW-FE)+(rhs_coefficient/delta_y)*(FS-FN)+dt*S[i][n_nodes_y-1];
			
		}
		//LEFT and RIGHT
		for(j=0;j<n_nodes_y;j++)
		{
			FW=Flux_calc(2*Tl,2*Told[0][j],delta_x);			
			FE=Flux_calc(Told[0][j],Told[1][j],delta_x);			
			if(j==0)
			{
				FS=0.0;
			}
			else
			{
				FS=Flux_calc(Told[0][j-1],Told[0][j],delta_y);
			}
			if(j==n_nodes_y-1)
			{
				FN=0.0;
			}
			else
			{
				FN=Flux_calc(Told[0][j],Told[0][j+1],delta_y);
			}			
			
			// LEFT WALL TEMPERATURE
			T[0][j] = Told[0][j] + (rhs_coefficient/delta_x)*(FW-FE)+(rhs_coefficient/delta_y)*(FS-FN)+dt*S[0][j];
			// printf("LEFT: FN=%e FS=%e FE=%e FW=%e\n Temperature=%e\n",FN,FS,FE,FW,T[0][j]);

			// #################################//
			// FE=Flux_calc(2*Told[n_nodes_x-1][j],2*Tr,delta_x);	
			FE=0.0;		
			FW=Flux_calc(Told[n_nodes_x-2][j],Told[n_nodes_x-1][j],delta_x);				;
			if(j==0)
			{
				FS=0.0;
			}
			else
			{
				FS=Flux_calc(Told[n_nodes_x-1][j-1],Told[n_nodes_x-1][j],delta_y);
			}
			if(j==n_nodes_y-1)
			{
				FN=0.0;
			}
			else
			{
				FN=Flux_calc(Told[n_nodes_x-1][j],Told[n_nodes_x-1][j+1],delta_y);
			}			
			// RIGHT WALL TEMPERATURE		
			T[n_nodes_x-1][j] = Told[n_nodes_x-1][j] + (rhs_coefficient/delta_x)*(FW-FE)+(rhs_coefficient/delta_y)*(FS-FN)+dt*S[n_nodes_x-1][j];
			
		}
		// Main Stencil
		for (i = 1; i < n_nodes_x-1; i++)
		{
			for(j=1; j<n_nodes_y-1;j++)
			{				
				FE=Flux_calc(Told[i][j],Told[i+1][j],delta_x);
				FW=Flux_calc(Told[i-1][j],Told[i][j],delta_x);
				FN=Flux_calc(Told[i][j],Told[i][j+1],delta_y);
				FS=Flux_calc(Told[i][j-1],Told[i][j],delta_y);				
				T[i][j] = Told[i][j] + (rhs_coefficient/delta_x)*(FW-FE)+(rhs_coefficient/delta_y)*(FS-FN)+dt*S[i][j];
			}
		}
		
		printf("Iteration: %d\n",iter);
		//Write output file
		//if(fmod(time_now,2.0)<5e-3)
		
   		sprintf(fname, "user_heat_diffusion_transient_explicit_%f.out", time_now);	
		
		if(iter%1000==0)
		{
			fp=fopen(fname,"w");
			for(j=0;j<n_nodes_y;j++)
			{
				for (i = 0;i < n_nodes_x;i++)
				{
					fprintf(fp,"%e \t",T[i][j]);
				}
				fprintf(fp,"\n");
			}
			fclose(fp);	
		}
		
		iter++;
		time_now=time_now+iter*dt;
	}
	
	free(T);
	free(Told);
	free(S);
	return(0);
}

double Flux_calc(double flux_in, double flux_out, double dx)
{
	double flux;
	flux= (flux_out-flux_in)/dx;
    return -1*flux;
}
