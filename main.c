#include <ctype.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif


void update_pos(double *rx,double *ry,double vx,double vy,double fx,double fy,double m,int dt,int K2);
void update_vel(double* vx,double* vy,double fx1,double fy1,double fx2,double fy2,double m,int dt,int K2);
void calculate_force(double* fx,double* fy,double rx1,double ry1,double rx2,double ry2,double m1,double m2,double K1);
void particle_force(double* fx,double* fy,double rx1,double ry1,double rx2,double ry2,double rxp,double ryp,double m1,double m2,double mp,double K1);
void verl_pos(double* rx,double* ry,double vx,double vy,double fx,double fy,double m, int dt,int K2);
void verl_vel(double* vx,double* vy,double fx1,double fy1,double fx2,double fy2,double m, int dt,int K2);
void vec_diff(double* rx, double* ry, double rx1, double ry1,double rx2,double ry2);
void swap(double *a,double *b);

int main(void) /* int argc, char *argv[] */
{
	/* General Constants and simulation parameters (timestep and number of steps) */
	const double G = 6.67430 * pow(10,(-11));
	const double dt = (60*60*12);
	const double N = 2000;

	/*
	if (argc < 2)
    {
        printf("Usage: ./main config_file\n");
        return 1;
    }

    FILE *f = fopen(argv[1], "r");
    {
        if (f == NULL)
        {
            printf("Unable to open a file\n");
            fclose(f);
            return 1;
        }
    }
    */

	/* Reference Quantities */
	double m_r = 1.988500 * pow(10,30); /* Solar mass [kg] */
	double r_r = 149597870000; /* Astronomic Unit [m] */
	double t_r = (60*24); /* half a day */
	double v_r = 29780; /* avg. velocity of earth */

	/* Data */ 

	double m_e = 5.9724 * pow(10,24); /* Eath's mass */
	double T = (60*60*24*365.25); /* Orbital period */
	double w = 2*M_PI/(T); /* Angluar velocity around the Z-axis */

	/* Scaled constants */
	double K1 = G*t_r*m_r/(v_r*(pow(v_r,2)));
	double K2 = 1;

	/* Masses*/
	double m1 = m_r;
	double m2 = m_e;
	double mp = 100000;


	/* Initial position and velocities */
	double rx1_0 = 0;
	double ry1_0 = 0;

	double rx2_0 = r_r;
	double ry2_0 = 0;

	double rxp_0 = 0.5*r_r*((m_r-m_e)/(m_r+m_e));
	double ryp_0 = 0.5*sqrt(3)*r_r;
	double rp_norm = sqrt(pow(rxp_0,2)+pow(ryp_0,2));

	double vx1_0 = 0;
	double vy1_0 = 0;

	double vx2_0 = 1;
	double vy2_0 = v_r;

	double vxp_0 = -v_r*(ryp_0/rp_norm);
	double vyp_0 = v_r*(rxp_0/rp_norm);


	/* Openning the output file */
    FILE *out; 
	out = fopen("out.txt", "w");
	{
        if (out == NULL)
        {
            printf("Unable to create a file\n");
            fclose(out);
            return 1;
        }
    }



    double rx1 = rx1_0;
    double ry1 = ry1_0;
    double rx2 = rx2_0;
    double ry2 = ry2_0;
    double rxp = rxp_0;
    double ryp = ryp_0;
    double vx1 = vx1_0;
    double vy1 = vy1_0;
    double vx2 = vx2_0;
    double vy2 = vy2_0;
    double vxp = vxp_0;
    double vyp = vyp_0;

    int i = 0;
    while(i<N)
    {
    	i++;
    	fprintf(out,"%f,%f,%f,%f,%f,%f",rx1,ry1,rx2,ry2,rxp,ryp);
    	putc('\n', out);
    	double fx1,fy1,fx2,fy2,fxp1,fyp1,fxp2,fyp2;
    	calculate_force(&fx1,&fy1,rx1,ry1,rx2,ry2,m1,m2,G);
    	particle_force(&fxp1,&fyp1,rx1,ry1,rx2,ry2,rxp,ryp,m1,m2,mp,G);
    	update_pos(&rx1,&ry1,vx1,vy1,fx1,fy1,m1,dt,K2);
    	update_pos(&rx2,&ry2,vx2,vy2,-fx1,-fy1,m2,dt,K2);
    	update_pos(&rxp,&ryp,vxp,vyp,fxp1,fyp1,mp,dt,K2);
    	calculate_force(&fx2,&fy2,rx1,ry1,rx2,ry2,m1,m2,G);
    	particle_force(&fxp2,&fyp2,rx1,ry1,rx2,ry2,rxp,ryp,m1,m2,mp,G);
    	update_vel(&vx1,&vy1,fx1,fy1,fx2,fy2,m1,dt,K2);
    	update_vel(&vx2,&vy2,-fx1,-fy1,-fx2,-fy2,m2,dt,K2);
    	update_vel(&vxp,&vyp,fxp1,fyp1,fxp2,fyp2,mp,dt,K2);

    }

    fclose(out);
   	return(0);


}

void update_pos(double* rx,double* ry,double vx,double vy,double fx,double fy,double m,int dt,int K2)
{
	double rx_new = *rx;
	double ry_new = *ry;
	verl_pos(&rx_new,&ry_new,vx,vy,fx,fy,m,dt,K2);
	* rx = rx_new;
	* ry = ry_new;
}

void update_vel(double* vx,double* vy,double fx1,double fy1,double fx2,double fy2,double m,int dt,int K2)
{
	double vx_new = *vx;
	double vy_new = *vy;
	verl_vel(&vx_new,&vy_new,fx1,fy1,fx2,fy2,m,dt,K2);
	*vx = vx_new;
	*vy = vy_new;

}


void calculate_force(double* fx,double* fy,double rx1,double ry1,double rx2,double ry2,double m1,double m2,double K1)
{
	assert(fx);
	assert(fy);
	double rx,ry;
	vec_diff(&rx,&ry,rx1,ry1,rx2,ry2);
	double norm  = sqrt(pow(rx,2)+pow(ry,2));
	double force_mag = K1*m1*m2/(pow(norm,3));
	*fx = force_mag*rx;
	*fy = force_mag*ry;
}

void particle_force(double* fx,double* fy,double rx1,double ry1,double rx2,double ry2,double rxp,double ryp,double m1,double m2,double mp,double K1)
{
	assert(fx);
	assert(fy);
	double rx1p,ry1p;
	double rx2p,ry2p;
	vec_diff(&rx1p,&ry1p,rxp,ryp,rx1,ry1);
	vec_diff(&rx2p,&ry2p,rx2,ry2,rxp,ryp);
	double norm1 = sqrt(pow(rx1p,2)+pow(ry1p,2));
	double norm2 = sqrt(pow(rx2p,2)+pow(ry2p,2));
	double force_mag1 = K1*mp*m1/(pow(norm1,3));
	double force_mag2 = K1*mp*m2/(pow(norm2,3));
	*fx = force_mag1*rx1p+force_mag2*rx2p;
	*fy = force_mag1*ry1p+force_mag2*ry2p;
}

/* Velocity verlet function that takes values at time t as input and returns values at t+dt */
void verl_pos(double* rx,double* ry,double vx,double vy,double fx,double fy,double m, int dt,int K2)
{
	assert(rx);
	assert(ry);
	*rx = (*rx) + vx*dt*K2 + (fx*pow(dt,2))/(2*m);
	*ry = (*ry) + vy*dt*K2 + (fy*pow(dt,2))/(2*m);
}

void verl_vel(double* vx,double* vy,double fx1,double fy1,double fx2,double fy2,double m, int dt,int K2)
{
	assert(vx);
	assert(vy);
	*vx = *vx*K2 + ((fx1+fx2)*dt)/(2*m);
	*vy = *vy*K2 + ((fy1+fy2)*dt)/(2*m);
}

void vec_diff(double* rx, double* ry, double rx1, double ry1,double rx2,double ry2)
{
	assert(rx);
	assert(ry);
	*rx = rx2 - rx1;
	*ry = ry2 - ry1;
}

void swap(double *a,double *b)
{
  double temp = *a;
  *a = *b;
  *b = temp;
}