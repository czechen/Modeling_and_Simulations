#include <ctype.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif


void update_pos(double *rx,double *ry,double vx,double vy,double fx,double fy,double m,int dt);
void update_vel(double* vx,double* vy,double fx1,double fy1,double fx2,double fy2,double m,int dt);
void calculate_force(double* fx,double* fy,double rx1,double ry1,double rx2,double ry2,double m1,double m2,double G);
void particle_force(double* fx,double* fy,double rx1,double ry1,double rx2,double ry2,double rxp,double ryp,double m1,double m2,double mp,double G);
void verl_pos(double* rx,double* ry,double vx,double vy,double fx,double fy,double m, int dt);
void verl_vel(double* vx,double* vy,double fx1,double fy1,double fx2,double fy2,double m, int dt);
void rot_update_pos(double* rx,double* ry,double vx,double vy,double ax,double ay,double dt);
void rot_update_vel(double* vx,double* vy,double ax1,double ay1,double ax2,double ay2,double dt);
void rot_acc(double* ax,double* ay,double vx,double vy,double x,double y,double rho1,double rho2,double a);
void vec_diff(double* rx, double* ry, double rx1, double ry1,double rx2,double ry2);
void rotverl_pos(double* rx,double* ry,double vx,double vy,double ax,double ay, double dt);
void rotverl_vel(double* vx,double* vy,double ax1,double ay1,double ax2,double ay2, double dt);
double pseudo_pot(double* omega, double x, double y,double a);


int main(void) /* int argc, char *argv[] */
{
	/* General Constants and simulation parameters (timestep and number of steps) */
	const double G = 6.67430 * pow(10,(-11));
	const double dt1 = (60*60*12);
	const double N = 10000;

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
	double t_r = (60*60*24); /* a day */
	double v_r = 29780; /* avg. velocity of earth */

	/* Data */ 
	double m_e = 5.9724 * pow(10,24); /* Eath's mass */
	double m_m = 0.07346 * pow(10,24);
	double v_m = 970;
	double r_m = 0.3844 * pow(10,9);

	/* Masses*/
	double m1 = m_r;
	double m2 = m_e;
	double M = m1+m2;
	double mp = 100000;


	/* Initial position and velocities */
	double rx1_0 = 0;
	double ry1_0 = 0;

	double rx2_0 = r_r;
	double ry2_0 = 0;

	double rxp_0 = 0.5*r_r*((m1-m2)/M);
	double ryp_0 = 0.5*sqrt(3)*r_r;
	double rp_norm = sqrt(pow(rxp_0,2)+pow(ryp_0,2));

	double vx1_0 = 0;
	double vy1_0 = 0;

	double vx2_0 = 0;
	double vy2_0 = v_r;

	double vxp_0 = -v_r*(ryp_0/rp_norm);
	double vyp_0 = v_r*(rxp_0/rp_norm);


	/* Openning the output file */
    
    FILE *out1; 
	out1 = fopen("out1.txt", "w");
	{
        if (out1 == NULL)
        {
            printf("Unable to create a file\n");
            fclose(out1);
            return 1;
        }
    }

    /* Calculation in inertial reference frame */

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
    	fprintf(out1,"%f,%f,%f,%f,%f,%f",rx1,ry1,rx2,ry2,rxp,ryp);
    	putc('\n', out1);
    	double fx1,fy1,fx2,fy2,fxp1,fyp1,fxp2,fyp2;
    	calculate_force(&fx1,&fy1,rx1,ry1,rx2,ry2,m1,m2,G);
    	particle_force(&fxp1,&fyp1,rx1,ry1,rx2,ry2,rxp,ryp,m1,m2,mp,G);
    	update_pos(&rx1,&ry1,vx1,vy1,fx1,fy1,m1,dt1);
    	update_pos(&rx2,&ry2,vx2,vy2,-fx1,-fy1,m2,dt1);
    	update_pos(&rxp,&ryp,vxp,vyp,fxp1,fyp1,mp,dt1);
    	calculate_force(&fx2,&fy2,rx1,ry1,rx2,ry2,m1,m2,G);
    	particle_force(&fxp2,&fyp2,rx1,ry1,rx2,ry2,rxp,ryp,m1,m2,mp,G);
    	update_vel(&vx1,&vy1,fx1,fy1,fx2,fy2,m1,dt1);
    	update_vel(&vx2,&vy2,-fx1,-fy1,-fx2,-fy2,m2,dt1);
    	update_vel(&vxp,&vyp,fxp1,fyp1,fxp2,fyp2,mp,dt1);

    }

    fclose(out1);

    /* Co-Rotating reference frame */

    /* Data */
    double T = (60*60*24*365.25); /* Orbital period */
	double w = 2*M_PI/(T); /* Angluar velocity around the Z-axis */
	
	/* Characteristic units */
	double l_c = r_r; /* characteristic length */
	double t_c = sqrt(l_c/(G*M)); /* characteristic time */
	double dt2 = 500*t_c; 
	double a = m2/(M);


	/* Center of mass */
	double rxCM = (rx1_0*m1 + rx2_0*m2)/M;
	double ryCM = (ry1_0*m1 + ry2_0*m2)/M;

	/* Distanes from the center of mass to massive bodies */
	double ux1 = (rx1_0-rxCM)/l_c; 
	double uy1 = (ry1_0-ryCM)/l_c;
	double ux2 = (rx2_0-rxCM)/l_c;
	double uy2 = (ry2_0-ryCM)/l_c;
	double uxp = (rxp_0-rxCM)/l_c;
	double uyp = (ryp_0-ryCM)/l_c;
	
	/* Distance from the bodies to the particle and velocities of the particle */
	double rhox1 = (ux1-uxp);
	double rhoy1 = (uy1-uyp);
	double rhox2 = (ux2-uxp);
	double rhoy2 = (uy2-uyp);
	double rho1 = (sqrt(pow(rhox1,2)+pow(rhoy1,2)));
	double rho2 = (sqrt(pow(rhox2,2)+pow(rhoy2,2)));
	double up  = (sqrt(pow(uxp,2)+pow(uyp,2)));

	double vrott_xp = 0;
	double vrott_yp = 0;
	


    FILE *out2; 
	out2 = fopen("out2.txt", "w");
	{
        if (out2 == NULL)
        {
            printf("Unable to create a file\n");
            fclose(out2);
            return 1;
        }
    }

    fprintf(out2,"%f,%f,%f,%f",ux1,uy1,ux2,uy2);
    putc('\n', out2);

    int j =0;
    while (j<N)
    {
    	j++;
    	fprintf(out2,"%f,%f",uxp,uyp);
    	putc('\n', out2);
    	double ax1,ay1,ax2,ay2;
    	rot_acc(&ax1,&ay1,vrott_xp,vrott_yp,uxp,uyp,rho1,rho2,a);
    	rot_update_pos(&uxp,&uyp,vrott_xp,vrott_yp,ax1,ay1,dt2);
    	rot_acc(&ax2,&ay2,vrott_xp,vrott_yp,uxp,uyp,rho1,rho2,a);
    	rot_update_vel(&vrott_xp,&vrott_yp,ax1,ay1,ax2,ay2,dt2);

    }
    fclose(out2);

    FILE *out3; 
	out3 = fopen("out3.txt", "w");
	{
        if (out3 == NULL)
        {
            printf("Unable to create a file\n");
            fclose(out3);
            return 1;
        }
    }
    double x,y,omega;
    double h=0.06;
    for(x = -2;x < 2.1;x=x+h)
    {
    	for(y = -2;y < 2.1;y=y+h)
    	{
    		pseudo_pot(&omega,x,y,a);
    		fprintf(out3,"%f,%f,%f",-2*omega,x,y);
    		putc('\n', out3);
    	}
    }
    fclose(out3);

   	return(0);


}

void update_pos(double* rx,double* ry,double vx,double vy,double fx,double fy,double m,int dt)
{
	double rx_new = *rx;
	double ry_new = *ry;
	verl_pos(&rx_new,&ry_new,vx,vy,fx,fy,m,dt);
	* rx = rx_new;
	* ry = ry_new;
}

void update_vel(double* vx,double* vy,double fx1,double fy1,double fx2,double fy2,double m,int dt)
{
	double vx_new = *vx;
	double vy_new = *vy;
	verl_vel(&vx_new,&vy_new,fx1,fy1,fx2,fy2,m,dt);
	*vx = vx_new;
	*vy = vy_new;

}


void calculate_force(double* fx,double* fy,double rx1,double ry1,double rx2,double ry2,double m1,double m2,double G)
{
	assert(fx);
	assert(fy);
	double rx,ry;
	vec_diff(&rx,&ry,rx1,ry1,rx2,ry2);
	double norm  = sqrt(pow(rx,2)+pow(ry,2));
	double force_mag = G*m1*m2/(pow(norm,3));
	*fx = force_mag*rx;
	*fy = force_mag*ry;
}

void particle_force(double* fx,double* fy,double rx1,double ry1,double rx2,double ry2,double rxp,double ryp,double m1,double m2,double mp,double G)
{
	assert(fx);
	assert(fy);
	double rx1p,ry1p;
	double rx2p,ry2p;
	vec_diff(&rx1p,&ry1p,rxp,ryp,rx1,ry1);
	vec_diff(&rx2p,&ry2p,rx2,ry2,rxp,ryp);
	double norm1 = sqrt(pow(rx1p,2)+pow(ry1p,2));
	double norm2 = sqrt(pow(rx2p,2)+pow(ry2p,2));
	double force_mag1 = G*mp*m1/(pow(norm1,3));
	double force_mag2 = G*mp*m2/(pow(norm2,3));
	*fx = force_mag1*rx1p+force_mag2*rx2p;
	*fy = force_mag1*ry1p+force_mag2*ry2p;
}

/* Velocity verlet function that takes values at time t as input and returns values at t+dt */
void verl_pos(double* rx,double* ry,double vx,double vy,double fx,double fy,double m, int dt)
{
	assert(rx);
	assert(ry);
	*rx = (*rx) + vx*dt + (fx*pow(dt,2))/(2*m);
	*ry = (*ry) + vy*dt + (fy*pow(dt,2))/(2*m);
}

void verl_vel(double* vx,double* vy,double fx1,double fy1,double fx2,double fy2,double m, int dt)
{
	assert(vx);
	assert(vy);
	*vx = *vx + ((fx1+fx2)*dt)/(2*m);
	*vy = *vy + ((fy1+fy2)*dt)/(2*m);
}

void rot_update_pos(double* rx,double* ry,double vx,double vy,double ax,double ay,double dt)
{
	double rx_new = *rx;
	double ry_new = *ry;
	rotverl_pos(&rx_new,&ry_new,vx,vy,ax,ay,dt);
	* rx = rx_new;
	* ry = ry_new;
}

void rot_update_vel(double* vx,double* vy,double ax1,double ay1,double ax2,double ay2,double dt)
{
	double vx_new = *vx;
	double vy_new = *vy;
	rotverl_vel(&vx_new,&vy_new,ax1,ay1,ax2,ay2,dt);
	*vx = vx_new;
	*vy = vy_new;

}

void rot_acc(double* ax,double* ay,double vx,double vy,double x,double y,double rho1,double rho2,double a)
{
	assert(ax);
	assert(ay);
	*ax = x+vy - (1-a)*(x+a)/(pow(rho1,3)) - a*(x+a-1)/(pow(rho2,3));
	*ay = y-vx - (1-a)*(y)/(pow(rho1,3)) - a*(y)/(pow(rho2,3));
}

double pseudo_pot(double* omega, double x, double y,double a)
{
	assert(omega);
	double d = sqrt(pow((x+a),2)+pow(y,2));
	double r = sqrt(pow((x+a-1),2)+pow(y,2));
	*omega = 0.5*(pow(x,2)+pow(y,2))+ (1-a)/d + a/r;
}


void rotverl_pos(double* rx,double* ry,double vx,double vy,double ax,double ay, double dt)
{
	assert(rx);
	assert(ry);
	*rx = (*rx) + vx*dt + (ax*pow(dt,2))/(2);
	*ry = (*ry) + vy*dt + (ay*pow(dt,2))/(2);
}

void rotverl_vel(double* vx,double* vy,double ax1,double ay1,double ax2,double ay2, double dt)
{
	assert(vx);
	assert(vy);
	*vx = *vx + ((ax1+ax2)*dt)/(2);
	*vy = *vy + ((ay1+ay2)*dt)/(2);
}


void vec_diff(double* rx, double* ry, double rx1, double ry1,double rx2,double ry2)
{
	assert(rx);
	assert(ry);
	*rx = rx2 - rx1;
	*ry = ry2 - ry1;
}
