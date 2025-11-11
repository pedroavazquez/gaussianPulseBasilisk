/*
 * =====================================================================================
	Simulation of jet excited with gaussian pulse

	Build a channel using multigrid with dimensions (nx=NX, ny=1, nz=1)

	Save snapshots in dump files. These file can be read with readLateralJet.c

	Usage:
		CC='mpicc -D_MPI=_ -DNX=_ -DRESTORE=_ -DMU1=_ -DMU2=_  -DTEND=_ -DDTOUT=_
		  -DDTDUMP=_  -DMOVIE=_ -DDTMOVIE=_ -DDTMAXMINE=_ -DLEVEL=_ 
		  -DWEBER=_ -DK0=_ -DA0=_ -DSIGMAZ=_ -DPHI=_'  make gaussianPulse.tst

		Parameters
			MPI 			= number of processes
			NX			= number of boxes (must be = MPI) (1)
			RESTORE			= 0 -> run from scratch (Default)
                       				  1 -> restart simulation fom file "dump"
			MU1			= Dynamic viscosity of inner fluid (0.001 Pa.s)
			MU2			= Dynamic viscosity of outer fluid (0.1 Pa.s)
			TEND			= Final simulation time
			DTOUT			= Snapshots saved every DTOUT
			DMOVIE			= if 1 -> snapshots every DMOVIE for mp4 movie
			DTDUMP			= dumps saved every DTDUMP
			DTMAXMINE		= maximum value of the time step
			LEVEL			= Refinement level
			WEBER			= Weber number
			K0			= Wavelength number
			A0			= Amplitud of the perturbation
			PHI 			= Phase


 * =====================================================================================
 */

#include "grid/multigrid.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/double-projection.h"
//#include "navier-stokes/perfs.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"
#include "tag.h"
#include "drop_stat.h"
#include "outputFacetsMine.h"

#define PI (3.14159265359)

#define R1  (1.0)
#define R2 (4.0)			// radius of the channel and size of the box
#define LBOX (R2)		  // size of one box
#ifndef NX						// number of boxes (default 1)
#define NX 1
#endif
#define LX (LBOX*NX)  // length of the channel

#ifndef RESTORE				// 0 -> run from scratch
#define RESTORE 0			// 1 -> restart simulation fom file "dump"
#endif

#ifndef WEBER
#define WEBER	8.
#endif

#ifndef K0
#define K0	0.69
#endif

#ifndef PHI
#define PHI 	0.
#endif

#ifndef A0
#define A0	0.01
#endif

#ifndef SIGMAZ
#define SIGMAZ	2.0
#endif


#ifndef MU1
#define MU1 1.000		// Dynamic viscosity of inner fluid
#endif
#ifndef MU2
#define MU2 1.e-6	// Dynamic viscosity of outer fluid
#endif

#ifndef TEND				// Final simulation time
#define TEND 165.
#endif
#ifndef DTOUT				// Snapshots saved every DTOUT
#define DTOUT 2.0
#endif

#ifndef MOVIE			// Snapshots for movie
#define MOVIE  0
#endif
#ifndef DTMOVIE			// Time intervalas for the movie
#define DTMOVIE  (TEND)
#endif
#ifndef DTDUMP		// dump file saved every DTDUMP
#define DTDUMP 10.
#endif
#ifndef DTMAXMINE
#define DTMAXMINE               0.01
#endif

#ifndef LEVEL				// Refinement level
#define LEVEL 8
#endif

#define OH1	1.e-2
#define OH2	1.e-5

//#define mu(f) ((clamp(f,0,1)*MU1 + (1.-clamp(f,0,1))*MU2))


// Profile for inlet velocity of the inner fluid
double Uin1 (double t) {
	double  We =  WEBER;
	//double  k  =  0.69;
	double  k  =  K0;
	//double  A  =  0.01;
	double  A  =  A0;
	//double  sigma_z =  2.;
	double  sigma_z =  SIGMAZ;
	double  sigma = sigma_z/sqrt(We);
	double  t_mid  = 5.*sigma;
	double  t_0  = 5.*sigma;
	double  phi  = PHI;
	double dt2 = (t-t_mid)*(t-t_mid);
	double sigma2 = sigma*sigma;
	double ft = dt2/sigma2;
        return  sqrt(We) + A*exp(-0.5*ft)*sin(k*sqrt(We)*(t-t_0)+phi);
}


/*
=================================================
Boundary conditions
=================================================
*/

u.n[left] = dirichlet(Uin1(t));
u.t[left] = dirichlet(0.0);
u.n[top] = dirichlet(0.0);
u.t[top] = neumann(0.0);
u.n[right] = neumann(0.0);
u.t[right] = neumann(0.0);
u.n[bottom] = neumann(0.0);
u.t[bottom] = neumann(0.0);

p[top] = dirichlet(0.);
pf[top] = dirichlet(0.);
p[right] = neumann(0.);
pf[right] = neumann(0.);

f[left] = (y<R1 ? 1.0 : 0.);
f[top] = 0.;

double xFront = 0.;
double delMax = 0.;

scalar del[];


int main (int argc, char * argv[])
{
  dimensions ( nx = NX, ny = 1, nz = 1);
  size ( LX ) ;
  init_grid ( 1 << LEVEL );
  origin ( 0., 0.);

  rho1 = 1., rho2 = 1.e-3; // densities of inner and outer fluids
  mu1 = OH1, mu2 = OH2;    // dynamic viscosities of inner and outer fluids
  f.sigma = 1.;            // surface tension coefficient
  stokes = false;          // stokes flow (true) or Navier--Stokes (false)

  TOLERANCE = 1.e-6 [*];
  DT = DTMAXMINE;

  run();

}

event init ( t = 0 )
{
	// Find Delta
	foreach()
		del[] = Delta;
	delMax = statsf(del).max;

	// Save information to stdout
	fprintf (stdout, "# NX    = %d\n", NX);
	fprintf (stdout, "# R2    = %g\n", R2);
	fprintf (stdout, "# WEBER = %g\n", WEBER);
	fprintf (stdout, "# K0    = %g\n", K0);
	fprintf (stdout, "# A0    = %g\n", A0);
	fprintf (stdout, "# PHI   = %g\n", PHI);
	fprintf (stdout, "# SIGMAZ= %g\n", SIGMAZ);
	fprintf (stdout, "# NX    = %d\n", NX);
	fprintf (stdout, "# LEVEL   = %d\n", LEVEL);
	fprintf (stdout, "# delMax  = %g  %g\n", delMax, LX/(1 << LEVEL));
	fprintf (stdout, "# DTMAX   = %g\n", DT);
	fprintf (stdout, "\n");

#if !RESTORE  // Start from scratch

	// Set initial values
	foreach() {
		u.x[] = sqrt(WEBER);
		u.y[] = 0.;
		f[] = ( (y<R1) ? 1.0 : 0);
		//f[] = 0.;
	}
#else  // Start from previous run

	// Read simulation from file dump
	restore(file = "dump");

	// Restore values of metric factors (essential in axi)
  foreach()
    cm[] = y;
  cm[top] = dirichlet(y);
  cm[bottom] = dirichlet(y);
  foreach_face()
    fm.x[] = max(y, 1./HUGE);
  fm.t[top] = dirichlet(y);
  fm.t[bottom] = dirichlet(y);
#endif
}

event defaults (i = 0) {

/**
If the acceleration vector *a* (defined by the Navier--Stokes solver) is constant, we make it variable. 
*/
  if (is_constant (a.x))
    a = new face vector;

  /**
  The restriction/refine attributes of the charge density are those of a tracer
  otherwise the conservation is not guaranteed. */

  /**
  By default the permittivity is unity and other quantities are zero. */

}

// Compute the length of the jet
event front ( i += 1 ) {
	scalar m[];
	double THR = 1.e-2;
	foreach()
		m[] = (f[]>THR);
	tag(m);
	
	scalar c[];
	foreach()
		c[] = m[]==1 ? f[] : 0;
	scalar pos[];
	position (c, pos, {1., 0.});
	xFront = statsf(pos).max;
	if (xFront < 0 )
		xFront = 0.;
}


event logfile (i++) {
  if (i == 0)
    fprintf (stdout,
	     "# 1:t  2:dt  3:i  4:xFront  5:mgp.i  6:mgpf.i  7:mgu.i\n");
  fprintf (stdout, "%10.8e %6.4e %d %g %d %d %d\n", 
	   t, dt, i, xFront, mgp.i, mgpf.i, mgu.i);
}

// remove bubbles
#if 1
event remove_bubbles_droplets ( i+=1 ) {

#if 0
	// remove bubles
  remove_droplets (f, minsize = 0, bubbles = true);
  //remove_droplets (f);

  scalar m[];
  double THR = 1e-2; //THRESHOLD
  foreach()
    m[] = f[] > THR;
  tag (m);
  foreach()
    f[] = m[] > 1 ? 0. : f[];
#endif

#if 0
  // remove droplets
  remove_droplets_vol (f, 1, true);
#endif

  // Remove droplets of volume less than 0.01
  xFront = remove_droplets_volmin (f, 0.01);
#if 0
  double hh = LX/(1 << LEVEL);
 if (t>=100.) 
   xFront = remove_droplets_volmin (f, 0.01);
 else
   xFront = remove_droplets_volmin (f, 0.0);
#endif
#if 0
  // remove bubbles
  remove_droplets_vol (f, 0, false);
#endif
}
#endif

event dumps( t += DTDUMP; t <= TEND )
{
	// Save the simulation
	char nameDump[200];
	sprintf(nameDump, "dump-%014.4f", t);
	dump ( file = nameDump);
	dump ( file = "dump");

}

#if MOVIE
event movie ( t += DTMOVIE; t <= TEND )
{

	// Save movie as mp4 files
	char label[200];
	view( tx = -1.5, ty = -0.3, sx = 3.0, sy = 20.0, width = 1200, height = 300 );

	// u.x
	//cells();
	box();
	sprintf(label, "Ux t:%2.1f", t);
	squares("u.x", spread = -1, linear=true, cbar = true, pos = {0.6, -0.2}, levels = 10, mid = true, format = "%6.4f", label = label);
	//squares("u.x", min = -2., max = 6., cbar = true, pos = {0.5, -0.2}, levels = 10, mid = true, format = "%6.4f", label = label);
	draw_vof("f", lw = 1, lc = {0,0,0} );
	//vectors("u", scale = 0.05);
	save("ux.mp4");

	// f
	//cells();
	clear();
	box();
	squares("f", min = 0., max = 1);
	draw_vof("f", lw = 1, lc = {1,1,1} );
	save("f.mp4");

}
#endif

event snapshots ( t += DTOUT; t <= TEND )
{

	// Save plots as png files
	char label[200];
	view( tx = -1.5, ty = -0.3, sx = 3.0, sy = 10.0, width = 1200, height = 300 );

	// u.x
	//cells();
	box();
	sprintf(label, "Ux t:%2.1f", t);
	squares("u.x", spread = -1,  linear=true, cbar = true, pos = {0.6, -0.2}, levels = 10, mid = true, format = "%6.4f", label = label);
	//squares("u.x", min = 0.8*U1, max = 2*U1, cbar = true, pos = {0.5, -0.2}, levels = 10, mid = true, format = "%6.4f", label = label);
	draw_vof("f", lw = 1, lc = {0,0,0} );
	//vectors("u", scale = 0.05);
	save("ux.png");

	// f
	clear();
	cells();
	box();
	squares("f", min = 0., max = 1., cbar = false);
	draw_vof("f", lw = 1, lc = {1,1,1} );
	save("f.png");

}

double radius ( double xx ) {

	double YY0 = 0.0;
	double YY1 = R2;
	double ddy = delMax/4.0;
	double yy = YY0;
	while ( yy < YY1) {
		if (interpolate(f, xx, yy) < 0.5) 
			return yy;
		yy += ddy;
	}
	return yy;

}

event saveProfile ( t += 0.1; t <= TEND ) {
	// Save profile of the velocity at the inlet
	char name[200];
	sprintf(name, "profile-%014.4f-%02d", t, pid());

	FILE * fp = fopen(name, "w");
	output_facetsMine(f, fp);
	fclose(fp);
#if 0
	int pid = pid();
	double XX0 = LBOX*pid;
	double XX1 = LBOX * (pid + 1);
	double ddx = 0.1;
	double ddy = LX/(1 << LEVEL);
	double yy = 0.0;
		FILE * fp = fopen(name, "w");
		double xx = XX0;
		while ( xx < XX1 ) {
			yy = 0.0;
			while ( yy < R2 ) {
				double ff = interpolate(f, xx, yy);
				if ( ff < 0.5) break;
				yy += ddy;
			}
			fprintf(fp, "%6.4f  %6.4f\n", xx, yy);
			xx += ddx;
		}
		fclose(fp);
	#endif
}
