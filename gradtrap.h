#include <stdio.h>
#include <math.h>

int gradtrap(float dk, float sm, float gm, float dt, float* h_ref, float* t_ramp_ref, float* t_plat_ref) {

	float gam = 4.258*2*M_PI; /* gyromagnetic ratio (rad/G/ms) */
	float A;
	float h, t_ramp, t_plat;

	/* calculate required area under gradient waveform */
	A = dk / gam;

	/* calculate height of waveform */
	h = sqrt(A * sm);
	if (h > gm) h = gm; /* limit height to max gradient */

	/* calculate the ramp and plateau time and round up to nearest sampling interval */
	t_ramp = h / sm;
	t_plat = A / h - t_ramp;

	/* round up to nearest sampling interview & correct h */
	t_ramp = dt * ceil(t_ramp / dt);
	t_plat = dt * ceil(t_plat / dt);
	h = A / (t_ramp + t_plat);
	
	/* assign reference values */
	*h_ref = h;
	*t_ramp_ref = t_ramp;
	*t_plat_ref = t_plat;

	return 1;
}
