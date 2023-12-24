/* See LICENSE file for copyright and license details.
  *
  * To understand liquid, start reading main().
  */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#define URAND ((double)rand() / (double)RAND_MAX)

typedef float** config_t;

typedef struct {
	int n;
	int nlateral;
	float fillfract;
	float density;
	float boxlenght;
	config_t config;
} system_t;

static void die(char *errstr, ...);
static config_t initialconfiguration(system_t *s);
static system_t* newsystem(int nlateral, float fillfract);
static void printsystem(system_t *s);
static void thermalize(system_t *s, int stepsbyparticle, float acceptratio);

void
die(char *errstr, ...)
{
	va_list ap;
	va_start(ap, errstr);
	vfprintf(stderr, errstr, ap);
	va_end(ap);
	exit(1);
}

config_t
initialconfiguration(system_t *s)
{
	config_t config = (config_t)calloc(s->n, sizeof(float *));
	if (!config) die("Cannot malloc!\n");
	for (int i = 0; i < s->n; i++)
	{
		config[i] = (float*)calloc(3, sizeof(float));
		if (!config[i]) die("Cannot malloc!\n");
	}

	float gapsize = cbrt(s->density);

	int p = 0;
	for (int i = 0; i < s->nlateral; i++)
		for (int j = 0; j < s->nlateral; j++)
			for (int k = 0; k < s->nlateral; k++)
			{
				config[p][0] = ((float)i + 0.5) * gapsize;
				config[p][1] = ((float)j + 0.5) * gapsize;
				config[p][2] = ((float)k + 0.5) * gapsize;
				p++;
			}

	return config;
}

void
torustopology(system_t *s)
{
	for (int i=0; i < s->n; i++)
	{
		s->config[i][0] = fmod(fabs(s->config[i][0]), s->boxlenght);
		s->config[i][1] = fmod(fabs(s->config[i][1]), s->boxlenght);
		s->config[i][2] = fmod(fabs(s->config[i][2]), s->boxlenght);
	}
}

system_t*
newsystem(int nlateral, float fillfract)
{
	system_t* s = (system_t*)malloc(sizeof(system_t));
	if (!s) die("error");

	int n = pow(nlateral, 3);
	float density = 6 * fillfract / M_PI;
	float boxlenght = cbrt((float)n / density);

	s->n = n;
	s->nlateral = nlateral;
	s->boxlenght = boxlenght;
	s->density = density;
	s->fillfract = fillfract;
	s->config = initialconfiguration(s);

	return s;
}

void
printsystem(system_t *s)
{
	printf("Number of particles: %i\n", s->n);
	printf("Box lenght: %f\n", s->boxlenght);
	printf("Configuration:\n");
	for (int i = 0; i < s->n; i++)
	{
		float x = s->config[i][0];
		float y = s->config[i][1];
		float z = s->config[i][2];
		printf("p_%i = (%.2f, %.2f, %.2f)\n", i, x, y, z);
	}
}

float
radialpotential(float r)
{
	float dl = 50;
	float dt = 1.4737;

	float da = dl * pow((dl / (dl - 1)), (dl - 1));
	float v = (da / dt) * (pow((1 / r), dl) - pow((1 / r), (dl - 1))) + 1 / dt;
	return v;
}

float
molecularpotential(system_t *s, int i, int j)
{
	float xi = fmod(fabs(s->config[i][0]), s->boxlenght);
	float yi = fmod(fabs(s->config[i][1]), s->boxlenght);
	float zi = fmod(fabs(s->config[i][2]), s->boxlenght);
	float xj = fmod(fabs(s->config[j][0]), s->boxlenght);
	float yj = fmod(fabs(s->config[j][1]), s->boxlenght);
	float zj = fmod(fabs(s->config[j][2]), s->boxlenght);
	float x2 = pow(fabs(xi - xj), 2);
	float y2 = pow(fabs(yi - yj), 2);
	float z2 = pow(fabs(zi - zj), 2);
	float r = sqrt(x2 + y2 + z2);

	if (r >= s->boxlenght / 2)
		return 0;

	return radialpotential(r);
}

float
particleenergy(system_t *s, int i)
{
	float acc = 0;
	for (int j = 0; j < s->n && j != i; j++)
		acc += molecularpotential(s, i, j);

	return acc;
}

float
systemenergy(system_t *s)
{
	float acc = 0;
	for (int i = 0; i < s->n; i++)
		for (int j = 0; j < i; j++)
			acc += molecularpotential(s, i, j);
	return acc;
}

float
trymoveaparticle(system_t *s, float *E, float *drmax, int *acc, int *att)
{
	int i = rand() % s->n;
	float energy = particleenergy(s, i);
	float ri[3] = {s->config[i][0], s->config[i][1], s->config[i][2]};
	float dx = *drmax * (URAND - 0.5);
	float dy = *drmax * (URAND - 0.5);
	float dz = *drmax * (URAND - 0.5);
	s->config[i][0] += dx;
	s->config[i][1] += dy;
	s->config[i][2] += dz;
	torustopology(s);
	
	float energynew = particleenergy(s, i);
	float denergy = energynew - energy;

	if (expf(-denergy) > URAND)
	{
		(*acc)++;
		*E += denergy; 
	} else {
		s->config[i][0] = ri[0];
		s->config[i][1] = ri[1];
		s->config[i][2] = ri[2];
	}
	(*att)++;
}

void
adjustdrmax(float *drmax, int acc, int att, float acceptratio)
{

	if ((float)acc / att < acceptratio)
		*drmax *= 0.95;
	else
		*drmax *= 1.05;
}

void
thermalize(system_t *s, int stepsbyparticle, float acceptratio)
{
	float drmax = 0.1;
	float E = systemenergy(s);
	int acc = 0;
	int att = 0;
	int cycles = stepsbyparticle * s->n / acceptratio;
	for (int i = 0;  i < cycles; i++)
	{
		trymoveaparticle(s, &E, &drmax, &acc, &att);
		adjustdrmax(&drmax, acc, att, acceptratio);
		printf("energy: %f\tdrmax: %f\tratio: \%f\n", E, drmax, (float)acc / att);
	}
}

void
savesystem(system_t *s)
{
	FILE *f = fopen("init.dat", "w");
	if (!f)
		die("Cannot open file.");
	for (int i = 0; i < s->n; i++)
	{
		float x = s->config[i][0];
		float y = s->config[i][1];
		float z = s->config[i][2];
		fprintf(f, "%f\t%f\t\%f\n", x, y, z);
	}
	fclose(f);
}

int
main(int argc, char *argv[])
{
	srand(time(0));
	system_t *s = newsystem(5, 0.35);
	thermalize(s, 5, 0.3);
	savesystem(s);

	return 0;
}
