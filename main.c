/* See LICENSE file for copyright and license details.
  *
  * To understand liquid, start reading main().
  */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

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

int
main(int argc, char *argv[])
{
	system_t *s = newsystem(3, 0.35);
	printsystem(s);

	return 0;
}
