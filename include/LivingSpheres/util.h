#pragma once
#include <cmath>


double GenerateRandomNumber(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

int GenerateIntegerRandomNumber(int fMin, int fMax)
{

	return (int)round(GenerateRandomNumber((double)fMin, (double)fMax));
}


double DotProduct(double* a, double* b, double size)
{
	double dot_product = 0.0;

	for (int i = 0; i < size; i++)
	{
		dot_product = a[i] * b[i];
	}

	return dot_product;
}

int sgn(double x)
{
	return (x > 0) - (x < 0);
}