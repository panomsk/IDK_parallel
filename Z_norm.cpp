#include "Z_norm.h"

void Z_norm(itemType* query, int lenQuery)
{
	itemType mean = 0;
	for (int i = 0; i < lenQuery; i++)
		mean = mean + query[i];
	mean = mean / lenQuery;

	itemType variance = 0;
	for (int i = 0; i < lenQuery; i++)
		variance = variance + (query[i] - mean)*(query[i] - mean);
	variance = variance / lenQuery;
	
	itemType sdiviation = sqrt(variance);
	for (int i = 0; i < lenQuery; i++)
		query[i] = (query[i] - mean) / sdiviation;

    //!! debug
	//printf("Mean = %f\n",mean);
	//printf("Variance = %f\n",variance);
	//printf("sdiviation = %f\n",sdiviation);
}


