#include "ReadData.h"

void read_timeseries(const char* filename, itemType *S, itemType **R)
{	
	FILE* fileS;
	// errno_t err;
	fileS = fopen(filename, "r"); 
	
	assert(fileS != NULL); 
	
	/*err = fopen_s(&fileS, filename, "r");
	assert(err == 0);*/
	char line[lengthsrt];

	int i = 0;
	while (fgets(line, lengthsrt, fileS))
	{
		char* tmp = strdup(line); // íà gcc
		
		S[i] = atof(getelem(tmp, 1));
		free(tmp);
		for (int j = 0; j < d_par; j++)
		{
			char* tmp = strdup(line);
			//char* tmp = _strdup(line);
			R[j][i] = atof(getelem(tmp, j + 2));
			free(tmp);
		}
		i++;
		if (i >= L_par + numimputevalue)
			break;
		
	}

}


const char* getelem(char* line, int num)
{
	// char* next_token1 = NULL;
	const char* tok;
	//!!for (tok = strtok(line, ";"); tok && *tok; tok = strtok(NULL, ";\n"))
	for (tok = strtok(line, ";"); tok && *tok; tok = strtok(NULL, ";\n"))
	// for (tok = strtok_s(line, ";", &next_token1); tok && *tok; tok = strtok_s(NULL, ";\n", &next_token1))
	{
		if (!--num)
		{
			return tok;
		}
	}
	return NULL;
}


