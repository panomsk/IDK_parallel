#include "DTW.h"

/*
//исходный вариант DTW не учитывает r_par
itemType DTW(itemType* Q, itemType* C)
{
	double DTWarray[l_par+1][l_par+1];
	
	for (int i = 0; i < l_par+1; i++)
	{
		for (int j = 0; j < l_par+1; j++)
			DTWarray[i][j] = inf_val;
	}
	DTWarray[0][0] = 0;
	for (int i = 1; i < l_par+1; i++)
	{
		for (int j = 1; j < l_par+1; j++)
		{
			itemType dt = dist(Q[i-1], C[j-1]);
			itemType last_min = min(DTWarray[i - 1][j - 1], DTWarray[i - 1][j]);
			last_min = min(last_min, DTWarray[i][j - 1]);
			DTWarray[i][j] = last_min + dt;
		}
	}

	return DTWarray[l_par][l_par];
}*/

/*itemType DTW_ver1(itemType* Q, itemType* C, int r)
{
	double DTWarray[l_par+1][l_par+1];
	
	for (int i = 0; i < l_par+1; i++)
	{
		for (int j = 0; j < l_par+1; j++)
			DTWarray[i][j] = inf_val;
	}

	DTWarray[0][0] = 0;

	for (int i = 1; i <= l_par; i++)
	{
		for (int j = max(1,i-r); j <= min(l_par,i+r); j++)
		{
			DTWarray[i][j]=0;
		}
	}

	for (int i = 1; i <= l_par; i++)
	{
		for (int j = max(1,i-r); j <= min(l_par,i+r); j++)
		{
			itemType dt = dist(Q[i-1], C[j-1]);
			itemType last_min = min(DTWarray[i - 1][j - 1], DTWarray[i - 1][j]);
			last_min = min(last_min, DTWarray[i][j - 1]);
			DTWarray[i][j] = last_min + dt;
		}
	}

	return DTWarray[l_par][l_par];
}*/

//!! где-то ошибка несимметрично отрабатывает функция
/*itemType DTW_ver2(itemType* Q, itemType* C, int r)
{
    itemType cost[l_par];
    itemType cost_prev[l_par];
    itemType temp;
    
    int i, j;

    for (i = 0; i <= (l_par-1); i++)
    {
        cost[i] = inf_val;
        cost_prev[i] = inf_val;
    }

    cost_prev[0] = dist(Q[0],C[0]);

    for (j = 1; j <= min(l_par-1,r); j++)
    {
        cost_prev[j] = cost_prev[j-1] + dist(Q[0],C[j]);
    }

    for (i = 1; i<=(l_par-1); i++)
    {
        for (j = max(0,i-r); j <= min(l_par-1,i+r); j++)
        {
            if ((j-1)<0)
                cost[j] = cost_prev[j];
            else
                cost[j] = min(cost_prev[j],cost_prev[j-1]);
        }

        for (j = max(0,i-r); j <= min(l_par-1,i+r); j++)
        {
            if ((j-1) < 0)
                cost[j] = dist(Q[i],C[j]) + cost[j];
            else
                cost[j] = dist(Q[i],C[j]) + min(cost[j],cost[j-1]);            
        }
        
        for (j = 0; j <= (l_par-1); j++)
        {
            temp = cost[j];
            cost[j] = cost_prev[j];
            cost_prev[j] = temp;
        }
    }
    return cost_prev[l_par-1];
}   
*/

//!! fast dtw возможно нужна ссылка - уточнить
itemType DTW(itemType* q, itemType* c, int r)		
{
	itemType dtw = 0;
	itemType *cost;
    itemType *cost_prev;
    itemType *cost_tmp;
    int i, j, k;
    itemType x, y, z, min_cost;

    // Instead of using matrix of size O(m^2) or O(mr), we will reuse two array of size O(r).
    cost = (itemType*)malloc(sizeof(itemType)*(2*r+1));
    for(k = 0; k < 2 * r + 1; k++)    
		cost[k] = inf_val;

    cost_prev = (itemType*)malloc(sizeof(itemType)*(2*r+1));
    for(k = 0; k < 2 * r + 1; k++)    
		cost_prev[k] = inf_val;

    for (i = 0; i < l_par; i++)
    {
        k = max(0, r-i); 
        min_cost = inf_val;

        for(j = max(0,i - r); j <= min(l_par - 1, i + r); j++, k++)
        {
            // Initialize all row and column
            if ((i == 0) && (j == 0))
            {
                cost[k] = dist(q[0], c[0]);
                min_cost = cost[k];
                continue;
            }

            if ((j - 1 < 0) || (k - 1 < 0))     
				y = inf_val;
            else                      
				y = cost[k - 1];
            if ((i - 1 < 0)||(k + 1 > 2 * r))   
				x = inf_val;
            else                      
				x = cost_prev[k + 1];
            if ((i - 1 < 0) || (j - 1 < 0))     
				z = inf_val;
            else                      
				z = cost_prev[k];

            // Classic DTW calculation
            cost[k] = min(min(x, y), z) + dist(q[i], c[j]);

            // Find minimum cost in row for early abandoning (possibly to use column instead of row).
            if (cost[k] < min_cost)
            {   
				min_cost = cost[k];
            }
        }

        cost_tmp = cost;
        cost = cost_prev;
        cost_prev = cost_tmp;
    }

    k--;

    // the DTW distance is in the last cell in the matrix of size O(m^2) or at the middle of our array.
    itemType final_dtw = cost_prev[k];
    free(cost);
    free(cost_prev);

	return final_dtw;
}