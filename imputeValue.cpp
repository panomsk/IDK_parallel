#include "imputeValue.h"

itemType Q[d_par][l_par];
itemType SetC[d_par][(L_par - (2 * l_par)) + numimputevalue][l_par];
itemType LB[d_par][(L_par - (2 * l_par)) + numimputevalue][num_lb];
itemType bsf[d_par];
itemType N[d_par][(L_par - (2 * l_par)) + numimputevalue];
int CandIndex[(L_par - (2 * l_par)) + numimputevalue];
int ind[k_par];
itemType RANK[(L_par - (2 * l_par)) + numimputevalue];
bool BitMap[d_par][(L_par - (2 * l_par)) + numimputevalue];
extern int num_of_threads;
extern int DTW_count;
extern int LB_count;

//!!проверка работы 
//itemType RANK_test[(L_par - (2 * l_par)) + numimputevalue];

struct maxRank{
	itemType Value;
	int Indax;
};
#pragma omp declare reduction(max : struct maxRank : \
	omp_out.Value = omp_in.Value < omp_out.Value ? omp_out.Value : omp_in.Value, \
	omp_out.Indax = omp_in.Value < omp_out.Value ? omp_out.Indax : omp_in.Indax ) \
	initializer( omp_priv = { 0, 0 } )


itemType impute(itemType* S, itemType** R, int h)
{

	itemType s = 0;
	PRF_START(start2);
	FillData(S, R, h);
	PRF_FINISH(finish2);
	PRF_START(start3);
	CslcZnorm(h);
	PRF_FINISH(finish3);

	//первое приближение bsf по произвольной серии опорного ряда после нормализации
	for (int i = 0; i < d_par; i++)
		bsf[i] = DTW(Q[i], SetC[i][(L_par - (3 * l_par))], r_par);
	
	//!! не используется
	//if (CalcLB(h) == -1) s = S[L_par + h - 1];
	CalcLB(h);
	
	PRF_START(start7);
	verticalOverlap(h, use_corr);
	int j = 0;
	int rank_counts = 0;
	while (j < k_par)
	{
		struct maxRank maxRANK = {0, 0};
		int i;
		#pragma omp parallel for num_threads(num_of_threads) reduction (max: maxRANK)
		for (i = 0; i < (L_par - (2 * l_par)) + h; i++)
		{
			if (maxRANK.Value < RANK[i])
			{ 
				maxRANK.Value = RANK[i];
				maxRANK.Indax = i;
			}
		}

		//!!проверка работы редукции
		/*struct maxRank test_rank = {0, 0};
		for (i = 0; i < (L_par - (2 * l_par)) + h; i++)
		{
			if (test_rank.Value < RANK[i])
			{ 
				test_rank.Value = RANK[i];
				test_rank.Indax = i;
			}
		}
		if (test_rank.Indax!=maxRANK.Indax)
		{
			printf("RANKS_ERROR reduction index = %d reduction value = %f for index = %d for value = %f\n", maxRANK.Indax, maxRANK.Value, test_rank.Indax, test_rank.Value);
		}*/




		//добавлена проверка отличия максимального ранка от 0. Если 0, значит больше нет обсчитанных DTW
		if(maxRANK.Value>0)
		{
			ind[j] = maxRANK.Indax;
			int t;
			// #pragma omp parallel for num_threads(num_of_threads)
			for (t = max(0, (ind[j] - l_par)); t < min ((ind[j] + l_par),((L_par - (2 * l_par)) + h)); t++)
			{
				RANK[t] = 0;
			}
			rank_counts+=1;
		}
		else break;
		j++;
	}

	//если нашли меньше рассчетных точек чем k_par
	assert(rank_counts==k_par);
	//ниже представлено более слабое условие, поставить при постоянной сработке условия выше
	//assert(rank_counts>0);

	#ifdef DEBUG
		//подсчет количества расчетов DTW
		DTW_count += d_par;
		LB_count += d_par*(L_par - (2 * l_par)) + h;	
		for (j = 0; j < d_par; j++)
		{
			for (int i = 0; i < (L_par - (2 * l_par)) + h; i++)
			{
				if (N[j][i] < inf_val)
				{
					DTW_count += 1;
				}
			}
		}
	#endif

	int n = 0;
	for (j = 0; j < k_par; j++) {
		if (ind[j] != -1) 
		{
			//приведено в соответствие с псевдокодом
			s = s + S[ind[j] + l_par -1];
			n++;
		}
	}
	//если нашли меньше рассчетных точек чем k_par
	assert((n==k_par)&&(k_par>=0));
	//ниже представлено более слабое условие, поставить при постоянной сработке условия выше
	//assert(n>0);
	
	PRF_FINISH(finish7);
	return s/n;
}

void FillData(itemType* S, itemType** R, int h)
{
	int i, j, q;
	for (j = 0; j < d_par; j++)
	{
		for (i = 0; i < l_par; i++)
			//исправлено - последняя точка запроса - восстанавливаемая точка в основном ряду
			Q[j][i] = R[j][L_par - l_par + i + h + 1];
	}
	for (i = 0; i < d_par; i++)
	{
		for (j = 0; j < (L_par - (2 * l_par)) + h; j++)
		{
			for (q = 0; q < l_par; q++)
				SetC[i][j][q] = R[i][j + q];
		}
	}
	for (i = 0; i < d_par; i++)
	{
		for (j = 0; j < (L_par - (2 * l_par)) + h; j++)
			N[i][j] = inf_val;
	}

	//не используется
	/*for (i = 0; i < d_par; i++)
	{
		for (j = 0; j < (L_par - (2 * l_par)) + h; j++)
			NInd[i][j] = -1;
	}*/

	for (i = 0; i < k_par; i++)
		ind[i] = -1;
	
	for (j = 0; j < (L_par - (2 * l_par)) + h; j++)
	{
		RANK[j] = 0;
	}

	for (i = 0; i < d_par; i++)
	{
		for (j = 0; j < (L_par - (2 * l_par)) + h; j++)
			BitMap[i][j] = true;
	}
	
}

void CslcZnorm(int h)
{
	int i ,j;
	#pragma omp parallel for num_threads(num_of_threads) collapse(2) 
	for (i = 0; i < d_par; i++)
	{
		for (j = 0; j < (L_par - (2 * l_par)) + h; j++) 
			Z_norm(SetC[i][j], l_par);
	}
	#pragma omp parallel for num_threads(num_of_threads)
	for (i = 0; i < d_par; i++)
	{
		Z_norm(Q[i], l_par);
	}
}

int CalcLB(int h)
{
	int i,j,s,cnt,left,right;
	itemType cur_dist;
	PRF_START(start4);
	#pragma omp parallel for num_threads(num_of_threads) collapse(2) private(j, i)
	for (i = 0; i < d_par; i++)
	{
		for (j = 0; j < ((L_par - (2 * l_par)) + h); j++)
		{
			LB[i][j][0] = LbKim(Q[i], SetC[i][j]);
			LB[i][j][1] = LbKeogh(Q[i], SetC[i][j], r_par);
			LB[i][j][2] = LbKeogh(SetC[i][j], Q[i], r_par);
		}
	}
	PRF_FINISH(finish4);
		
	for (i = 0; i < d_par; i++)
	{
		while(true)
		{
			PRF_START(start5);
			#pragma omp parallel for num_threads(num_of_threads)
			for (j = 0; j < ((L_par - (2 * l_par)) + h); j++)
			{
				BitMap[i][j] = BitMap[i][j] && (bsf[i] > LB[i][j][0]) && (bsf[i] > LB[i][j][1]) && (bsf[i] > LB[i][j][2]);
			}
			cnt = 0;
			for (j = 0; j < ((L_par - (2 * l_par)) + h); j++)
			{
				if (BitMap[i][j])
				{
					CandIndex[cnt] = j;
					cnt += 1;
				}
			}

			//!!printf("h = %d i = %d, cnt = %d bsf[i] = %f\n", h, i, cnt, bsf[i]);

			PRF_FINISH(finish5);
			if (cnt == 0) break;
			

			PRF_START(start6);
			s = 1;
			while(true)
			{
				left = num_of_threads*(s-1);
				right = min(cnt, num_of_threads*s);
				cur_dist = bsf[i];
				
				#pragma omp parallel for num_threads(num_of_threads)
				for (j = left; j < right; j++)
				{
					if (BitMap[i][CandIndex[j]])
					{
						N[i][CandIndex[j]] = DTW(Q[i], SetC[i][CandIndex[j]], r_par);
						BitMap[i][CandIndex[j]] = false;

						#pragma omp critical 
						{
							if (cur_dist > N[i][CandIndex[j]]) 
							{
								cur_dist = N[i][CandIndex[j]];
								//!!printf("cur_dist = %f  bsf[i] = %f  position = %d\n", cur_dist, bsf[i], CandIndex[j]);
							}
						}
					}
				}
				s+=1;
				if ((right==cnt)||(bsf[i]>cur_dist)) 
				{
					bsf[i] = min(bsf[i], cur_dist);
					break;
				}
			}
			PRF_FINISH(finish6);
		}
	}
	return 0;
}

void verticalOverlap(int h, bool corr)
{
	int i ,j;
		
	//!!исправлена ошибка параллельной работы вложенного цикла #pragma omp parallel for num_threads(num_of_threads) collapse(2) private(j, i)
	for (j = 0; j < d_par; j++)
	{
		#pragma omp parallel for num_threads(num_of_threads)
		for (i = 0; i < (L_par - (2 * l_par)) + h; i++)
		{
			if (N[j][i] < inf_val)
			{
				//RANK[i] = RANK[i] + ((1 / N[j][i]) * (((d_par - j) + 1) / d_par));
				if(corr == true)
				{
					//вес зависит от номера опорного ряда (ряды отсортированы в порядке убывания корреляции с основным рядом)
					RANK[i] = RANK[i] + ((1 / (N[j][i]+itemtype_epsilon)) * ((d_par - j) / d_par));
				}
				else
				{
					//веса всех опорных рядов одинаковые
					RANK[i] = RANK[i] + (1 / (N[j][i]+itemtype_epsilon));
				}
			}
		}
	}

	//проверка корректности параллельного вычисления Rank
	/*
	for (j = 0; j < (L_par - (2 * l_par)) + h; j++)
	{
		RANK_test[j] = 0;
	}
	for (j = 0; j < d_par; j++)
	{
		for (i = 0; i < (L_par - (2 * l_par)) + h; i++)
		{
			if (N[j][i] < inf_val)
			{
				RANK_test[i] = RANK_test[i] + ((1 / N[j][i]) * (((d_par - j) + 1) / d_par));
			}
		}
	}

	for (i = 0; i < (L_par - (2 * l_par)) + h; i++)
	{
		if(RANK_test[i]!=RANK[i])
		{
			printf("RANKS_ERROR_2 h = %d index = %d parallel value = %f for value = %f\n", h, i, RANK[i], RANK_test[i]);			
		}
	}*/
	

}