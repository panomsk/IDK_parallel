#ifndef CONFIG_H
#define CONFIG_H
#include <math.h>

#define L_par (50000) //стартовая точка основного ряда с которой начинаем восстановление
#define lengthsrt (5000)
#define numimputevalue (1000) // количество точек для восстановления

//Инициализация bsf: 
//0 - отступ в 3*l_par, 
//1 - в качестве начального приближения берем подпоследовательность, соответствующую точке основного ряда, 
//    предшествующей восстанавливаемой
//2 - идея Яны (ищем min LB[][3])
#define init_bsf (1)
#define cut_off_neighbors (true) 

#define d_par (3) // количество опорных рядов, используемых для восстановления
#define l_par (100) // ширина поискового окна
#define k_par (3) // количество точек для восстановления значения основного ряда

#define r_par (int(floor(l_par/4))) // ширина полосы Сако-Чиба (используется при расчете DTW и LB)
#define inf_val (1e+10)
#define itemtype_epsilon (2.2204460492503131e-16) // машинный epsilon для типа double
#define num_lb (4)
#define use_corr (false)//(false)

//#define dist(x,y) fabs((x)-(y))
#define dist(x,y) (((x)-(y))*((x)-(y)))
#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)>(y)?(x):(y))  

#define subseq_count ((L_par) - (l_par))

typedef double itemType;

#endif // CONFIG_H