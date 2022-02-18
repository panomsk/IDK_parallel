/*
Parallel Dynamic Itemset Counting for Intel Many Integrarted Core architecture
(Original algorithm in
Brin S., Rajeev Motwani R., Ullman J.D., Tsur S.
Dynamic Itemset Counting and Implication Rules for Market Basket Data //
Proc. of the 1997 ACM SIGMOD Int. Conf. on Management of Data, May 13-15, 1997, Tucson, Arizona, USA. ACM, 1997. P. 255-264.).

Simple profiler.

(c) 2016 Mikhail Zymbler
*/

#ifndef PROFILER_H
#define PROFILER_H

#include <omp.h>

#define PROFILE

extern double start1, finish1, start2, finish2, start3, finish3, start4, finish4, start5, finish5, start6, finish6, start7, finish7, start8, finish8, start9, finish9, start10, finish10;

// Place START and FINISH tags
#ifdef NPROFILE
#define PRF_START(name)
#define PRF_FINISH(name)
#else
#define PRF_THREAD	(0)
#define PRF_START(startpoint) do { if (omp_get_thread_num() == PRF_THREAD) startpoint += omp_get_wtime(); } while (0);
#define PRF_FINISH(finishpoint) do { if (omp_get_thread_num() == PRF_THREAD) finishpoint += omp_get_wtime(); } while (0);
#endif

// Init profile tags
#ifdef NPROFILE
#define PRF_INIT
#else
#define PRF_INIT do { start1=finish1=start2=finish2=start3=finish3=start4=finish4=start5=finish5=start6=finish6=start7=finish7=start8=finish8=start9=finish9=start10=finish10=0.0; } while (0);
#endif

#endif
