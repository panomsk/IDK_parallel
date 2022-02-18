#ifndef DTW_H
#define DTW_H

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "config.h"

/*itemType DTW(itemType* Q, itemType* C);
itemType DTW_ver1(itemType* Q, itemType* C, int r);
itemType DTW_ver2(itemType* Q, itemType* C, int r);*/
itemType DTW(itemType* q, itemType* c, int r);

#endif // DTW_H