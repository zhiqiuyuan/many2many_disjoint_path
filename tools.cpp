#include "tools.h"

long long Rand(long long possible_max)
{
    double min_interval_width = possible_max / (double)RAND_MAX;
    long long re = rand();
    while (min_interval_width >= RAND_MAX)
    {
        re = re * RAND_MAX + rand();
        min_interval_width /= RAND_MAX;
    }
    return re % possible_max;
}