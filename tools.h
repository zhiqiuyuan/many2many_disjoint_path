#ifndef _TOOLS_H
#define _TOOLS_H

#include <iostream>
#include <fstream>
#include <iomanip>

#include <functional>

#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <deque>
#include <cstring>
#include <algorithm>

#include <chrono>

#include <emmintrin.h>
#include <immintrin.h>

#include "assert.h"
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

#include "config.h"

#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

// return [0,possible_max)
long long Rand(long long possible_max);

#endif