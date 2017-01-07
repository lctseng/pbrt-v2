#pragma once

#include <chrono>

#define MEASURE_TIMING 1

#if MEASURE_TIMING
#define BEGIN_TIMING(FUNCTION_NAME) \
	static long long total_time = 0.0; \
	printf("\n" #FUNCTION_NAME " started\n"); \
	auto start = std::chrono::high_resolution_clock::now(); \

#define END_TIMING(FUNCTION_NAME) \
	auto delta = std::chrono::high_resolution_clock::now() - start; \
	auto time = std::chrono::duration_cast<std::chrono::milliseconds>(delta).count(); \
	total_time += time; \
	printf("\n" #FUNCTION_NAME " used: %lld ms (total %lld ms)\n", time, total_time);

#else
#define BEGIN_TIMING(FUNTION_NAME)
#define END_TIMING(FUNTION_NAME)
#endif