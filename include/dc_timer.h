#ifndef DC_TIMER_H
#define DC_TIMER_H

/* 为实验报告提供毫秒级高分辨率计时；Windows 使用性能计数器。 */
#ifdef _WIN32
#include <windows.h>

static inline double dc_timer_now_milliseconds(void)
{
    LARGE_INTEGER frequency;
    LARGE_INTEGER counter;
    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&counter);
    return 1000.0 * (double)counter.QuadPart / (double)frequency.QuadPart;
}
#else
#include <time.h>

static inline double dc_timer_now_milliseconds(void)
{
    struct timespec value;
    timespec_get(&value, TIME_UTC);
    return 1000.0 * (double)value.tv_sec + 1.0e-6 * (double)value.tv_nsec;
}
#endif

#endif
