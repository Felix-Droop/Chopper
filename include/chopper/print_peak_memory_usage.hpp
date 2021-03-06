#pragma once

#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>

void print_peak_memory_usage()
{
    rusage usage;

    if (getrusage(RUSAGE_SELF, &usage) == 0)
    {
        std::cerr << "peak memory usage: " << usage.ru_maxrss << " kilobytes" << std::endl;
    }
    else 
    {
        std::cerr << "couldn't determine peak memory usage" << std::endl;
    }
}