#pragma once

#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>

void print_peak_memory_usage()
{
    rusage usage;
    int ret = getrusage(RUSAGE_SELF, &usage); 

    if (ret == 0)
    {
        std::cout << "peak memory usage: " << usage.ru_maxrss << '\n';
    }
    else 
    {
        std::cout << "Couldn't determine peak memory usage.\n";
    }
}