#pragma once

#ifdef ENABLE_APEX_PROFILING
#include "apex_api.hpp"
#endif

class Profiler {
public:
#ifdef ENABLE_APEX_PROFILING
    Profiler(const char* func) {
       p = apex::start(func);
    }

    ~Profiler(){
        apex::stop(p);
    }

    apex::profiler* p;
#else
    Profiler(const char* func){

    }
#endif
};
