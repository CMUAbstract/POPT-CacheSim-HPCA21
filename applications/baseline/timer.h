// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef TIMER_H_
#define TIMER_H_

#include <sys/time.h>
#include <ctime>
#include <chrono>


/*
GAP Benchmark Suite
Class:  Timer
Author: Scott Beamer

Simple timer that wraps gettimeofday
*/

#define COMPILER_BARRIER() asm volatile("" ::: "memory")
 
class Timer {
 public:
  Timer() { elapsed_time = 0.0;}

  void Start() {
    //gettimeofday(&start_time_, NULL);
    //clock_gettime(CLOCK_MONOTONIC, &start_time_);
    COMPILER_BARRIER(); 
    start_time_ = std::chrono::steady_clock::now();
    COMPILER_BARRIER(); 
  }

  void Stop() {
    #if 0
    //gettimeofday(&elapsed_time_, NULL);
    clock_gettime(CLOCK_MONOTONIC, &elapsed_time_);
    elapsed_time_.tv_sec  -= start_time_.tv_sec;
    //elapsed_time_.tv_usec -= start_time_.tv_usec;
    elapsed_time_.tv_nsec -= start_time_.tv_nsec;
    #endif
    COMPILER_BARRIER(); 
    stop_time_ = std::chrono::steady_clock::now();
    COMPILER_BARRIER(); 
    auto diff = stop_time_ - start_time_;
    elapsed_time += std::chrono::duration<double>(diff).count();
  }

  double Seconds() const {
    //return elapsed_time_.tv_sec + elapsed_time_.tv_usec/1e6 + elapsed_time_.tv_nsec/1e9;
    //return elapsed_time_.tv_sec + elapsed_time_.tv_nsec/1e9;
    return elapsed_time;
  }

  double Millisecs() const {
    //return 1000*elapsed_time_.tv_sec + elapsed_time_.tv_usec/1000 + elapsed_time_.tv_nsec/1e6;
    //return 1000*elapsed_time_.tv_sec + elapsed_time_.tv_nsec/1e6;
    return elapsed_time * 1000;
  }

  double Microsecs() const {
    //return 1e6*elapsed_time_.tv_sec + elapsed_time_.tv_usec + elapsed_time_.tv_nsec/1000;
    //return 1e6*elapsed_time_.tv_sec + elapsed_time_.tv_nsec/1000;
    return elapsed_time * 1000000;
  }
  
  double Nanosecs() const {
    //return 1e9*elapsed_time_.tv_sec + 1e3*elapsed_time_.tv_usec + elapsed_time_.tv_nsec;
    //return 1e9*elapsed_time_.tv_sec + elapsed_time_.tv_nsec;
    return elapsed_time * 1000000000;
  }

 private:
  //struct timeval start_time_;
  //struct timeval elapsed_time_;
  //struct timespec start_time_;
  //struct timespec elapsed_time_;
  std::chrono::steady_clock::time_point start_time_;
  std::chrono::steady_clock::time_point stop_time_;
  double elapsed_time;
};

// Times op's execution using the timer t
#define TIME_OP(t, op) { t.Start(); (op); t.Stop(); }

#endif  // TIMER_H_
