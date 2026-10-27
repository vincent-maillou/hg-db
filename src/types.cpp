// types.cpp - Implementation of utility functions
#include "hypergraph_reorder/types.hpp"

#include <chrono>

namespace hypergraph_reorder {

Timer::Timer() { start_time = get_time(); }

void Timer::reset() { start_time = get_time(); }

double Timer::elapsed_ms() const { return (get_time() - start_time) * 1000.0; }

double Timer::get_time() {
  // steady_clock: monotonic (immune to NTP steps / VM migrations).
  // high_resolution_clock aliases system_clock on libstdc++ (is_steady == 0),
  // which makes timings corruptible by wall-clock adjustments (spec D8).
  using namespace std::chrono;
  return duration_cast<duration<double>>(steady_clock::now().time_since_epoch())
      .count();
}

}  // namespace hypergraph_reorder
