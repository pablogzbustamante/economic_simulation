#include "sim/thread_pool.hpp"

namespace sim {

ThreadPool::ThreadPool(int threads)
  : nthreads(threads > 0 ? threads : 1), stop(false), inflight(0) {
  workers.reserve((std::size_t)nthreads);
  for (int t = 0; t < nthreads; ++t) {
    workers.emplace_back([this, t] { worker_loop(t); });
  }
}

ThreadPool::~ThreadPool() {
  stop.store(true);
  cv.notify_all();
  for (auto& w : workers) {
    if (w.joinable()) w.join();
  }
}

int ThreadPool::size() const {
  return nthreads;
}

void ThreadPool::worker_loop(int) {
  for (;;) {
    std::function<void()> job;
    {
      std::unique_lock<std::mutex> lk(mtx);
      cv.wait(lk, [&] { return stop.load() || !q.empty(); });
      if (stop.load() && q.empty()) return;
      job = std::move(q.front());
      q.pop();
    }
    job();
    if (inflight.fetch_sub(1) == 1) {
      cv_done.notify_one();
    }
  }
}

void ThreadPool::parallel_for(int begin, int end, const std::function<void(int,int,int)>& fn) {
  const int n = end - begin;
  if (n <= 0) return;

  const int T = nthreads;
  const int chunk = (n + T - 1) / T;

  inflight.store(0);

  for (int tid = 0; tid < T; ++tid) {
    const int b = begin + tid * chunk;
    const int e = std::min(end, b + chunk);
    if (b >= e) continue;

    inflight.fetch_add(1);
    {
      std::lock_guard<std::mutex> lk(mtx);
      q.push([=] { fn(b, e, tid); });
    }
    cv.notify_one();
  }

  std::unique_lock<std::mutex> lk(mtx);
  cv_done.wait(lk, [&] { return inflight.load() == 0; });
}

}