#pragma once
#include <vector>
#include <thread>
#include <functional>
#include <condition_variable>
#include <mutex>
#include <queue>
#include <atomic>

namespace sim {

class ThreadPool {
public:
  explicit ThreadPool(int threads);
  ~ThreadPool();

  int size() const;

  void parallel_for(int begin, int end, const std::function<void(int,int,int)>& fn);

private:
  void worker_loop(int tid);

  int nthreads;
  std::vector<std::thread> workers;

  std::mutex mtx;
  std::condition_variable cv;
  std::queue<std::function<void()>> q;

  std::atomic<bool> stop;
  std::atomic<int> inflight;

  std::condition_variable cv_done;
};

}