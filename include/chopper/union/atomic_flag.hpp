#pragma once

#include <vector>
#include <thread>
#include <atomic>
#include <functional>

// whether threads should keep waiting for work
std::atomic_flag stay_alive;

// whether there is work available
std::atomic_flag work_available;

// after a working session has started, whether all threads are done with their work
// and the main thread has processed has set the required flags
std::atomic_flag all_done;

class thread_provider {
    // pool of threads that stay alive for multiple working sessions
    std::vector<std::thread> pool;

    // whether the respective thread has started working
    std::vector<std::atomic_flag> started;

    // whether the respective thread is working at the moment
    std::vector<std::atomic_flag> busy;

    // work the respectiev thread is supposed to do
    std::vector<std::function<void()>> work;

    // results of a threads work
    std::vector<size_t> results;
    
    // number of worker threads
    size_t num_threads;

public:
    thread_provider(size_t num_threads_) :
        num_threads{num_threads_}
    {
        // threads should stay alive from now on
        stay_alive.test_and_set(std::memory_order_release);

        // there is not work available at the beginning
        work_available.clear(std::memory_order_release);

        // there is nothing running, all threads are in sync
        all_done.test_and_set(std::memory_order_release);

        for (size_t i = 0; i < num_threads; ++i) {
            // the thread has not yet started working
            started.emplace_back();
            started.back().clear(std::memory_order_release);

            // the thread is not busy at the beginning
            busy.emplace_back();
            busy.back().clear(std::memory_order_release);

            pool.emplace_back([this, i] () {
                // only keep alive as long as this flag is swithed on
                while (stay_alive.test(std::memory_order_relaxed)) {
                    // wait until there is something to do
                    work_available.wait(false, std::memory_order_acquire));

                    // notify the main thread that this thread has started working
                    started[i].test_and_set(std::memory_order_acquire);
                    started[i].notify_one();

                    // set this threads status to busy
                    busy[i].test_and_set(std::memory_order_acquire);

                    // do the work
                    work[i]();

                    // notify the main thread that this thread is done
                    busy[i].clear(std::memory_order_release);
                    busy[i].notify_one();

                    // wait until all threads are done with their work
                    all_done.wait(false, std::memory_order_acquire);
                }
            });
        }
    }

    size_t start_working_session(size_t size, bool last_session) {
        // set work
        work.clear();
        for (size_t i = 0; i < num_threads; ++i) {
            work.push_back([this, size, i] () {
                
                // do some random task
                size_t result = 0;
                for (size_t j = 0; j < size; ++j) {
                    if (j % 42 == 7) {
                        result += (j * 3) % 5;
                    }
                }
                results[i] = result;
            });
        }

        // signal start of asychronous working session
        all_done.test_and_set(std::memory_order_acq_rel);

        // notify all worker threads that there is work
        work_available.test_and_set(std::memory_order_acq_rel);
        work_available.notify_all();

        // make sure worker threads have started working
        for (auto & s : started) {
            s.wait(false, std::memory_order_acq_rel);
            s.clear();
        }

        // current work has started, there is no new work available
        work_available.clear(std::memory_order_release);

        // wait for all threads to finish
        for (auto & b : busy) {
            b.wait(true, std::memory_order_acquire);
        }

        // collect results
        size_t result = 0;
        for (size_t i = 0; i < num_threads; ++i) {
            result += results[i];
        }

        // order the threads to stop looking for work and die if its the last session
        if (last_session) {
            stay_alive.clear(std::memory_order_release);
        }

        // notify threads that the working session is over
        all_done.test_and_set(std::memory_order_acq_rel);
        all_done.notify_all();

        return result;
    }
};
