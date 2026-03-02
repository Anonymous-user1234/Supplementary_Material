#include<sys/time.h>

#define OUTPUT_TIME() {std::cout << "\t" << getCurrentSystemTime() << std::endl << std::endl;}

#define OUTPUT_TT_TIME() {std::cout << "\t\t" << getCurrentSystemTime() << std::endl << std::endl;}



/* convert timeval to miliseconds */
#define TIMEVAL2F(stamp) ((stamp).tv_sec * 1000.0 + (stamp).tv_usec / 1000.0)

/* get timestamp to the precision of miliseconds since the program starts */
inline double get_timestamp() {
    static double __init_stamp = -1;
    static struct timeval __cur_time;
    if (-1 == __init_stamp) {
        gettimeofday(&__cur_time, NULL);
        __init_stamp = TIMEVAL2F(__cur_time);
    }
    
    gettimeofday(&__cur_time, NULL);
    return ((TIMEVAL2F(__cur_time) - __init_stamp) / 1000.0);
}

/* print msg with timestamp */
#define PRINTF_STAMP(format, ...) \
do { \
    flockfile(stdout); \
    printf("%12.2f - ", get_timestamp()); \
    printf(format, ##__VA_ARGS__); \
    fflush(stdout); \
    funlockfile(stdout); \
} while(0)


inline std::string getCurrentSystemTime() {
    auto tt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    struct tm* ptm = localtime(&tt);
    char date[60] = { 0 };
    sprintf(date, "%d-%02d-%02d-%02d:%02d:%02d", (int)ptm->tm_year + 1900, (int)ptm->tm_mon + 1, (int)ptm->tm_mday, (int)ptm->tm_hour, (int)ptm->tm_min, (int)ptm->tm_sec);
    return std::string(date);
}


class Timer {
   
private:
    double start_time, end_time, durable_time, overall_time;
    bool initialized = false;

public:

    inline void clear(){
        overall_time = 0.;
        initialized = true;
    }

    inline void start() { start_time = get_timestamp(); }
    inline void end() {
        end_time = get_timestamp();
        durable_time = end_time - start_time;
        if (initialized) overall_time += durable_time;
        else {
            initialized = true;
            overall_time = durable_time;
        }
    }    
    inline double used_time() { return get_timestamp() - start_time; }

    inline void show_time() { std::cout << get_timestamp() - start_time << " sec. used\n"; }

    inline double duration() { return durable_time; }
    inline double overall() { return overall_time; }
};
