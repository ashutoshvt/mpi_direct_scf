double wall_time() {
    static bool first_call = true;
    static double start_time;

    struct timeval tv;
    gettimeofday(&tv,0);
    double now = tv.tv_sec + 1e-6*tv.tv_usec;

    if (first_call) {
        first_call = false;
        start_time = now;
    }
    return now - start_time;
}

