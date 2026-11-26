/*
 * matmul_bench.c
 * Measures how loop order and blocking change the runtime of N x N matrix
 * multiplication (C = A * B, double precision, row-major storage).
 *
 * Versions:
 *   ijk     : textbook order; inner loop walks down a column of B
 *   ikj     : reordered; inner loop walks along rows of B and C
 *   blocked : tiled version (block size BS), loop order ii-kk-jj / i-k-j
 *
 * Build:   gcc -O2 -o matmul_bench matmul_bench.c -lm
 * Run:     ./matmul_bench                 (runs=10, sizes 64..1024, block 32)
 *          ./matmul_bench 10 1024 32      (runs, maxN, block)
 *          ./matmul_bench 10 1024 blocks  (block-size sweep at N = maxN)
 * Output:  CSV on stdout:  version,N,block,run,seconds
 *          Information lines start with '#'. Save output with:  > results.csv
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef _WIN32
#include <windows.h>
/* Windows: high-resolution performance counter */
static double now_seconds(void) {
    LARGE_INTEGER freq, cnt;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&cnt);
    return (double)cnt.QuadPart / (double)freq.QuadPart;
}
#else
#include <time.h>
/* Linux / Mac: monotonic clock */
static double now_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}
#endif

static void mm_ijk(const double *A, const double *B, double *C, int n) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++)
                sum += A[(size_t)i * n + k] * B[(size_t)k * n + j];
            C[(size_t)i * n + j] = sum;
        }
}

static void mm_ikj(const double *A, const double *B, double *C, int n) {
    memset(C, 0, (size_t)n * n * sizeof(double));
    for (int i = 0; i < n; i++)
        for (int k = 0; k < n; k++) {
            double a = A[(size_t)i * n + k];
            for (int j = 0; j < n; j++)
                C[(size_t)i * n + j] += a * B[(size_t)k * n + j];
        }
}

static void mm_blocked(const double *A, const double *B, double *C, int n, int bs) {
    memset(C, 0, (size_t)n * n * sizeof(double));
    for (int ii = 0; ii < n; ii += bs)
        for (int kk = 0; kk < n; kk += bs)
            for (int jj = 0; jj < n; jj += bs) {
                int imax = ii + bs < n ? ii + bs : n;
                int kmax = kk + bs < n ? kk + bs : n;
                int jmax = jj + bs < n ? jj + bs : n;
                for (int i = ii; i < imax; i++)
                    for (int k = kk; k < kmax; k++) {
                        double a = A[(size_t)i * n + k];
                        for (int j = jj; j < jmax; j++)
                            C[(size_t)i * n + j] += a * B[(size_t)k * n + j];
                    }
            }
}

static double max_abs_diff(const double *X, const double *Y, size_t len) {
    double m = 0.0;
    for (size_t i = 0; i < len; i++) {
        double d = fabs(X[i] - Y[i]);
        if (d > m) m = d;
    }
    return m;
}

static void fill(double *M, size_t len, unsigned seed) {
    unsigned s = seed;
    for (size_t i = 0; i < len; i++) {
        s = s * 1664525u + 1013904223u;            /* simple repeatable generator */
        M[i] = (double)(s >> 8) / 16777216.0;      /* values in [0,1) */
    }
}

/* run one version `runs` times and print one CSV line per run */
static void run_version(const char *name, int n, int bs, int runs,
                        const double *A, const double *B, double *C, int which) {
    for (int r = 0; r < runs; r++) {
        double t0 = now_seconds();
        if (which == 0)      mm_ijk(A, B, C, n);
        else if (which == 1) mm_ikj(A, B, C, n);
        else                 mm_blocked(A, B, C, n, bs);
        double t1 = now_seconds();
        printf("%s,%d,%d,%d,%.6f\n", name, n, which == 2 ? bs : 0, r + 1, t1 - t0);
        fflush(stdout);
    }
}

int main(int argc, char **argv) {
    int runs = argc > 1 ? atoi(argv[1]) : 10;
    int maxn = argc > 2 ? atoi(argv[2]) : 1024;
    int sweep = (argc > 3 && strcmp(argv[3], "blocks") == 0);
    int bs = (argc > 3 && !sweep) ? atoi(argv[3]) : 32;
    if (runs < 1 || maxn < 8 || bs < 1) { fprintf(stderr, "bad arguments\n"); return 1; }

    printf("# matmul_bench: runs=%d maxN=%d mode=%s block=%d\n", runs, maxn,
           sweep ? "block-sweep" : "size-sweep", bs);
    printf("version,N,block,run,seconds\n");

    int sizes_all[] = {64, 128, 256, 512, 1024, 2048};
    int blocks_all[] = {8, 16, 32, 64, 128, 256};

    if (!sweep) {
        for (size_t s = 0; s < sizeof(sizes_all) / sizeof(int); s++) {
            int n = sizes_all[s];
            if (n > maxn) break;
            size_t len = (size_t)n * n;
            double *A = malloc(len * sizeof(double)), *B = malloc(len * sizeof(double));
            double *C = malloc(len * sizeof(double)), *R = malloc(len * sizeof(double));
            if (!A || !B || !C || !R) { fprintf(stderr, "out of memory at N=%d\n", n); return 1; }
            fill(A, len, 1u); fill(B, len, 2u);

            /* correctness check: the other versions must match ijk */
            mm_ijk(A, B, R, n);
            mm_ikj(A, B, C, n);
            double d1 = max_abs_diff(R, C, len);
            mm_blocked(A, B, C, n, bs);
            double d2 = max_abs_diff(R, C, len);
            printf("# N=%d correctness: max|ijk-ikj|=%.3e max|ijk-blocked|=%.3e\n", n, d1, d2);
            if (d1 > 1e-9 || d2 > 1e-9) { fprintf(stderr, "RESULT MISMATCH at N=%d\n", n); return 2; }

            run_version("ijk", n, bs, runs, A, B, C, 0);
            run_version("ikj", n, bs, runs, A, B, C, 1);
            run_version("blocked", n, bs, runs, A, B, C, 2);
            free(A); free(B); free(C); free(R);
        }
    } else {
        int n = maxn;
        size_t len = (size_t)n * n;
        double *A = malloc(len * sizeof(double)), *B = malloc(len * sizeof(double));
        double *C = malloc(len * sizeof(double));
        if (!A || !B || !C) { fprintf(stderr, "out of memory\n"); return 1; }
        fill(A, len, 1u); fill(B, len, 2u);
        for (size_t b = 0; b < sizeof(blocks_all) / sizeof(int); b++)
            run_version("blocked", n, blocks_all[b], runs, A, B, C, 2);
        free(A); free(B); free(C);
    }
    return 0;
}
