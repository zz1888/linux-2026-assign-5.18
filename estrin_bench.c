#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_DEGREE    10000
#define DEGREE_STEP   10
#define TEST_TIMES    100
#define MEASURE_TIMES 100

/* ── reference: our best (split=1, unroll=10) ── */
double poly_split1_unroll10(double a[], double x, long degree)
{
    long i;
    double result = a[0];
    double xpwr = x;
    double x1  = x;
    double x2  = x1 * x;
    double x3  = x2 * x;
    double x4  = x3 * x;
    double x5  = x4 * x;
    double x6  = x5 * x;
    double x7  = x6 * x;
    double x8  = x7 * x;
    double x9  = x8 * x;
    double x10 = x9 * x;
    for (i = 1; i <= degree - 9; i += 10) {
        result += (a[i+0]
                 + a[i+1]*x1 + a[i+2]*x2 + a[i+3]*x3
                 + a[i+4]*x4 + a[i+5]*x5 + a[i+6]*x6
                 + a[i+7]*x7 + a[i+8]*x8 + a[i+9]*x9) * xpwr;
        xpwr *= x10;
    }
    for (; i <= degree; i++) {
        result += a[i] * xpwr;
        xpwr *= x;
    }
    return result;
}

/* ── tree reduction: split=1, unroll=10, explicit tree + no x_pow_0 mul ── */
double poly_tree10(double a[], double x, long degree)
{
    long i;
    double result = a[0];
    double xpwr = x;
    double x1  = x;
    double x2  = x1 * x;
    double x3  = x2 * x;
    double x4  = x3 * x;
    double x5  = x4 * x;
    double x6  = x5 * x;
    double x7  = x6 * x;
    double x8  = x7 * x;
    double x9  = x8 * x;
    double x10 = x9 * x;
    for (i = 1; i <= degree - 9; i += 10) {
        /* explicit balanced tree, a[i+0] needs no multiply */
        double t0 = a[i+0]      + a[i+1]*x1;
        double t1 = a[i+2]*x2   + a[i+3]*x3;
        double t2 = a[i+4]*x4   + a[i+5]*x5;
        double t3 = a[i+6]*x6   + a[i+7]*x7;
        double t4 = a[i+8]*x8   + a[i+9]*x9;
        double s0 = t0 + t1;
        double s1 = t2 + t3;
        double s2 = s0 + s1;
        double chunk = s2 + t4;
        result += chunk * xpwr;
        xpwr *= x10;
    }
    for (; i <= degree; i++) {
        result += a[i] * xpwr;
        xpwr *= x;
    }
    return result;
}

/* ── tree reduction with 2 outer accumulators (unroll=20) ── */
double poly_tree10_split2(double a[], double x, long degree)
{
    long i;
    double r0 = a[0], r1 = 0.0;
    double x1 = x;
    double x2 = x1*x, x3 = x2*x, x4 = x3*x, x5 = x4*x;
    double x6 = x5*x, x7 = x6*x, x8 = x7*x, x9 = x8*x, x10 = x9*x;
    double x20 = x10 * x10;
    double xp0 = x, xp1 = x*x10;   /* xpwr for block A and block B */
    for (i = 1; i <= degree - 19; i += 20) {
        /* block A: a[i..i+9] */
        double ta0 = a[i+0]    + a[i+1]*x1;
        double ta1 = a[i+2]*x2 + a[i+3]*x3;
        double ta2 = a[i+4]*x4 + a[i+5]*x5;
        double ta3 = a[i+6]*x6 + a[i+7]*x7;
        double ta4 = a[i+8]*x8 + a[i+9]*x9;
        r0 += (((ta0+ta1)+(ta2+ta3))+ta4) * xp0;
        /* block B: a[i+10..i+19] */
        double tb0 = a[i+10]    + a[i+11]*x1;
        double tb1 = a[i+12]*x2 + a[i+13]*x3;
        double tb2 = a[i+14]*x4 + a[i+15]*x5;
        double tb3 = a[i+16]*x6 + a[i+17]*x7;
        double tb4 = a[i+18]*x8 + a[i+19]*x9;
        r1 += (((tb0+tb1)+(tb2+tb3))+tb4) * xp1;
        xp0 *= x20;
        xp1 *= x20;
    }
    double result = r0 + r1;
    double xpwr = xp0;
    for (; i <= degree; i++) {
        result += a[i] * xpwr;
        xpwr *= x;
    }
    return result;
}

/* ── tree reduction: unroll=16 balanced tree (depth 4, like Estrin B=16) ── */
double poly_tree16(double a[], double x, long degree)
{
    long i;
    double result = a[0];
    double xpwr = x;
    double x1=x, x2=x1*x, x3=x2*x, x4=x3*x, x5=x4*x, x6=x5*x, x7=x6*x;
    double x8=x7*x, x9=x8*x, x10=x9*x, x11=x10*x, x12=x11*x, x13=x12*x;
    double x14=x13*x, x15=x14*x, x16=x15*x;
    for (i = 1; i <= degree - 15; i += 16) {
        double t0 = a[i+ 0]      + a[i+ 1]*x1;
        double t1 = a[i+ 2]*x2   + a[i+ 3]*x3;
        double t2 = a[i+ 4]*x4   + a[i+ 5]*x5;
        double t3 = a[i+ 6]*x6   + a[i+ 7]*x7;
        double t4 = a[i+ 8]*x8   + a[i+ 9]*x9;
        double t5 = a[i+10]*x10  + a[i+11]*x11;
        double t6 = a[i+12]*x12  + a[i+13]*x13;
        double t7 = a[i+14]*x14  + a[i+15]*x15;
        double s0 = t0 + t1, s1 = t2 + t3, s2 = t4 + t5, s3 = t6 + t7;
        double u0 = s0 + s1, u1 = s2 + s3;
        double chunk = u0 + u1;
        result += chunk * xpwr;
        xpwr *= x16;
    }
    for (; i <= degree; i++) {
        result += a[i] * xpwr;
        xpwr *= x;
    }
    return result;
}

/* ── Estrin block=4 ── */
double poly_estrin4(double a[], double x, long degree)
{
    double x2 = x * x;
    double x4 = x2 * x2;
    double result = 0.0;
    double xpwr = 1.0;
    long i;
    for (i = 0; i <= degree - 3; i += 4) {
        double p0 = a[i+0] + a[i+1] * x;
        double p1 = a[i+2] + a[i+3] * x;
        double r  = p0 + p1 * x2;
        result += r * xpwr;
        xpwr *= x4;
    }
    for (; i <= degree; i++) {
        result += a[i] * xpwr;
        xpwr *= x;
    }
    return result;
}

/* ── Estrin block=8 ── */
double poly_estrin8(double a[], double x, long degree)
{
    double x2 = x * x;
    double x4 = x2 * x2;
    double x8 = x4 * x4;
    double result = 0.0;
    double xpwr = 1.0;
    long i;
    for (i = 0; i <= degree - 7; i += 8) {
        double p0 = a[i+0] + a[i+1] * x;
        double p1 = a[i+2] + a[i+3] * x;
        double p2 = a[i+4] + a[i+5] * x;
        double p3 = a[i+6] + a[i+7] * x;
        double q0 = p0 + p1 * x2;
        double q1 = p2 + p3 * x2;
        double r  = q0 + q1 * x4;
        result += r * xpwr;
        xpwr *= x8;
    }
    for (; i <= degree; i++) {
        result += a[i] * xpwr;
        xpwr *= x;
    }
    return result;
}

/* ── Estrin block=16 ── */
double poly_estrin16(double a[], double x, long degree)
{
    double x2  = x  * x;
    double x4  = x2 * x2;
    double x8  = x4 * x4;
    double x16 = x8 * x8;
    double result = 0.0;
    double xpwr = 1.0;
    long i;
    for (i = 0; i <= degree - 15; i += 16) {
        double p0 = a[i+ 0] + a[i+ 1] * x;
        double p1 = a[i+ 2] + a[i+ 3] * x;
        double p2 = a[i+ 4] + a[i+ 5] * x;
        double p3 = a[i+ 6] + a[i+ 7] * x;
        double p4 = a[i+ 8] + a[i+ 9] * x;
        double p5 = a[i+10] + a[i+11] * x;
        double p6 = a[i+12] + a[i+13] * x;
        double p7 = a[i+14] + a[i+15] * x;
        double q0 = p0 + p1 * x2;
        double q1 = p2 + p3 * x2;
        double q2 = p4 + p5 * x2;
        double q3 = p6 + p7 * x2;
        double s0 = q0 + q1 * x4;
        double s1 = q2 + q3 * x4;
        double r  = s0 + s1 * x8;
        result += r * xpwr;
        xpwr *= x16;
    }
    for (; i <= degree; i++) {
        result += a[i] * xpwr;
        xpwr *= x;
    }
    return result;
}

/* ── Estrin block=8, two independent outer accumulators ── */
double poly_estrin8_split2(double a[], double x, long degree)
{
    double x2  = x  * x;
    double x4  = x2 * x2;
    double x8  = x4 * x4;
    double x16 = x8 * x8;
    double r0 = 0.0, r1 = 0.0;
    double xpwr0 = 1.0, xpwr1 = x8;
    long i;
    for (i = 0; i <= degree - 15; i += 16) {
        /* block A (elements i..i+7) */
        double p0 = a[i+0] + a[i+1]*x;
        double p1 = a[i+2] + a[i+3]*x;
        double p2 = a[i+4] + a[i+5]*x;
        double p3 = a[i+6] + a[i+7]*x;
        double q0 = p0 + p1*x2;
        double q1 = p2 + p3*x2;
        r0 += (q0 + q1*x4) * xpwr0;

        /* block B (elements i+8..i+15) */
        double p4 = a[i+ 8] + a[i+ 9]*x;
        double p5 = a[i+10] + a[i+11]*x;
        double p6 = a[i+12] + a[i+13]*x;
        double p7 = a[i+14] + a[i+15]*x;
        double q2 = p4 + p5*x2;
        double q3 = p6 + p7*x2;
        r1 += (q2 + q3*x4) * xpwr1;

        xpwr0 *= x16;
        xpwr1 *= x16;
    }
    double result = r0 + r1;
    double xpwr = (i == 0) ? 1.0 : xpwr0;
    for (; i <= degree; i++) {
        result += a[i] * xpwr;
        xpwr *= x;
    }
    return result;
}

/* ── Estrin block=32 ── */
double poly_estrin32(double a[], double x, long degree)
{
    double x2  = x   * x;
    double x4  = x2  * x2;
    double x8  = x4  * x4;
    double x16 = x8  * x8;
    double x32 = x16 * x16;
    double result = 0.0;
    double xpwr = 1.0;
    long i;
    for (i = 0; i <= degree - 31; i += 32) {
        double p0  = a[i+ 0] + a[i+ 1]*x;
        double p1  = a[i+ 2] + a[i+ 3]*x;
        double p2  = a[i+ 4] + a[i+ 5]*x;
        double p3  = a[i+ 6] + a[i+ 7]*x;
        double p4  = a[i+ 8] + a[i+ 9]*x;
        double p5  = a[i+10] + a[i+11]*x;
        double p6  = a[i+12] + a[i+13]*x;
        double p7  = a[i+14] + a[i+15]*x;
        double p8  = a[i+16] + a[i+17]*x;
        double p9  = a[i+18] + a[i+19]*x;
        double p10 = a[i+20] + a[i+21]*x;
        double p11 = a[i+22] + a[i+23]*x;
        double p12 = a[i+24] + a[i+25]*x;
        double p13 = a[i+26] + a[i+27]*x;
        double p14 = a[i+28] + a[i+29]*x;
        double p15 = a[i+30] + a[i+31]*x;
        double q0  = p0  + p1  * x2;
        double q1  = p2  + p3  * x2;
        double q2  = p4  + p5  * x2;
        double q3  = p6  + p7  * x2;
        double q4  = p8  + p9  * x2;
        double q5  = p10 + p11 * x2;
        double q6  = p12 + p13 * x2;
        double q7  = p14 + p15 * x2;
        double s0  = q0  + q1  * x4;
        double s1  = q2  + q3  * x4;
        double s2  = q4  + q5  * x4;
        double s3  = q6  + q7  * x4;
        double t0  = s0  + s1  * x8;
        double t1  = s2  + s3  * x8;
        double r   = t0  + t1  * x16;
        result += r * xpwr;
        xpwr *= x32;
    }
    for (; i <= degree; i++) {
        result += a[i] * xpwr;
        xpwr *= x;
    }
    return result;
}

/* ──────────────────────────────────────── */

double tvgetf()
{
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

long read_cpu_freq()
{
    FILE *f = fopen("/sys/devices/system/cpu/cpu0/cpufreq/cpuinfo_max_freq", "r");
    char buf[32] = {0};
    if (fgets(buf, sizeof(buf), f) == NULL) { fclose(f); return -1; }
    fclose(f);
    return atol(buf) * 1000;
}

/* measure CPE at a single degree, same method as test_poly.c */
double bench_once(double (*fn)(double*, double, long),
                  double a[], long degree, long cpu_freq)
{
    double cyc[TEST_TIMES];
    volatile double ans = 0;
    for (int t = 0; t < TEST_TIMES; t++) {
        double t1 = tvgetf();
        for (int j = 0; j < MEASURE_TIMES; j++)
            ans += fn(a, -1.0, degree);
        double t2 = tvgetf();
        cyc[t] = (t2 - t1) * cpu_freq / MEASURE_TIMES;
    }
    for (int i = 0; i < TEST_TIMES-1; i++)
        for (int j = i+1; j < TEST_TIMES; j++)
            if (cyc[j] < cyc[i]) { double tmp=cyc[i]; cyc[i]=cyc[j]; cyc[j]=tmp; }
    int trim = TEST_TIMES / 20;
    double sum = 0;
    for (int i = trim; i < TEST_TIMES-trim; i++) sum += cyc[i];
    return sum / (TEST_TIMES - 2*trim) / degree;
}

/* average CPE over degree DEGREE_STEP .. MAX_DEGREE-DEGREE_STEP (same as default_test) */
double bench(double (*fn)(double*, double, long), long cpu_freq)
{
    double a[MAX_DEGREE + 1];
    for (int i = 0; i <= MAX_DEGREE; i++) a[i] = (double)i;

    double total = 0.0;
    int count = 0;
    for (int d = DEGREE_STEP; d < MAX_DEGREE; d += DEGREE_STEP) {
        total += bench_once(fn, a, d, cpu_freq);
        count++;
    }
    return total / count;
}

int main()
{
    long cpu_freq = read_cpu_freq();
    printf("CPU max freq = %ld Hz\n\n", cpu_freq);

    /* correctness check at degree=MAX_DEGREE */
    double a[MAX_DEGREE + 1];
    for (int i = 0; i <= MAX_DEGREE; i++) a[i] = (double)i;
    double expected = (double)MAX_DEGREE / 2;
    double v;
    v = poly_split1_unroll10(a, -1, MAX_DEGREE);
    printf("split1_unroll10  ans=%.1f %s\n", v, v==expected?"OK":"WRONG");
    v = poly_estrin4(a, -1, MAX_DEGREE);
    printf("estrin4          ans=%.1f %s\n", v, v==expected?"OK":"WRONG");
    v = poly_estrin8(a, -1, MAX_DEGREE);
    printf("estrin8          ans=%.1f %s\n", v, v==expected?"OK":"WRONG");
    v = poly_estrin16(a, -1, MAX_DEGREE);
    printf("estrin16         ans=%.1f %s\n", v, v==expected?"OK":"WRONG");
    v = poly_estrin8_split2(a, -1, MAX_DEGREE);
    printf("estrin8_split2   ans=%.1f %s\n", v, v==expected?"OK":"WRONG");
    v = poly_tree10(a, -1, MAX_DEGREE);
    printf("tree10           ans=%.1f %s\n", v, v==expected?"OK":"WRONG");
    v = poly_tree10_split2(a, -1, MAX_DEGREE);
    printf("tree10_split2    ans=%.1f %s\n", v, v==expected?"OK":"WRONG");
    v = poly_tree16(a, -1, MAX_DEGREE);
    printf("tree16           ans=%.1f %s\n", v, v==expected?"OK":"WRONG");
    printf("\n");

    double pa[MAX_DEGREE + 1];
    for (int i = 0; i <= MAX_DEGREE; i++) pa[i] = (double)i;

    /* focused comparison: chain vs tree vs Estrin */
    static const long probes[] = {100, 500, 1000, 2000, 5000, 10000};
    int np = sizeof(probes) / sizeof(probes[0]);

    printf("%-8s  %-12s  %-12s  %-12s  %-12s\n",
           "degree", "split_chain", "tree10", "tree16", "estrin_B16");
    printf("%-8s  %-12s  %-12s  %-12s  %-12s\n",
           "------", "-----------", "----------", "----------", "----------");
    for (int pi = 0; pi < np; pi++) {
        long d = probes[pi];
        printf("%-8ld  %-12.6f  %-12.6f  %-12.6f  %-12.6f\n", d,
               bench_once(poly_split1_unroll10, pa, d, cpu_freq),
               bench_once(poly_tree10,          pa, d, cpu_freq),
               bench_once(poly_tree16,          pa, d, cpu_freq),
               bench_once(poly_estrin16,        pa, d, cpu_freq));
    }
    printf("\n");

    /* plot data: dense sweep degree 200~10000 step 200 */
    printf("\nWriting plot data to estrin_plot.dat ...\n");
    FILE *dat = fopen("estrin_plot.dat", "w");
    fprintf(dat, "# degree split1_u10 estrin_B8 estrin_B16 estrin_B32\n");
    for (long d = 200; d <= MAX_DEGREE; d += 200) {
        double c0 = bench_once(poly_split1_unroll10, pa, d, cpu_freq);
        double c1 = bench_once(poly_estrin8,         pa, d, cpu_freq);
        double c2 = bench_once(poly_estrin16,        pa, d, cpu_freq);
        double c3 = bench_once(poly_estrin32,        pa, d, cpu_freq);
        fprintf(dat, "%ld %.6f %.6f %.6f %.6f\n", d, c0, c1, c2, c3);
        printf("  degree=%5ld  split=%.3f  B8=%.3f  B16=%.3f  B32=%.3f\n",
               d, c0, c1, c2, c3);
        fflush(stdout);
    }
    fclose(dat);
    printf("Done. Run: gnuplot estrin_plot.gp\n");
    return 0;
}
