#include "gen_poly.h"

static int func_cnt = 0;

int main()
{
    int split_num = 4;
    int unrol_num = 24;

    gen_init(split_num * unrol_num);
    for (int s = 1; s <= split_num; s++) {
        for (int u = 1; u <= unrol_num; u++) {
            gen_append_poly(s, u);
        }
    }
    return 0;
}

void gen_init(int func_num)
{
    /* header file */
    FILE *hf = fopen("dynamic_poly.h", "w");
    fprintf(hf, "typedef double (*PolyFunc)");
    fprintf(hf, "(double a[], double x, long degree);\n");
    fprintf(hf, "extern const PolyFunc func_arr[];\n");
    fclose(hf);

    /* c file */
    FILE *cf = fopen("dynamic_poly.c", "w");
    fprintf(cf, "#include \"dynamic_poly.h\"\n");
    fprintf(cf, "const PolyFunc func_arr[] = {\n");
    for (int i = 0; i < func_num; i++) {
        fprintf(cf, "poly_%d,\n", i);
    }
    fprintf(cf, "};\n");
    fclose(cf);
}

/* emit a balanced-tree sum of terms a[i+lo]*x_pow_lo .. a[i+hi]*x_pow_hi.
 * Balanced parenthesization keeps the add-reduction depth at O(log n)
 * instead of the O(n) left-associative chain, exposing more ILP.
 * NOTE: this intentionally changes the floating-point addition order, so
 * results are not guaranteed bit-identical to the chain version for a
 * general x. The benchmark uses x = -1 (integer terms), so it stays exact. */
static void emit_tree(FILE *cf, int lo, int hi)
{
    if (lo == hi) {
        fprintf(cf, "a[i + %d] * x_pow_%d", lo, lo);
        return;
    }
    int mid = (lo + hi) / 2;
    fprintf(cf, "(");
    emit_tree(cf, lo, mid);
    fprintf(cf, " + ");
    emit_tree(cf, mid + 1, hi);
    fprintf(cf, ")");
}

void gen_append_poly(int split_num, int unrol_num)
{
    /* header file */
    FILE *hf = fopen("dynamic_poly.h", "a");
    fprintf(hf, "double poly_%d(double a[], double x, long degree);\n",
            func_cnt);
    fclose(hf);

    /* c file */
    FILE *cf = fopen("dynamic_poly.c", "a");
    /* function name */
    fprintf(cf, "double poly_%d(double a[], double x, long degree)", func_cnt);
    fprintf(cf, "{ //(%d, %d)\n", split_num, unrol_num);
    fprintf(cf, "long i;\n");

    /* variables used for splitting */
    fprintf(cf, "double result_0 = a[0];\n");
    for (int k = 1; k < split_num; k++) {
        fprintf(cf, "double result_%d = 0;\n", k);
    }

    /* compute power to the x */
    int pw_times = split_num * unrol_num;
    fprintf(cf, "double xpwr = x;\n");
    fprintf(cf, "double x_pow_0 = 1;\n");
    fprintf(cf, "double x_pow_1 = x;\n");
    for (int k = 2; k <= pw_times; k++) {
        fprintf(cf, "double x_pow_%d = x_pow_%d * x;\n", k, k - 1);
    }

    /* main for loop */
    fprintf(cf, "for (i = 1; i <= degree - %d; i += %d) {\n", pw_times - 1,
            pw_times);
    for (int k = 0; k < split_num; k++) {
        int base = k * unrol_num;
        fprintf(cf, "result_%d += (", k);
        emit_tree(cf, base, base + unrol_num - 1);
        fprintf(cf, ") * xpwr;\n");
    }
    fprintf(cf, "xpwr = xpwr * x_pow_%d;\n", pw_times);
    fprintf(cf, "}\n");

    /* remain element */
    fprintf(cf, "for (; i <= degree; i++) {\n");
    fprintf(cf, "result_0 += a[i] * xpwr;\n");
    fprintf(cf, "xpwr = xpwr * x;\n}\n");

    /* sum the all splitting variables */
    fprintf(cf, "return result_0");
    for (int k = 1; k < split_num; k++) {
        fprintf(cf, " + result_%d", k);
    }
    fprintf(cf, ";\n}\n");
    fclose(cf);

    func_cnt++;
}
