#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/********* generated code snippet *********/

#define N %N%

static inline void fdct_1d(float *dst, const float *src,
                           int dst_stridea, int dst_strideb,
                           int src_stridea, int src_strideb)
{
    int i;

    for (i = 0; i < N; i++) {
%CODE_FDCT%
        dst += dst_strideb;
        src += src_strideb;
    }
}

static void fdct(float *dst, const float *src)
{
    float tmp[N*N];
    fdct_1d(tmp, src, 1, N, 1, N);
    fdct_1d(dst, tmp, N, 1, N, 1);
}

static inline void idct_1d(float *dst, const float *src,
                           int dst_stridea, int dst_strideb,
                           int src_stridea, int src_strideb)
{
    int i;

    for (i = 0; i < N; i++) {
%CODE_IDCT%
        dst += dst_strideb;
        src += src_strideb;
    }
}

static void idct(float *dst, const float *src)
{
    float tmp[N*N];
    idct_1d(tmp, src, 1, N, 1, N);
    idct_1d(dst, tmp, N, 1, N, 1);
}


/********* slow reference dct code ********/

static float dct_matrix    [N*N];
static float dct_trp_matrix[N*N];

void init_dct(void)
{
    int i, j;

    // dct matrix
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            if (i == 0)
                dct_matrix[i*N + j] = 1 / sqrt(N);
            else
                dct_matrix[i*N + j] = sqrt(2./N) * cos(((2*j+1)*i*M_PI) / (2*N));

    // dct matrix transposed
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            dct_trp_matrix[i*N + j] = dct_matrix[j*N + i];
}

static void dct_1d_ref(float *dst, const float *src,
                       int stridea, int strideb,
                       const float *matrix)
{
    int x;
    for (x = 0; x < N; x++) {
        int i, j;
        for (j = 0; j < N; j++) {
            float sum = 0.;
            for (i = 0; i < N; i++)
                sum += matrix[j*N + i] * src[i*stridea];
            dst[j*stridea] = sum;
        }
        dst += strideb;
        src += strideb;
    }
}

static void fdct_ref(float *dst, const float *src)
{
    float tmp[N*N];
    dct_1d_ref(tmp, src, 1, N, dct_matrix);
    dct_1d_ref(dst, tmp, N, 1, dct_matrix);
}

static void idct_ref(float *dst, const float *src)
{
    float tmp[N*N];
    dct_1d_ref(tmp, src, 1, N, dct_trp_matrix);
    dct_1d_ref(dst, tmp, N, 1, dct_trp_matrix);
}

/********** test **************************/

static int check_output(const char *name, const float *ref, const float *out)
{
    int i, ret = 0;

    for (i = 0; i < N*N; i++) {
        const int ok = fabs(ref[i] - out[i]) < 0.0005;
        if (!ok) {
            printf("%s ref:%9.3f out:%9.3f diff:%9.3f\n",
                   name, ref[i], out[i], ref[i] - out[i]);
            ret = -1;
        }
    }
    return ret;
}

int main()
{
    int i;
    float src[N*N];
    float ref_fdct[N*N], ref_idct[N*N];
    float out_fdct[N*N], out_idct[N*N];

    for (i = 0; i < N*N; i++)
        src[i] = random() % 256;

    init_dct();

    fdct_ref(ref_fdct, src);
    fdct(out_fdct, src);
    if (check_output("FDCT", ref_fdct, out_fdct) < 0)
        return 1;

    idct_ref(ref_idct, ref_fdct);
    idct(out_idct, out_fdct);
    if (check_output("IDCT", ref_idct, out_idct) < 0)
        return 1;

    return 0;
}
