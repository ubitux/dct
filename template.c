#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/********* generated code snippet *********/

#define N %N%

#define COS_1PI8  1.3065629648763766 /* sqrt(2)*cos(  pi/8) */
#define COS_3PI8  0.5411961001461971 /* sqrt(2)*cos(3*pi/8) */
#define COS_9PI8 -1.3065629648763768 /* sqrt(2)*cos(9*pi/8) */

static void fdct_1d(float *dst, const float *src,
                    int stridea, int strideb)
{
    int i;

    for (i = 0; i < N; i++) {
%CODE%
        dst += strideb;
        src += strideb;
    }
}

static void fdct(float *dst, const float *src)
{
    float tmp[N*N];
    fdct_1d(tmp, src, 1, N);
    fdct_1d(dst, tmp, N, 1);
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

static void fdct_1d_ref(float *dst, const float *src,
                        int stridea, int strideb)
{
    int x;
    for (x = 0; x < N; x++) {
        int i, j;
        for (j = 0; j < N; j++) {
            float sum = 0.;
            for (i = 0; i < N; i++)
                sum += dct_matrix[j*N + i] * src[i*stridea];
            dst[j*stridea] = sum;
        }
        dst += strideb;
        src += strideb;
    }
}

static void fdct_ref(float *dst, const float *src)
{
    float tmp[N*N];
    fdct_1d_ref(tmp, src, 1, N);
    fdct_1d_ref(dst, tmp, N, 1);
}

/********** test **************************/

int main()
{
    int i, ret = 0;
    float src[N*N];
    float dstref[N*N];
    float dstdct[N*N];

    for (i = 0; i < N*N; i++)
        src[i] = random() % 256;

    init_dct();
    fdct_ref(dstref, src);

    fdct(dstdct, src);

    for (i = 0; i < N*N; i++) {
        const int ok = fabs(dstref[i] - dstdct[i]) < 0.0001;
        if (!ok) {
            printf("src:%9.3f ref:%9.3f out:%9.3f diff:%9.3f\n",
                   src[i], dstref[i], dstdct[i], dstref[i]-dstdct[i]);
            ret = 1;
        }
    }

    return ret;
}
