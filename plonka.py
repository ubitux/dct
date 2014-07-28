import sys, random
import numpy as np
from math import cos, sin, sqrt, pi

np.set_printoptions(precision=3, linewidth=999)

# HACK: Overriden by gen_c.py
sqrt2 = sqrt(2)
sqrt1_2 = 1./sqrt(2)
cos_k_pi_n = lambda k, n: cos(k*pi/n)
sin_k_pi_n = lambda k, n: sin(k*pi/n)

# a b
# c d
quad = lambda a_b, c_d: np.vstack((np.hstack(a_b),
                                   np.hstack(c_d)))

# a,b,c,d matrices, x scalar, output:
#   a   b
#     x
#   c   d
def quint(a, b, c, d, x):
    ah, aw = a.shape
    bh, bw = b.shape
    ch, cw = c.shape
    dh, dw = d.shape
    assert ah == bh and ch == dh and \
           aw == cw and bw == dw
    w = aw + 1 + bw
    h = ah + 1 + ch
    m = np.zeros((h, w))
    m[ah,aw] = x
    m[0:ah,  0:aw]  = a
    m[0:bh,  aw+1:] = b
    m[ah+1:, 0:cw]  = c
    m[bh+1:, cw+1:] = d
    return m

# a 0
# 0 b
def diag(a, b=None):
    if b is None:
        return a * I(len(a))
    if isinstance(a, int) or isinstance(a, float):
        a = np.array([[a]])
    if isinstance(b, int) or isinstance(b, float):
        b = np.array([[b]])
    return quad((a, np.zeros((a.shape[0], b.shape[1]))),
                (np.zeros((b.shape[0], a.shape[1])), b))

def npf(f, n):
    m = []
    for j in range(n):
        m.append([f(j, k) for k in range(n)])
    return np.array(m)

def E(x, n):
    return sqrt1_2 if x in (0, n) else 1.

def C_I(np1):
    n = np1 - 1
    f = lambda j, k: sqrt(2./n) * E(j, n) * E(k, n) * cos_k_pi_n(j*k, n)
    return npf(f, np1)

def C_II(n):
    f = lambda j, k: sqrt(2./n) * E(j, n) * cos_k_pi_n(j*(2*k+1), 2*n)
    return npf(f, n)

def C_III(n):
    return C_II(n).T

def C_IV(n):
    f = lambda j, k: sqrt(2./n) * cos_k_pi_n((2*j+1)*(2*k+1), 4*n)
    return npf(f, n)

def S_I(nm1):
    n = nm1 + 1
    f = lambda j, k: sqrt(2./n) * sin_k_pi_n((j+1)*(k+1), n)
    return npf(f, nm1)

I = lambda n: np.identity(n, dtype=int)
J = lambda n: np.fliplr(I(n))
D = lambda n: diag(np.array([(-1)**k for k in range(n)]))

def permute_m(n):
    '''Permutation matrix (not transposed), defined as P in the paper'''
    #assert n >= 4
    m = np.zeros((n, n))
    col = 0
    for i in range(n):
        m[i, col] = 1
        col += 2
        if col >= n:
            col = 1
    return m

def add_m(n, b, modified):
    '''Addition matrices'''
    if not modified:
        '''Addition matrix, defined as A in the paper'''
        n1 = n / 2
        if b == 0:
            return I(n)
        elif b == 1:
            isq2_n1m1 = I(n1-1) * sqrt1_2
            m = sqrt1_2 * quad((I(n1-1),  I(n1-1)),
                               (I(n1-1), -I(n1-1)))
            op1 = diag(diag(1, m), -1)
            op2 = diag(I(n1), D(n1).dot(J(n1)))
            return op1.dot(op2)
    else:
        '''Second addition matrix, defined as A~ in the paper'''
        n1 = n / 2
        if b == 1:
            return I(n + 1)
        elif b == 0:
            op1 = diag(I(n1), J(n1))
            op2m = quint(I(n1-1),  J(n1-1),
                         J(n1-1), -I(n1-1),
                         sqrt2)
            op2 = diag(1, sqrt1_2 * op2m)
            op3 = diag(I(n1+1), (-1)**n1 * D(n1-1))
            return op1.dot(op2).dot(op3)
        elif b == -1:
            return diag(D(n1), I(n1 - 1))
    assert False

def twiddle_m(n, b, modified):
    '''Twiddle matrices'''
    if not modified:
        '''Twiddle matrix, defined as T in the paper'''
        n1 = n / 2
        if b == 0:
            return sqrt1_2 * quad((I(n1),  J(n1)),
                                  (I(n1), -J(n1)))
        elif b == 1:
            m00 = diag(c(n, n1))
            m01 = diag(s(n, n1)).dot(J(n1))
            m10 = (-J(n1)).dot(diag(s(n, n1)))
            m11 = diag(J(n1).dot(c(n, n1)))
            op1 = diag(I(n1), D(n1))
            op2 = quad((m00, m01), (m10, m11))
            return op1.dot(op2)
    else:
        '''Second twiddle matrix, defined as T~ in the paper'''
        n1 = n / 2
        if b == 1:
            return sqrt1_2 * quint(I(n1),  J(n1),
                                   I(n1), -J(n1),
                                   sqrt2)
        elif b == 0:
            op1 = diag(I(n1+1), D(n1-1))
            m00 = diag(ct(n, n1-1))
            m01 = diag(st(n, n1-1)).dot(J(n1-1))
            m10 = (-J(n1-1)).dot(diag(st(n, n1-1)))
            m11 = diag(J(n1-1).dot(ct(n, n1-1)))
            m = quint(m00, m01, m10, m11, 1)
            op2 = diag(1, m)
            return op1.dot(op2)
        elif b == -1:
            return diag(J(n1), I(n1-1)).dot(twiddle_m(n-1, b=1, modified=1))
    assert False

c = lambda n, n1: [cos_k_pi_n(2*k+1, 4*n) for k in range(n1)]
s = lambda n, n1: [sin_k_pi_n(2*k+1, 4*n) for k in range(n1)]
# XXX: page 6/37, paper has a mistake (range is set from 1 to n1-1 instead of 1
# to n1)
ct = lambda n, n1m1: [cos_k_pi_n(k, 2*n) for k in range(1, n1m1+1)]
st = lambda n, n1m1: [sin_k_pi_n(k, 2*n) for k in range(1, n1m1+1)]

cosII_neq2_mat = np.array([[1,  1],
                           [1, -1]])
cosIV_neq2_mat = sqrt2 * C_IV(2)
cosI_neq2_mat  = 1./2 * (np.array([[1,1,0],[0,0,sqrt2],[1,-1,0]]).dot(
                         np.array([[1,0,1],[0,sqrt2,0],[1,0,-1]])))
cosIII_neq2_mat = sqrt1_2 * np.array([[1,  1],
                                      [1, -1]])
sinI_neq2_mat = I(1)

def tfm_run(name, x, first_call=True):

    tfm_props = {
        'cosII':  (cosII_neq2_mat,   0,  0, 'cosII',  'cosIV',  0, 0),
        'cosIV':  (cosIV_neq2_mat,   0,  1, 'cosII',  'cosII',  0, 0),
        'cosI':   (cosI_neq2_mat,   -1,  1, 'cosI',   'cosIII', 1, 1),
        'cosIII': (cosIII_neq2_mat,  0,  0, 'cosI',   'sinI',   1, 1),
        'sinI':   (sinI_neq2_mat,    1, -1, 'cosIII', 'sinI',   0, 1),
    }

    neq2_mat, n_delay, b, v1_tfm, v2_tfm, hvec_size_add, modified_matrix = tfm_props[name]
    n_orig = len(x)
    n = n_orig + n_delay
    n1 = n / 2
    if n == 2:
        y = neq2_mat.dot(x)
    elif n >= 4:
        u = sqrt2 * twiddle_m(n, b, modified_matrix).dot(x)
        v1 = tfm_run(v1_tfm, u[:n1+hvec_size_add], False)
        v2 = tfm_run(v2_tfm, u[n1+hvec_size_add:], False)
        w = add_m(n, b, modified_matrix).dot(np.hstack((v1, v2)))
        y = permute_m(n_orig).T.dot(w)
    else:
        assert False
    return 1./sqrt(n>>modified_matrix) * y if first_call else y

if __name__ == '__main__':
    tests = {
        'cosI':   (C_I,   1),
        'cosII':  (C_II,  0),
        'cosIII': (C_III, 0),
        'cosIV':  (C_IV,  0),
        'sinI':   (S_I,  -1),
    }

    if len(sys.argv) != 3:
        print 'Usage: %s <%s> <N>' % (sys.argv[0], '|'.join(tests.keys()))
        sys.exit(0)

    name = sys.argv[1]
    t = int(sys.argv[2])

    n = 1<<t
    ref, delay = tests[name]

    samples = [random.randint(0x00, 0xff) for x in range(n + delay)]

    fdct_ref = ref(n + delay).dot(samples)
    fdct_out = tfm_run(name, samples)
    diffs = fdct_ref - fdct_out
    ok = True
    for diff in diffs:
        if abs(diff) > 1e-10:
            ok = False
            print
            print 'ref', fdct_ref
            print 'out', fdct_out
            print 'cmp', diffs
            sys.exit(1)
