import sys, random
import numpy as np
from math import cos, sin, sqrt, pi

np.set_printoptions(precision=3, linewidth=999)

#def dct_II_split_radix_factor(s):
#    b = [0]
#    for i in range(s):
#        b += b
#        b[-1] ^= 1
#    return b

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
    return 1./sqrt(2) if x in (0, n) else 1.

def C_I(np1):
    n = np1 - 1
    f = lambda j, k: sqrt(2./n) * E(j, n) * E(k, n) * cos(j*k*pi/n)
    return npf(f, np1)

def C_II(n):
    f = lambda j, k: sqrt(2./n) * E(j, n) * cos(j*(2*k+1)*pi/(2*n))
    return npf(f, n)

def C_III(n):
    return C_II(n).T

def C_IV(n):
    f = lambda j, k: sqrt(2./n) * cos((2*j+1)*(2*k+1)*pi/(4*n))
    return npf(f, n)

def S_I(nm1):
    n = nm1 + 1
    f = lambda j, k: sqrt(2./n) * sin((j+1)*(k+1)*pi/n)
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

def add_m(n, b):
    '''Addition matrix, defined as A in the paper'''
    n1 = n / 2
    if b == 0:
        return I(n)
    elif b == 1:
        isq2_n1m1 = I(n1-1) * 1./sqrt(2)
        m = 1./sqrt(2) * quad((I(n1-1),  I(n1-1)),
                              (I(n1-1), -I(n1-1)))
        op1 = diag(diag(1, m), -1)
        op2 = diag(I(n1), D(n1).dot(J(n1)))
        return op1.dot(op2)
    assert False

def add_m2(n, b):
    '''Second addition matrix, defined as A~ in the paper'''
    n1 = n / 2
    if b == 1:
        return I(n + 1)
    elif b == 0:
        op1 = diag(I(n1), J(n1))
        op2m = quint(I(n1-1),  J(n1-1),
                     J(n1-1), -I(n1-1),
                     sqrt(2))
        op2 = diag(1, 1./sqrt(2) * op2m)
        op3 = diag(I(n1+1), (-1)**n1 * D(n1-1))
        return op1.dot(op2).dot(op3)
    elif b == -1:
        return diag(D(n1), I(n1 - 1))
    assert False

def twiddle_m(n, b):
    '''Twiddle matrix, defined as T in the paper'''
    n1 = n / 2
    if b == 0:
        return 1./sqrt(2) * quad((I(n1),  J(n1)),
                                 (I(n1), -J(n1)))
    elif b == 1:
        m00 = diag(c(n, n1))
        m01 = diag(s(n, n1)).dot(J(n1))
        m10 = (-J(n1)).dot(diag(s(n, n1)))
        m11 = diag(J(n1).dot(c(n, n1)))
        op1 = diag(I(n1), D(n1))
        op2 = quad((m00, m01), (m10, m11))
        return op1.dot(op2)
    assert False

def twiddle_m2(n, b):
    '''Second twiddle matrix, defined as T~ in the paper'''
    n1 = n / 2
    if b == 1:
        return 1./sqrt(2) * quint(I(n1),  J(n1),
                                  I(n1), -J(n1),
                                  sqrt(2))
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
        return diag(J(n1), I(n1-1)).dot(twiddle_m2(n-1, 1))
    assert False

c = lambda n, n1: [cos((2*k+1)*pi/(4*n)) for k in range(n1)]
s = lambda n, n1: [sin((2*k+1)*pi/(4*n)) for k in range(n1)]
# XXX: page 6/37, paper has a mistake (range is set from 1 to n1-1 instead of 1
# to n1)
ct = lambda n, n1m1: [cos(k*pi/(2*n)) for k in range(1, n1m1+1)]
st = lambda n, n1m1: [sin(k*pi/(2*n)) for k in range(1, n1m1+1)]

def rectpl(x, neq2_mat, n_delay,
           twiddle_mat, add_mat, b,
           v1_func, v2_func, hvec_size_add):
    n_orig = len(x)
    n = n_orig + n_delay
    n1 = n / 2
    if n == 2:
        return neq2_mat.dot(x)
    if n >= 4:
        u = sqrt(2) * twiddle_mat(n, b).dot(x)
        v1 = v1_func(u[:n1+hvec_size_add])
        v2 = v2_func(u[n1+hvec_size_add:])
        w = add_mat(n, b).dot(np.hstack((v1, v2)))
        return permute_m(n_orig).T.dot(w)
    assert False

cosII_neq2_mat = np.array([[1,  1],
                           [1, -1]])
cosIV_neq2_mat = sqrt(2) * C_IV(2)
cosI_neq2_mat  = 1./2 * (np.array([[1,1,0],[0,0,sqrt(2)],[1,-1,0]]).dot(
                         np.array([[1,0,1],[0,sqrt(2),0],[1,0,-1]])))
cosIII_neq2_mat = 1./sqrt(2) * np.array([[1,  1],
                                         [1, -1]])
sinI_neq2_mat = I(1)

cosII  = lambda x: rectpl(x, cosII_neq2_mat,   0, twiddle_m,  add_m,   0, cosII,  cosIV,  0)
cosIV  = lambda x: rectpl(x, cosIV_neq2_mat,   0, twiddle_m,  add_m,   1, cosII,  cosII,  0)
cosI   = lambda x: rectpl(x, cosI_neq2_mat,   -1, twiddle_m2, add_m2,  1, cosI,   cosIII, 1)
cosIII = lambda x: rectpl(x, cosIII_neq2_mat,  0, twiddle_m2, add_m2,  0, cosI,   sinI,   1)
sinI   = lambda x: rectpl(x, sinI_neq2_mat,    1, twiddle_m2, add_m2, -1, cosIII, sinI,   0)

if __name__ == '__main__':
    t = int(sys.argv[1])
    n = 1<<t
    print '# t=%d -> n=%d' % (t, n)

    tests = [
        ('cos I',   cosI,   C_I,   1, lambda n: sqrt(n/2)),
        ('cos II',  cosII,  C_II,  0, lambda n: sqrt(n)),
        ('cos III', cosIII, C_III, 0, lambda n: sqrt(n/2)),
        ('cos IV',  cosIV,  C_IV,  0, lambda n: sqrt(n)),
        ('sin I',   sinI,   S_I,  -1, lambda n: sqrt(n/2)),
    ]

    for name, func, ref, delay, scale in tests:
        samples = [random.randint(0x00, 0xff) for x in range(n + delay)]

        print '%s(%d)...' % (name, n),
        fdct_ref = ref(n + delay).dot(samples)
        fdct_out = 1./scale(n) * func(samples)
        diffs = fdct_ref - fdct_out
        ok = True
        for diff in diffs:
            if abs(diff) > 1e-10:
                ok = False
                print
                print 'ref', fdct_ref
                print 'out', fdct_out
                print 'cmp', diffs
                print 'FAIL'
                break
        if ok:
            print 'OK'
