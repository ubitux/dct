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

I = lambda n: np.identity(n, dtype=int)
J = lambda n: np.fliplr(I(n))

T_0_no_scale = lambda n: quad((I(n / 2),  J(n / 2)),
                              (I(n / 2), -J(n / 2)))

diag = lambda a, b: quad((a, np.zeros((a.shape[0], b.shape[1]))),
                         (np.zeros((b.shape[0], a.shape[1])), b))

def T_1(n):
    i_n1 = I(n / 2)
    j_n1 = J(n / 2)
    cos_vec = [cos((2*k+1)*pi/(4*n)) for k in range(n/2)]
    sin_vec = [sin((2*k+1)*pi/(4*n)) for k in range(n/2)][::-1]
    rev_cos_sign_vec = [(-x if i & 1 else x) for i, x in enumerate(cos_vec[::-1])]
    rev_sin_sign_vec = [(-x if i & 1 else x) for i, x in enumerate(sin_vec[::-1])]
    m00 = np.array(cos_vec) * i_n1
    m01 = np.array(sin_vec) * j_n1
    m10 = np.array(rev_sin_sign_vec) * j_n1
    m11 = np.array(rev_cos_sign_vec) * i_n1
    return quad((m00, m01),
                (m10, m11))

def D(n):
    m = np.zeros((n, n), dtype=int)
    for i in range(n):
        m[i, i] = -1 if i & 1 else 1
    return m

def A_1(n):
    isq2_n1m1 = I(n/2-1) * 1./sqrt(2)
    m = quad((isq2_n1m1,  isq2_n1m1),
             (isq2_n1m1, -isq2_n1m1))
    op1 = diag(diag(np.array([[1]]), m), np.array([[-1]]))
    op2 = diag(I(n / 2), np.dot(D(n / 2), J(n / 2)))
    return np.dot(op1, op2)

def permute(v):
    n = len(v)
    while n > 2:
        v = np.hstack((v[::2], v[1::2]))
        n /= 2
    return v

def cosII(x):
    n = len(x)
    #print 'cosII n=%d' % n
    if n == 2:
        return np.dot(np.array([[1,  1],
                                [1, -1]]), x)
    if n >= 4:
        u = np.dot(T_0_no_scale(n), x)
        v1 = cosII(u[:n/2])
        v2 = cosIV(u[n/2:])
        w = np.hstack((v1, v2))
        return permute(w)
    assert False

def cosIV(x):
    n = len(x)
    #print 'cosIV n=%d' % n
    if n == 2:
        return sqrt(2) * np.dot(C_IV(2), x)
    if n >= 4:
        u = np.dot(sqrt(2) * T_1(n), x)
        v1 = cosII(u[:n/2])
        v2 = cosII(u[n/2:])
        w = np.dot(A_1(n), np.hstack((v1, v2)))
        return permute(w)
    assert False

def C(f, n):
    m = []
    for j in range(n):
        m.append([f(j, k) for k in range(n)])
    return np.array(m)

def C_II(n):
    f = lambda j, k: sqrt(2./n) * (1./sqrt(2) if j == 0 else 1) * cos(j*(2*k+1)*pi/(2*n))
    return C(f, n)

def C_IV(n):
    f = lambda j, k: sqrt(2./n) * cos((2*j+1)*(2*k+1)*pi/(4*n))
    return C(f, n)

t = int(sys.argv[1])
n = 1<<t
print '# t=%d -> n=%d' % (t, n)

samples = random.sample(range(0x100), n)

tests = [
    ('cos II', cosII, C_II),
    ('cos IV', cosIV, C_IV),
]

for name, func, ref in tests:
    print
    print '=== %s ===' % name
    fdct_ref = np.dot(ref(n), samples)
    fdct_out = 1./sqrt(n) * func(samples)
    print 'ref', fdct_ref
    print 'out', fdct_out
    print 'cmp', fdct_ref - fdct_out
