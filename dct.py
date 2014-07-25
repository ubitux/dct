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

def diag(a, b=None):
    if b is None:
        return a * I(len(a))
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

def P(n):
    #assert n >= 4
    m = np.zeros((n, n))
    col = 0
    for i in range(n):
        m[i, col] = 1
        col += 2
        if col >= n:
            col = 1
    return m.T

def A_1(n):
    n1 = n / 2
    isq2_n1m1 = I(n1-1) * 1./sqrt(2)
    m = 1./sqrt(2) * quad((I(n1-1),  I(n1-1)),
                          (I(n1-1), -I(n1-1)))
    op1 = diag(diag(np.array([[1]]), m), np.array([[-1]]))
    op2 = diag(I(n1), np.dot(D(n1), J(n1)))
    return np.dot(op1, op2)

def At_0(n):
    n1 = n / 2
    op1 = diag(I(n1), J(n1))
    op2m = quint(I(n1-1),  J(n1-1),
                 J(n1-1), -I(n1-1),
                 sqrt(2))
    op2 = diag(np.array([[1]]), 1./sqrt(2) * op2m)
    op3 = diag(I(n1+1), (-1)**n1 * D(n1-1))
    return np.dot(op1, op2, op3)

def At_m1(nm1):
    n = nm1 + 1
    n1 = n / 2
    return diag(D(n1), I(n1 - 1))

def T_0(n):
    n1 = n / 2
    return 1./sqrt(2) * quad((I(n1),  J(n1)),
                             (I(n1), -J(n1)))

def Tt_1(np1):
    n = np1 - 1
    n1 = n / 2
    return 1./sqrt(2) * quint(I(n1),  J(n1),
                              I(n1), -J(n1),
                              sqrt(2))

def Tt_m1(nm1):
    n = nm1 + 1
    n1 = n / 2
    return np.dot(diag(J(n1), I(n1-1)), Tt_1(nm1))

def T_1(n):
    n1 = n / 2
    m00 = diag(c(n, n1))
    m01 = np.dot(diag(s(n, n1)), J(n1))
    m10 = np.dot(-J(n1), diag(s(n, n1)))
    m11 = diag(np.dot(J(n1), c(n, n1)))
    op1 = diag(I(n1), D(n1))
    op2 = quad((m00, m01), (m10, m11))
    return np.dot(op1, op2)

def Tt_0(n):
    n1 = n / 2
    op1 = diag(I(n1+1), D(n1-1))
    m00 = diag(ct(n, n1-1))
    m01 = np.dot(diag(s(n, n1-1)), J(n1-1))
    m10 = np.dot(-J(n1-1), diag(s(n, n1-1)))
    m11 = diag(np.dot(J(n1-1), ct(n, n1-1)))
    m = quint(m00, m01, m10, m11, 1)
    op2 = diag(np.array([[1]]), m)
    return np.dot(op1, op2)

c = lambda n, n1: [cos((2*k+1)*pi/(4*n)) for k in range(n1)]
s = lambda n, n1: [sin((2*k+1)*pi/(4*n)) for k in range(n1)]
ct = lambda n, n1m1: [cos(k*pi/(2*n)) for k in range(n1m1)]
st = lambda n, n1m1: [sin(k*pi/(2*n)) for k in range(n1m1)]

def cosII(x, indent=0):
    n = len(x)
    print (' '*indent*4) + 'cosII n=%d (x=%s)' % (n, x)
    n1 = n / 2
    if n == 2:
        return np.dot(np.array([[1,  1],
                                [1, -1]]), x)
    if n >= 4:
        u = np.dot(sqrt(2) * T_0(n), x)
        v1 = cosII(u[:n1], indent + 1)
        v2 = cosIV(u[n1:], indent + 1)
        return np.dot(P(n), np.hstack((v1, v2)))
    assert False

def cosIV(x, indent=0):
    n = len(x)
    print (' '*indent*4) + 'cosIV n=%d (x=%s)' % (n, x)
    n1 = n / 2
    if n == 2:
        return sqrt(2) * np.dot(C_IV(2), x)
    if n >= 4:
        u = np.dot(sqrt(2) * T_1(n), x)
        v1 = cosII(u[:n1], indent + 1)
        v2 = cosII(u[n1:], indent + 1)
        w = np.dot(A_1(n), np.hstack((v1, v2)))
        return np.dot(P(n), w)
    assert False

def cosI(x, indent=0):
    np1 = len(x)
    print (' '*indent*4) + 'cosI n+1=%d (x=%s)' % (np1, x)
    n = np1 - 1
    n1 = n / 2
    if n == 2:
        return 1./2 * np.dot(np.dot(np.array([[1,1,0],[0,0,sqrt(2)],[1,-1,0]]),
                                    np.array([[1,0,1],[0,sqrt(2),0],[1,0,-1]])), x)
    if n >= 4:
        u = np.dot(sqrt(2) * Tt_1(n+1), x)
        v1 = cosI(u[:n1+1], indent + 1)
        v2 = cosIII(u[n1+1:], indent + 1)
        return np.dot(P(n + 1), np.hstack((v1, v2)))
    assert False

def cosIII(x, indent=0):
    n = len(x)
    print (' '*indent*4) + 'cosIII n=%d (x=%s)' % (n, x)
    n1 = n / 2
    if n == 2:
        return 1./sqrt(2) * np.dot(np.array([[1,  1],
                                             [1, -1]]), x)
    if n >= 4:
        u = np.dot(sqrt(2) * Tt_0(n), x)
        v1 = cosI(u[:n1+1], indent + 1)
        v2 = sinI(u[n1+1:], indent + 1)
        w = np.dot(At_0(n), np.hstack((v1, v2)))
        return np.dot(P(n), w)

def sinI(x, indent=0):
    nm1 = len(x)
    print (' '*indent*4) + 'sinI n-1=%d (x=%s)' % (nm1, x)
    n = nm1 + 1
    n1 = n / 2
    if n == 2:
        return np.array(x)
    if n >= 4:
        u = np.dot(sqrt(2) * Tt_m1(nm1), x)
        print u
        assert len(u) == n-1
        v1 = cosIII(u[:n1], indent + 1)
        v2 = sinI(u[n1:], indent + 1)
        w = np.dot(At_m1(n - 1), np.hstack((v1, v2)))
        return np.dot(P(n - 1), w)
    assert False

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
    samples = random.sample(range(0x100), n + delay)

    print
    print '=== %s ===' % name
    fdct_ref = np.dot(ref(n + delay), samples)
    fdct_out = 1./scale(n) * func(samples)
    print 'ref', fdct_ref
    print 'out', fdct_out
    diffs = fdct_ref - fdct_out
    ok = True
    for diff in diffs:
        if abs(diff) > 1e-12:
            ok = False
            print 'cmp', diffs
            print 'FAIL'
            break
    if ok:
        print 'OK'
