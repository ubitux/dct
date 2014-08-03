import sys, random
import numpy as np
import sympy as sp
from math import cos, sin, sqrt, pi

import plonka

def dotx(y, matrix, x, code):
    calc = sp.Matrix(matrix).dot(x)
    if not (isinstance(calc, list) or isinstance(calc, sp.Matrix)):
        calc = [calc]
    for ij in zip(y, calc):
        code.append(ij)
    return y[:len(calc)]

def next_syms(level, x, symlevel=0):
    n = len(x)
    x0s = str(x[0])
    prefix = 'a' if x0s.startswith('src') else chr(ord(x0s[0]) + 1)*level
    return sp.Matrix([sp.Symbol('%s%d' % (prefix, i)) for i in range(n)])

def rectpl_expr(y, x, neq2_mat, n_delay,
                modified_matrix, b,
                v1_func, v2_func, hvec_size_add,
                level, scale, code):
    n_orig = len(x)
    n = n_orig + n_delay
    n1 = n / 2
    if n == 2:
        return dotx(y, neq2_mat, x, code)
    if n >= 4:
        u = dotx(next_syms(level, x), plonka.sqrt2 * plonka.twiddle_m(n, b, modified_matrix), x, code)
        vsyms = next_syms(level, u)
        v1 = v1_func(vsyms[:n1+hvec_size_add], u[:n1+hvec_size_add], level=level + 1, code=code)
        v2 = v2_func(vsyms[n1+hvec_size_add:], u[n1+hvec_size_add:], level=level + 1, code=code)
        v = np.hstack((v1, v2))
        w = plonka.add_m(n, b, modified_matrix).dot(v)
        scale_factor = 1./scale(n) if level == 1 else 1
        return dotx(y, scale_factor * plonka.permute_m(n_orig).T, w, code)
    assert False

cosII_expr  = lambda y, x, code, level=1: rectpl_expr(y, x, plonka.cosII_neq2_mat,   0, 0,  0, cosII_expr,  cosIV_expr,  0, level=level, scale=lambda z:sqrt(z), code=code)
cosIV_expr  = lambda y, x, code, level=1: rectpl_expr(y, x, plonka.cosIV_neq2_mat,   0, 0,  1, cosII_expr,  cosII_expr,  0, level=level, scale=lambda z:sqrt(z), code=code)
cosI_expr   = lambda y, x, code, level=1: rectpl_expr(y, x, plonka.cosI_neq2_mat,   -1, 1,  1, cosI_expr,   cosIII_expr, 1, level=level, scale=lambda z:sqrt(z/2), code=code)
cosIII_expr = lambda y, x, code, level=1: rectpl_expr(y, x, plonka.cosIII_neq2_mat,  0, 1,  0, cosI_expr,   sinI_expr,   1, level=level, scale=lambda z:sqrt(z/2), code=code)
sinI_expr   = lambda y, x, code, level=1: rectpl_expr(y, x, plonka.sinI_neq2_mat,    1, 1, -1, cosIII_expr, sinI_expr,   0, level=level, scale=lambda z:sqrt(z/2), code=code)

def get_code(n, fn):
    x = sp.Matrix([sp.Symbol('src[%*d*stridea]' % (len(str(n)), i)) for i in range(n)])
    y = sp.Matrix([sp.Symbol('dst[%*d*stridea]' % (len(str(n)), i)) for i in range(n)])
    code = []
    fn(y, x, code)
    indent = 8 * ' '
    outcode = []
    syms = []
    for (i, j) in code:

        j = str(j).replace('1.0*','')

        line = '%s = %s;' % (i, j)
        if '[' not in str(i) and i not in syms:
            line = 'float ' + line
            syms.append(i)
        line = indent + line
        outcode.append(line)
    return '\n'.join(outcode)

def write_dct_code(n):
    outsrc = open('template.c').read() #% tpldata
    outsrc = outsrc.replace('%N%', str(n))
    outsrc = outsrc.replace('%CODE_FDCT%', get_code(n, cosII_expr))
    outsrc = outsrc.replace('%CODE_IDCT%', get_code(n, cosIII_expr))
    open('dct%d.c' % n, 'w').write(outsrc)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Usage: %s <N>' % sys.argv[0]
        sys.exit(0)
    write_dct_code(1<<int(sys.argv[1]))
