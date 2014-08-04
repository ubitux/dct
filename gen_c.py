import sys, random, re
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
    return y

call_num = 0

def next_syms(x):
    n = len(x)
    x0s = str(x[0])
    global call_num
    prefix = 'x%x_' % call_num
    call_num += 1
    return sp.Matrix([sp.Symbol('%s%xx' % (prefix, i)) for i in range(n)])

def tfm_run(name, y, x, code, scale_factor=None):
    neq2_mat, n_delay, b, v1_tfm, v2_tfm, hvec_size_add, modified_matrix = plonka.tfm_props[name]
    n_orig = len(x)
    n = n_orig + n_delay
    n1 = n / 2
    if scale_factor is None:
        scale_factor = 1./sqrt(n>>modified_matrix)
    if n == 2:
        return dotx(y, scale_factor * neq2_mat, x, code)
    elif n >= 4:
        u = dotx(next_syms(x), plonka.sqrt2 * plonka.twiddle_m(n, b, modified_matrix), x, code)
        vsyms = next_syms(u)
        v1 = tfm_run(v1_tfm, vsyms[:n1+hvec_size_add], u[:n1+hvec_size_add], code, scale_factor)
        v2 = tfm_run(v2_tfm, vsyms[n1+hvec_size_add:], u[n1+hvec_size_add:], code, scale_factor)
        v = np.hstack((v1, v2))
        w = plonka.add_m(n, b, modified_matrix).dot(v)
        return dotx(y, plonka.permute_m(n_orig).T, w, code)
    assert False

def get_code(n, fn):
    global call_num
    call_num = 0
    x = sp.Matrix([sp.Symbol('src[%*d*src_stridea]' % (len(str(n)), i)) for i in range(n)])
    y = sp.Matrix([sp.Symbol('dst[%*d*dst_stridea]' % (len(str(n)), i)) for i in range(n)])
    code = []
    tfm_run(fn, y, x, code)
    indent = 8 * ' '
    outcode = []
    aliases = {}
    for (dst, src) in code:
        dst = str(dst)
        src = str(src).replace('1.0*','')

        s = src.split()
        if len(s) == 3 and s[1] in ('-', '+'):
            a = s[0].split('*')
            b = s[2].split('*')
            if len(a) == 2 and len(b) == 2:
                cst1, xval1 = a
                cst2, xval2 = b
                if cst1.startswith('x'): xval1, cst1 = a
                if cst2.startswith('x'): xval2, cst2 = b
                if cst1 == cst2:
                    src = '%s * (%s %s %s)' % (cst1, xval1, s[1], xval2)

        line = '%s = %s;' % (dst, src)
        if '[' not in dst:
            line = 'const float ' + line
        if re.match(r'^const float x[0-9a-f]+_[0-9a-f]+x = x[0-9a-f]+_[0-9a-f]+x;$', line):
            #outcode.append(indent + '//' + line)
            aliases[dst] = aliases.get(src, src)
            continue
        for var, rep in aliases.items():
            line = line.replace(var, rep)
        line = indent + line
        outcode.append(line)
    ret = '\n'.join(outcode)
    # suffixes were added so previous replaces don't break; we don't need them anymore
    return re.sub(r'(x[0-9a-f]+_[0-9a-f]+)x', r'\1', ret)

def write_dct_code(n):
    outsrc = open('template.c').read() #% tpldata
    outsrc = outsrc.replace('%N%', str(n))
    outsrc = outsrc.replace('%CODE_FDCT%', get_code(n, 'cosII'))
    outsrc = outsrc.replace('%CODE_IDCT%', get_code(n, 'cosIII'))
    open('dct%d.c' % n, 'w').write(outsrc)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Usage: %s <N>' % sys.argv[0]
        sys.exit(0)
    write_dct_code(1<<int(sys.argv[1]))
