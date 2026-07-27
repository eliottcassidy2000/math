"""Exact structure of ker D4 on the conic-cap strata (kind-pasteur)."""
import sympy as sp
from g1_conic_cap_hunt_kps import (x, y, A_C1, A_C2, kernelB, coords_of_B,
                                   in_span, monoms_upto, dvec, log)

for tag, A in (('C1', A_C1), ('C2', A_C2)):
    Ax, Ay = dvec(A, x), dvec(A, y)
    degB = 3
    basis, bs, M = kernelB(A, degB)
    # module span: h*A (deg h <=1), f*A_x (deg f <=2), g*A_y (deg g <=2), + all linear B
    gens = []
    names = []
    for h in [x**i*y**j for (i, j) in monoms_upto(1)]:
        gens.append([sp.expand(h*a) for a in A]); names.append('A*%s' % h)
    for f in [x**i*y**j for (i, j) in monoms_upto(2)]:
        gens.append([sp.expand(f*a) for a in Ax]); names.append('Ax*%s' % f)
        gens.append([sp.expand(f*a) for a in Ay]); names.append('Ay*%s' % f)
    lin = []
    for comp in range(3):
        for m in (x, y):
            v = [sp.Integer(0)]*3
            v[comp] = m
            lin.append(v)
    genc = [coords_of_B(g, degB) for g in gens]
    linc = [coords_of_B(v, degB) for v in lin]
    kerc = [coords_of_B(b, degB) for b in basis]
    # 1) is every kernel element in module+linear?
    allgen = genc + linc
    ok1 = all(in_span(allgen, k) for k in kerc)
    # 2) is every module element + linear in the kernel?
    ok2 = all(in_span(kerc, g) for g in allgen)
    # 3) module alone (no linear): which kernel elements escape?
    only_mod = sum(1 for k in kerc if in_span(genc, k))
    log('%s: ker D4 == tangent-module{hA+fA_x+gA_y} + linearB : %s (both inclusions %s/%s);'
        ' kernel elements inside pure module: %d/%d' % (tag, ok1 and ok2, ok1, ok2, only_mod, len(kerc)))
    # constants never in kernel?
    const = []
    for comp in range(3):
        v = [sp.Integer(0)]*3
        v[comp] = sp.Integer(1)
        const.append(coords_of_B(v, degB))
    ok3 = [in_span(kerc, c) for c in const]
    log('%s: constant B_e1,e2,e3 in kernel? %s' % (tag, ok3))
log('DONE')
