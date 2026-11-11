#!/usr/bin/env python3
"""Scoped suffix/carry controls; no unrestricted target coverage census."""
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
HELPER = ROOT/'05-knowledge/results/chanpak_55_105_golden_20261110/h5/check.py'
spec = spec_from_file_location('inherited_fixed_core', HELPER)
corelib = module_from_spec(spec)
spec.loader.exec_module(corelib)

def need(q, message):
    if not q:
        raise RuntimeError(message)

@lru_cache(None)
def g(n):
    if n < 2:
        return 1
    q, b = divmod(n, 2)
    return g(q) if b else g(q)+q+(q&1)

def suffix_coefficients(k, r):
    a = b = c = 0
    scale, offset = 1, 0
    for i in reversed(range(k)):
        bit = (r >> i) & 1
        if not bit:
            a += scale
            c += offset
            if scale == 1:
                b += 1
            else:
                c += offset & 1
        scale *= 2
        offset = 2*offset+bit
    need(offset == r and scale == 1 << k, 'suffix address')
    return a, b, c

def digit_matrix(m, digit):
    dim = m+3
    matrix = [[0]*dim for _ in range(dim)]
    for r in range(m):
        j, bit = divmod(m*digit+r, 2)
        matrix[r][j] += 1
        if not bit:
            matrix[r][m] += m
            matrix[r][m+1] += j+(j&1)
            matrix[r][m+2] += j+1-(j&1)
    matrix[m][m] = 2
    matrix[m][m+1] = matrix[m][m+2] = digit
    matrix[m+1+digit][m+1] = matrix[m+1+digit][m+2] = 1
    need(all(x >= 0 for row in matrix for x in row), 'nonnegative digit matrix')
    return matrix

def apply(matrix, state):
    return tuple(sum(x*y for x, y in zip(row, state)) for row in matrix)

def direct_state(m, n):
    return tuple(g(m*n+r) for r in range(m))+(n, int(n%2 == 0), n&1)

def read(m, n, matrices):
    state = direct_state(m, 0)
    for bit in bin(n)[2:]:
        state = apply(matrices[int(bit)], state)
    return state

def direct_h(rows):
    """Independent predecessor DP, unlike imported endpoint-forward DP."""
    n = len(rows)
    incoming = [sum(1 << u for u in range(n) if rows[u] >> v & 1) for v in range(n)]
    @lru_cache(None)
    def ending(mask, v):
        rest = mask ^ (1 << v)
        if not rest:
            return 1
        return sum(ending(rest, u) for u in range(n) if rest & incoming[v] & (1 << u))
    return 1 if not n else sum(ending((1 << n)-1, v) for v in range(n))

def main():
    checks = Counter()
    for k in range(6):
        for r in range(1 << k):
            a, b, c = suffix_coefficients(k, r)
            for n in range(128):
                need(g((1 << k)*n+r) == g(n)+a*n+b*(n&1)+c, 'exact suffix identity')
                checks['suffix_identity'] += 1
    for q in range(256):
        need(g(4*q)-g(4*q+1)-g(4*q+2)+g(4*q+3) == -2, 'translated stencil')
        checks['translated_g_stencils'] += 1
    # Fixed cylinder before a further suffix: the coefficient remains one.
    for k in range(3):
        for r in range(1 << k):
            for ell in range(3):
                for s in range(1 << ell):
                    for q in range(16):
                        args = [((1 << k)*(4*q+i)+r)*(1 << ell)+s for i in range(4)]
                        val = g(args[0])-g(args[1])-g(args[2])+g(args[3])
                        need(val == -2, 'composed cylinder stencil')
                        checks['cylinder_compositions'] += 1
    cores = ((0,), (2,4,1), (4,1,10,3), (2,4,9,19,7))
    for core in cores:
        cf = corelib.coefficients(core)
        for cut in cf['cuts']:
            for b in range(8):
                for q in range(8):
                    f = [corelib.value(cf, cut, 4*q+i, b) for i in range(4)]
                    need(f[0]-f[1]-f[2]+f[3] == -2*cf['h'], 'core-count detector')
                    checks['core_stencils'] += 1
    carry_receipts = {}
    for m in (1,3,5,7,9,15):
        matrices = (digit_matrix(m, 0), digit_matrix(m, 1))
        need(apply(matrices[0], direct_state(m, 0)) == direct_state(m, 0), 'leading-zero invariant')
        for n in range(256):
            state = read(m, n, matrices)
            need(state == direct_state(m, n), 'all affine-address channels')
            for bit in (0,1):
                need(apply(matrices[bit], state) == direct_state(m, 2*n+bit), 'digit transition')
                checks['digit_transitions'] += 1
            checks['carry_input_addresses'] += 1
        carry_receipts[m] = {'dimension': m+3, 'sample_n': 37,
                             'sample_state': read(m, 37, matrices)}
        if m == 3:
            carry_receipts[m]['matrices'] = matrices
    # Literal graph readouts for both one-halved odd affine maps.
    core = (2,4,1)
    cf = corelib.coefficients(core)
    cut = cf['cuts'][1]
    graph_receipts = []
    for n in range(6):
        state = read(3, n, (digit_matrix(3,0), digit_matrix(3,1)))
        for r in (1,2):
            address = 3*n+r
            value = 3*state[r]+15*address+10+(address&1)
            rows = corelib.graph(core, 1, address, 1)
            need(corelib.tournament(rows), 'actual tournament')
            need(value == corelib.value(cf,cut,1,address) == direct_h(rows), 'literal address readout')
            checks['literal_graph_readouts'] += 1
            if n <= 1:
                graph_receipts.append({'n':n,'channel':r,'address':address,'H':value,'rows':rows})
    old = corelib.value(cf,cut,1,1)
    new = corelib.value(cf,cut,1,2)
    need((old,new) == (29,49) and new not in (3*old+2,3*old+4), 'address/count separation hostile')
    dependencies = [HELPER,
        ROOT/'05-knowledge/results/chanpak_55_105_golden_20261110/h5/PROOF.md',
        ROOT/'05-knowledge/results/session253_cross_problem_frontiers_20261111/integral/PROOF.md',
        ROOT/'05-knowledge/results/h13_paired_descent_20261111/square_graceful/RESIDUE_PHASE.md']
    print(json.dumps({'status':'PASS','scope':'fixed suffixes<=5 bits; specified small core cuts; six odd multipliers; literal graph controls',
                      'checks':dict(checks),'carry_lifts':carry_receipts,'graph_controls':graph_receipts,
                      'hostile':{'old_address':1,'new_address':2,'old_H':old,'new_H':new},
                      'source_sha256':{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in dependencies},
                      'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest()},indent=2,sort_keys=True))

if __name__ == '__main__':
    main()
