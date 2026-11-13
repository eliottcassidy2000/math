"""Exact independent controls for the uniform-module root-phase theorem.

Run python -B root_phase.py or python -O -B root_phase.py. No asserts.
The universe is all labelled contexts with 0..3 exterior vertices, every
module size 1..6, and every root position; plus 16 deterministic 4-exterior
contexts with sizes 2..5. Literal permutations independently check n<=6.
"""
from itertools import combinations, permutations, product
from math import comb, gcd
import json
from pathlib import Path
import sys


def require(test, message):
    if not test:
        raise RuntimeError(message)


def tournament(e, code, m, phase):
    """Bits: exterior pairs first, then module beats each exterior vertex."""
    n = e + m
    out = [0] * n
    bit = 0
    def edge(a, b):
        out[a] |= 1 << b
    for a, b in combinations(range(e), 2):
        edge(m+a, m+b) if code >> bit & 1 else edge(m+b, m+a)
        bit += 1
    for a in range(e):
        for v in range(m):
            edge(v, m+a) if code >> bit & 1 else edge(m+a, v)
        bit += 1
    order = list(range(1, phase+1)) + [0] + list(range(phase+1, m))
    for a, b in combinations(order, 2):
        edge(a, b)
    for a, b in combinations(range(n), 2):
        require(bool(out[a] >> b & 1) != bool(out[b] >> a & 1), 'tournament')
    return out


def state(out):
    n = len(out)
    total = 1 << n
    dp = [[0]*n for _ in range(total)]
    h = [0]*total
    h[0] = 1
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, total):
        for last, count in enumerate(dp[mask]):
            if not count:
                continue
            nxt = out[last] & (total-1-mask)
            while nxt:
                one = nxt & -nxt
                dp[mask | one][one.bit_length()-1] += count
                nxt -= one
        h[mask] = sum(dp[mask])
    full = total-1
    # A separate start-at-root DP avoids dependence on a reversal convention.
    rooted = [[0]*n for _ in range(total)]
    rooted[1][0] = 1
    for mask in range(1, total, 2):
        for last, count in enumerate(rooted[mask]):
            if not count:
                continue
            nxt = out[last] & (full-mask)
            while nxt:
                one = nxt & -nxt
                rooted[mask | one][one.bit_length()-1] += count
                nxt -= one
    F = sum(dp[mask][0] * h[full ^ mask] for mask in range(1, total, 2))
    return (h[full], sum(rooted[full]), h[full ^ 1], F)


def literal_state(out):
    n = len(out)
    def paths(vertices):
        return [p for p in permutations(vertices)
                if all(out[a] >> b & 1 for a, b in zip(p, p[1:]))]
    full = (1 << n)-1
    hp = paths(range(n))
    F = 0
    for mask in range(1, full+1, 2):
        first = [v for v in range(n) if mask >> v & 1]
        second = [v for v in range(n) if not (mask >> v & 1)]
        F += sum(p[-1] == 0 for p in paths(first)) * len(paths(second))
    return (len(hp), sum(p[0] == 0 for p in hp), len(paths(range(1,n))), F)


def S(j, l, r):
    return sum((-1)**(j-1-i)*comb(j-1,i)*i**l*(i+1)**r for i in range(j))


def kernel_assignments(j, l, r):
    # Every assignment gives ordered increasing paths; bin 0 contains root.
    return sum(set(a) | {0} == set(range(j)) and all(x != 0 for x in a[:l])
               for a in product(range(j), repeat=l+r))


def differences(values):
    return [b-a for a,b in zip(values, values[1:])]


def controller(st, k):
    h,s,q,F = st
    d = F-q
    D = d-(h-s)
    t = 2**k*d+q-D
    b = 3*t+2*D
    return b,t


def native(b,t):
    return b > t > 0 and gcd(b,t) == 1 and t % 3 == 0 and b % 3 != 0


def main():
    counts = dict(contexts=0, tournaments=0, literal_checks=0, adjacent_steps=0,
                  strict_even_differences=0, kernel_checks=0,
                  contraction_differences=0, odd_source_equalities=0,
                  binomial_reconstructions=0)
    for m in range(1,7):
        for l in range(m):
            r = m-1-l
            for j in range(1,m+1):
                require(S(j,l,r) == kernel_assignments(j,l,r), 'kernel assignments')
                if r:
                    require(S(j,l+1,r-1)-S(j,l,r) == -S(j,l,r-1), 'S derivative')
                    require(S(j,r-1,l+1)-S(j,r,l) == S(j,r-1,l), 'E derivative')
                counts['kernel_checks'] += 1
    for e in range(5):
        codes = range(1 << (e*(e-1)//2 + e)) if e <= 3 else [i*67 for i in range(16)]
        sizes = range(1,7) if e <= 3 else range(2,6)
        for code in codes:
            for m in sizes:
                sts, graphs = [], []
                for l in range(m):
                    graph = tournament(e, code, m, l)
                    st = state(graph)
                    if e+m <= 6:
                        require(st == literal_state(graph), 'literal independent state')
                        counts['literal_checks'] += 1
                    graphs.append(graph)
                    sts.append(st)
                    counts['tournaments'] += 1
                require(len({s[0] for s in sts}) == 1, 'fixed h')
                require(len({s[2] for s in sts}) == 1, 'fixed q')
                require(all(x <= 0 for x in differences([s[1] for s in sts])), 's monotone')
                require(all(x > 0 for x in differences([s[3] for s in sts])), 'F strict')
                Ds = [F+s-h-q for h,s,q,F in sts]
                first_states = {a:state(tournament(e,code,a,0)) for a in range(1,m+1)}
                for l in range(m):
                    binomial_s = sum((-1)**t*comb(l,t)*first_states[m-t][1] for t in range(l+1))
                    binomial_F = sum(comb(l,t)*first_states[m-t][3] for t in range(l+1))
                    require((binomial_s,binomial_F) == (sts[l][1],sts[l][3]),
                            'binomial root response reconstruction')
                    counts['binomial_reconstructions'] += 1
                for l in range(m-1):
                    changed = sum((graphs[l][v] ^ graphs[l+1][v]).bit_count() for v in range(e+m))
                    require(changed == 2, 'one arc step')
                    counts['adjacent_steps'] += 1
                for order in range(1,m):
                    Ds = differences(Ds)
                    for l,value in enumerate(Ds):
                        contracted = state(tournament(e,code,m-order,l))
                        require(value == contracted[3]+(-1)**order*contracted[1],
                                'all-order contracted response')
                        source = all(not (row & 1) for row in tournament(e,code,m-order,l))
                        require(value >= 0 and ((value == 0) == (order % 2 == 1 and source)),
                                'all-order exact equality boundary')
                        counts['contraction_differences'] += 1
                        counts['odd_source_equalities'] += value == 0
                    if order % 2 == 0:
                        require(all(x > 0 for x in Ds), 'strict even difference')
                        counts['strict_even_differences'] += len(Ds)
                counts['contexts'] += 1
    # Minimal hostile: exterior vertex distinguishes the two module vertices.
    nonuniform_before = [0b010, 0b100, 0b001]  # directed cycle
    nonuniform_after = [0, 0b101, 0b001]      # reverse only arc 0->1
    require(sum((a^b).bit_count() for a,b in zip(nonuniform_before,nonuniform_after)) == 2,
            'single-arc hostile')
    require([state(g)[0] for g in [nonuniform_before,nonuniform_after]] == [3,1],
            'uniformity essential')
    # Homogeneous-pair flips do not exhaust all H-neutral arc flips.
    other_before, other_after = [8,1,3,6], [12,1,2,6]
    require(sum((a^b).bit_count() for a,b in zip(other_before,other_after)) == 2,
            'non-twin neutral is one arc')
    require(((other_before[0]^other_before[2]) & ~5) != 0,
            'non-twin neutral has distinguishable endpoints')
    require(state(other_before)[0] == state(other_after)[0] == 5,
            'non-twin neutral H')
    require(state(other_before) == literal_state(other_before) and
            state(other_after) == literal_state(other_after),
            'non-twin neutral independent states')
    # C3[TT3,1,1] exact arithmetic hostile: exterior 3->4, C->3, 4->C.
    book = [state(tournament(2,3,3,l)) for l in range(3)]
    require(book == [(9,4,5,10),(9,2,5,16),(9,1,5,25)], 'H9 states')
    b,t = controller(book[1],2)
    require((b,t) == (143,45) and native(b,t), 'interior native entry')
    require(all(not native(*controller(book[l],k)) for l in (0,2) for k in range(1,41)),
            'endpoint hostile finite control')
    report = {'status':'FINITE-EXACT', 'python_optimized':not __debug__, 'counts':counts,
              'H9_root_states':book, 'H9_D':[F+s-h-q for h,s,q,F in book],
              'H9_interior_R2NR':{'b':b,'t':t},
              'nonuniform_single_arc_H':[3,1],
              'non_twin_neutral_arc':{'before':other_before,'after':other_after,'H':5},
              'scope':'All labelled homogeneous exterior contexts e<=3, m<=6; 16 deterministic e=4 contexts m=2..5. General theorem is in PROOF.md.'}
    name = 'verification.optimized.json' if not __debug__ else 'verification.json'
    Path(__file__).with_name(name).write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
