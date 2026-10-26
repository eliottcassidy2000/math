# Which dense graph sums can carry new arithmetic in the k-vertex Paley census (q = 3 mod 4)?  For k = 4..7: all
# graphs with minimum degree >= 3 (geng -d3), keep even edge counts (odd ones vanish under x -> -x), drop graphs with
# an automorphism reversing an odd number of arcs (sum = -sum), and flag the all-even-degree graphs, which reduce to
# k - 1 vertices by PGL(2) invariance (T_G = (q + 1 - k) T_(G - u), as for the octahedron).  For the survivors, the
# signed multiplicity m_H = sum over labelled copies in K_k (arcs i -> j for i < j) of the copy's sign is the exact
# coefficient of T_H in 2^C(k,2) TT_k (elimination never creates k-vertex graphs).  mac-mini-2026-10-25-pyth.
import itertools, subprocess, sys

def g6_edges(s):
    n = ord(s[0]) - 63
    bits = []
    for ch in s[1:]:
        v = ord(ch) - 63
        bits.extend((v >> (5 - i)) & 1 for i in range(6))
    E = []; k = 0
    for j in range(1, n):
        for i in range(j):
            if bits[k]: E.append((i, j))
            k += 1
    return n, E

def analyse(n, E):
    Es = set(E)
    perms = list(itertools.permutations(range(n)))
    auts = [s for s in perms if all(((min(s[i], s[j]), max(s[i], s[j])) in Es) for i, j in E)]
    def rev(s): return sum(1 for i, j in E if s[i] > s[j])
    odd = any(rev(s) % 2 for s in auts)
    deg = [sum(1 for e in E if v in e) for v in range(n)]
    signed = sum((-1) ** rev(s) for s in perms) // len(auts) if not odd else 0
    return odd, deg, signed, len(auts)

for k in range(4, int(sys.argv[1]) + 1):
    out = subprocess.run(['geng', '-q', '-d3', str(k)], capture_output=True, text=True).stdout.split()
    tot = len(out); even = vanish = pgl = 0; surv = []
    for s in out:
        n, E = g6_edges(s)
        if len(E) % 2: continue
        even += 1
        odd, deg, signed, na = analyse(n, E)
        if odd: vanish += 1; continue
        if all(d % 2 == 0 for d in deg): pgl += 1; kind = 'PGL(2)-reducible (all degrees even)'
        else: kind = 'genuine'
        surv.append((s, len(E), sorted(deg), na, signed, kind))
    print(f"k = {k}: {tot} graphs with min degree >= 3; {even} with even edge count; {vanish} vanish by an odd automorphism; "
          f"{len(surv)} survive ({pgl} PGL(2)-reducible, {len(surv) - pgl} genuine)")
    for s, ne, deg, na, signed, kind in surv:
        comp_note = ''
        print(f"    {s:10s} |E| = {ne:2d}  degrees {deg}  |Aut| = {na:4d}  signed multiplicity in K{k} = {signed:+d}  [{kind}]")
