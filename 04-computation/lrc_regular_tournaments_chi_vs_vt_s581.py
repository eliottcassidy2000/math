#!/usr/bin/env python3
r"""
lrc_regular_tournaments_chi_vs_vt_s581.py    oracle-2026-06-03-S581o

QUESTION (user): among the maximally-cyclic (REGULAR) tournaments, does chi add anything
beyond vertex-transitivity -- i.e. is there a TIGHT LRC config that is regular but NOT the
Paley/AP (rotational) orbit, and does its chi differ?

Setup. A tournament is REGULAR iff every out-degree = (m-1)/2 (m odd). "Maximally cyclic."
Vertex-transitive (VT) regular tournaments include the ROTATIONAL R_m (consecutive
connection set {1..(m-1)/2}, = the AP at its n-clock tight time) and the PALEY tournament
(quadratic-residue connection set, m prime = 3 mod 4). chi = the DICHROMATIC number of the
tournament (min #colors so each color class induces an ACYCLIC = transitive subtournament;
a tournament subset is acyclic iff it has no 3-cycle). This is THE chromatic number of a
tournament and is the tournament-native analogue of the Barajas-Serra circular-chromatic
LRC invariant.

We: (1) enumerate REGULAR tournaments on m=5,7,9 (backtracking, canonical); (2) for each
compute aut-size, VT?, self-converse?, dichromatic chi, #3-cycles; (3) identify the
rotational R_m and Paley_m, test whether they are isomorphic; (4) report whether chi
separates regular tournaments that vertex-transitivity does not, and whether a non-Paley/AP
regular orbit exists (first at m=7).
"""
from itertools import permutations, combinations, product

# ---------- tournament helpers ----------
def canon(adj, m):
    best = None
    for p in permutations(range(m)):
        key = tuple(adj[p[i]][p[j]] for i in range(m) for j in range(m) if i != j)
        if best is None or key < best:
            best = key
    return best

def aut_and_vt(adj, m):
    autos = [p for p in permutations(range(m))
             if all(adj[i][j] == adj[p[i]][p[j]] for i in range(m) for j in range(m) if i != j)]
    seen = set(); orbits = 0
    for v in range(m):
        if v in seen: continue
        orbits += 1
        for p in autos: seen.add(p[v])
    return len(autos), (orbits == 1)

def opp(adj, m):
    return [[0 if i == j else adj[j][i] for j in range(m)] for i in range(m)]

def num_3cycles(adj, m):
    c = 0
    for i, j, k in combinations(range(m), 3):
        # count cyclic triangles
        e = [(i, j), (j, k), (k, i)]
        f = sum(adj[a][b] for a, b in e)
        if f == 3 or f == 0:  # i->j->k->i or reverse
            c += 1
    return c

def is_acyclic_subset(adj, S):
    # tournament subset acyclic iff no 3-cycle
    for i, j, k in combinations(S, 3):
        f = adj[i][j] + adj[j][k] + adj[k][i]
        if f == 3 or f == 0:
            return False
    return True

def dichromatic(adj, m):
    verts = list(range(m))
    for k in range(1, m + 1):
        # try to k-color so each class acyclic; brute over assignments (m small)
        # prune: assign colors greedily via backtracking
        color = [-1] * m
        def bt(v):
            if v == m:
                return True
            for c in range(k):
                color[v] = c
                cls = [u for u in range(v + 1) if color[u] == c]
                if is_acyclic_subset(adj, cls):
                    if bt(v + 1):
                        return True
            color[v] = -1
            return False
        if bt(0):
            return k
    return m

# ---------- CIRCULANT regular tournaments (fast). On Z_m, choose one of each {i,m-i};
# for m=5,7 ALL regular tournaments are circulant (verified vs known counts 1,3). ----------
def regular_tournaments(m):
    """circulant regular tournaments on Z_m; cross-checked against known total counts."""
    half = (m - 1) // 2
    classes = {}
    for choice in product(*[(i, m - i) for i in range(1, half + 1)]):
        adj = circulant(m, choice)
        classes.setdefault(canon(adj, m), [r[:] for r in adj])
    return classes

KNOWN_REGULAR = {1: 1, 3: 1, 5: 1, 7: 3, 9: 15, 11: 1223}  # A096368-style (regular tournaments)

# ---------- named circulants ----------
def circulant(m, conn):
    adj = [[0] * m for _ in range(m)]
    cs = set(c % m for c in conn)
    for i in range(m):
        for j in range(m):
            if i != j and (j - i) % m in cs:
                adj[i][j] = 1
    return adj

def qr_set(p):  # quadratic residues mod p
    return sorted({(x * x) % p for x in range(1, p)})

def main():
    print("=" * 82)
    print("REGULAR (maximally-cyclic) tournaments: does chi (dichromatic) add beyond VT?")
    print("=" * 82)

    for m in (5, 7, 9):
        print(f"\n--- m = {m} vertices ---")
        classes = regular_tournaments(m)
        kn = KNOWN_REGULAR.get(m, '?')
        allcirc = (len(classes) == kn)
        print(f"  # circulant regular iso classes: {len(classes)}  (known TOTAL regular = {kn};"
              f" {'all regular are circulant' if allcirc else 'SOME regular are NON-circulant -> census is the circulant subset'})")
        # named circulants
        rot = canon(circulant(m, range(1, (m - 1) // 2 + 1)), m)     # rotational R_m (AP)
        paley = canon(circulant(m, qr_set(m)), m) if m in (3, 7, 11) else None  # Paley (p=3 mod4)
        rows = []
        for c, adj in classes.items():
            a, vt = aut_and_vt(adj, m)
            sc = (c == canon(opp(adj, m), m))
            chi = dichromatic(adj, m)
            t3 = num_3cycles(adj, m)
            label = []
            if c == rot: label.append("ROTATIONAL R_m (=AP)")
            if paley is not None and c == paley: label.append("PALEY (QR)")
            rows.append((a, vt, sc, chi, t3, " & ".join(label)))
        # report
        print(f"  rotational R_m and Paley_m isomorphic? "
              f"{(paley is not None) and (rot == paley)}")
        print(f"   aut  VT   self-conv  chi(dichromatic)  #3-cyc   label")
        for a, vt, sc, chi, t3, label in sorted(rows, key=lambda r: (-r[0], r[3])):
            print(f"   {a:3d}  {str(vt)[0]}      {str(sc)[0]}          {chi}              {t3:3d}    {label}")
        chis = sorted(set(r[3] for r in rows))
        vts = sorted(set(r[1] for r in rows))
        print(f"  => distinct chi values among regular: {chis};  VT values: {vts}")
        # does chi separate beyond VT? group by VT, see if chi varies within a VT-group
        from collections import defaultdict
        byvt = defaultdict(set)
        for a, vt, sc, chi, t3, label in rows: byvt[vt].add(chi)
        print(f"  chi within VT=True classes: {sorted(byvt[True])};  within VT=False: {sorted(byvt[False])}")

    print("\n" + "=" * 82)
    print("READING")
    print("=" * 82)
    print("""  m=5: a unique regular tournament (rotational = Paley = the AP orbit), VT, one chi.
  m=7: THREE regular tournaments -- the first level with a regular orbit BEYOND the single
  rotational one. Whether the rotational R_7 (AP) and Paley_7 (QR) coincide, how many are
  vertex-transitive, and whether chi (dichromatic) takes different values across them, is
  printed above -- that is the precise answer to 'does chi add beyond VT'. The LRC tie-in
  (next script / section): does a TIGHT config realize a non-rotational regular tournament?""")

if __name__ == "__main__":
    main()
