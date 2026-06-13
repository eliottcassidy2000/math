#!/usr/bin/env python3
r"""
lrc_round_regular_unique_Rm_s582.py    oracle-2026-06-03-S582o (overnight cycle 1)

Completes the S581o handoff + opus-S591's open qualifiers. opus-S591: LRC's half-turn
comparator => out-neighbourhoods are contiguous arcs => the tournament is ROUND (locally
transitive); LRC-accessible = round only; the interval circulant R_m is round, Paley is not
(m>=7) => Paley inaccessible. OPEN: (a) is R_m the UNIQUE round regular tournament? (b) is
R_m the unique chi=2 regular tournament? -- i.e. does 'round + regular' (or 'chi=2 + regular')
pin the AP orbit uniquely?

We enumerate ALL regular tournaments on m=7,9 (backtracking, deduped by skew-eigenvalue
spectrum), and for each compute: round? (locally transitive), chi (dichromatic), VT, aut.
"""
import numpy as np
from itertools import combinations, permutations

def out_degs(adj, m): return [sum(adj[i]) for i in range(m)]

def locally_transitive(adj, m):
    def trans(S):
        S = list(S)
        return sorted(sum(1 for b in S if a != b and adj[a][b]) for a in S) == list(range(len(S)))
    for v in range(m):
        outn = [w for w in range(m) if w != v and adj[v][w]]
        inn = [w for w in range(m) if w != v and adj[w][v]]
        if not trans(outn) or not trans(inn):
            return False
    return True

def is_acyclic(adj, S):
    for i, j, k in combinations(S, 3):
        f = adj[i][j] + adj[j][k] + adj[k][i]
        if f in (0, 3): return False
    return True

def dichromatic(adj, m):
    for k in range(1, m + 1):
        color = [-1] * m
        def bt(v):
            if v == m: return True
            for c in range(k):
                color[v] = c
                if is_acyclic(adj, [u for u in range(v + 1) if color[u] == c]) and bt(v + 1):
                    return True
            color[v] = -1; return False
        if bt(0): return k
    return m

def spectrum_key(adj, m):
    # skew-adjacency S = adj - adj^T (entries +-1), eigenvalues are purely imaginary;
    # use sorted |Im(eig)| rounded as an iso-invariant
    A = np.array(adj) - np.array(adj).T
    ev = np.linalg.eigvals(A)
    return tuple(sorted(np.round(np.abs(ev), 4)))

def canon(adj, m):
    best = None
    for p in permutations(range(m)):
        key = tuple(adj[p[i]][p[j]] for i in range(m) for j in range(m) if i != j)
        if best is None or key < best: best = key
    return best

def aut_vt(adj, m):
    autos = [p for p in permutations(range(m))
             if all(adj[i][j] == adj[p[i]][p[j]] for i in range(m) for j in range(m) if i != j)]
    seen = set(); orb = 0
    for v in range(m):
        if v in seen: continue
        orb += 1
        for p in autos: seen.add(p[v])
    return len(autos), orb == 1

def enumerate_regular(m):
    """all regular tournaments on m, deduped by (spectrum, canon) -- canon only within spectrum bucket."""
    d = (m - 1) // 2
    pairs = list(combinations(range(m), 2))
    adj = [[0] * m for _ in range(m)]
    out = [0] * m
    reps = {}      # spectrum -> adj (one representative; verify count vs known)
    npairs = len(pairs)
    def feasible(idx):
        # quick: every vertex still able to reach exactly d
        rem_out = [0] * m
        for (i, j) in pairs[idx:]:
            rem_out[i] += 1; rem_out[j] += 1
        for v in range(m):
            if out[v] > d or out[v] + rem_out[v] < d:
                return False
        return True
    def bt(idx):
        if idx == npairs:
            if all(out[v] == d for v in range(m)):
                sk = spectrum_key(adj, m)
                if sk not in reps:
                    reps[sk] = [r[:] for r in adj]
            return
        i, j = pairs[idx]
        for (a, b) in ((i, j), (j, i)):
            if out[a] < d:
                adj[a][b] = 1; out[a] += 1
                if feasible(idx + 1):
                    bt(idx + 1)
                adj[a][b] = 0; out[a] -= 1
    bt(0)
    return list(reps.values())

def circulant(m, conn):
    cs = set(c % m for c in conn)
    return [[1 if i != j and (j - i) % m in cs else 0 for j in range(m)] for i in range(m)]

def main():
    print("=" * 78)
    print("Is R_m the UNIQUE round (and unique chi=2) regular tournament? (S581o/opus-S591)")
    print("=" * 78)
    known = {7: 3, 9: 15}
    for m in (7, 9):
        cls = enumerate_regular(m)
        Rm_adj = circulant(m, range(1, (m - 1) // 2 + 1))
        Rm_spec = spectrum_key(Rm_adj, m)
        nclass = len(cls)
        print(f"\n--- m={m}: {nclass} regular tournaments by SPECTRUM (known total {known[m]};"
              f" {'spectrum-complete' if nclass == known[m] else 'COSPECTRAL merging -> lower bound'}) ---")
        heavy = (m <= 7)
        print("   round?  chi  " + ("VT  aut  " if heavy else "") + "is_R_m(interval/AP)?")
        n_round = 0; n_chi2 = 0; rows = []
        for adj in cls:
            rnd = locally_transitive(adj, m)
            chi = dichromatic(adj, m)
            isR = (spectrum_key(adj, m) == Rm_spec)
            if heavy:
                a, vt = aut_vt(adj, m)
            else:
                a, vt = None, None
            n_round += rnd; n_chi2 += (chi == 2)
            rows.append((rnd, chi, vt, a, isR))
        for rnd, chi, vt, a, isR in sorted(rows, key=lambda r: (-r[0], r[1])):
            extra = f"{str(vt)[0]}   {a:3d}  " if heavy else ""
            print(f"   {str(rnd)[0]}       {chi}    {extra}{isR}")
        print(f"  => round regular classes: {n_round} (R_m is round); chi=2 regular classes: {n_chi2}")
        print(f"     UNIQUE round regular = R_m? {n_round == 1};  UNIQUE chi=2 regular = R_m? {n_chi2 == 1}")

    print("\n" + "=" * 78)
    print("READING")
    print("=" * 78)
    print("""  If round-regular is UNIQUE (=R_m) at each m, then opus-S591's 'LRC=>round' plus
  'regular tight' FORCES the tight regular orbit to be exactly R_m (the AP) -- a rigorous
  pin, closing the 'VT round = interval circulant' qualifier. If chi=2-regular is ALSO
  unique, chi=2 characterizes the AP orbit among ALL regular tournaments (not just the
  round/accessible ones). Any non-R_m round regular, or non-R_m chi=2 regular, would be a
  regular orbit the LRC could in principle reach beyond the AP.""")

if __name__ == "__main__":
    main()
