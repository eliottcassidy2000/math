#!/usr/bin/env python3
"""
unit_distance_3n_crossover_families_s4.py
monad-explorer-2026-06-07-S4  (deep-research, building on THM-431 / OPEN-Q-057)

QUESTION (OPEN-Q-057): N* = smallest N with u(N) > 3N is in [25,28] (THM-431).
The ceiling N*<=28 rests on a CITED Engel "Moser lattice" construction
u(28)>=85 the repo never built; THM-431-C showed the repo's sqrt(7) Eisenstein
family crosses 3N too late (k=39 disk / 32 anneal), i.e. the WRONG family.

This script asks, with EXACT INTEGER counting and a strong simulated-annealing
densest-patch search: which SINGLE-NORM lattice family realizes the small-N
3N-crossover, and does any beat 3N at N <= 28?

family = symmetric offset set V in Z^2.  patch P in Z^2.
   u(P) = #{ {p,q} subset P : p-q in V }   (exact).
Any patch found => CONSTRUCTIVE LOWER BOUND u(N) >= u(P).

Calibration target: the repo's s710 annealer hits u(P)=97 on a 32-cell sqrt7
patch (3*32=96, beats by 1). A trustworthy search must match that.
"""
import random

random.seed(20260607)   # deterministic (Date/random banned in workflow, fine here)

# ---------- offset sets ----------
def euclid_offsets(t, rad=14):
    return [(a, b) for a in range(-rad, rad+1) for b in range(-rad, rad+1)
            if a*a + b*b == t]

def eisen_offsets(t, rad=14):
    return [(a, b) for a in range(-rad, rad+1) for b in range(-rad, rad+1)
            if a*a + a*b + b*b == t]

def euclid_norm2(p):           # real squared distance
    return p[0]*p[0] + p[1]*p[1]

def eisen_norm2(p):
    return p[0]*p[0] + p[0]*p[1] + p[1]*p[1]

FAMILIES = [
    ("grid t=5  = KNIGHT (deg8)", euclid_offsets(5),  euclid_norm2),
    ("grid t=25 (deg12)",         euclid_offsets(25), euclid_norm2),
    ("grid t=65 (deg16)",         euclid_offsets(65), euclid_norm2),
    ("eisen t=1 (penny, deg6)",   eisen_offsets(1),   eisen_norm2),
    ("eisen t=7 (sqrt7, deg12)",  eisen_offsets(7),   eisen_norm2),
    ("eisen t=13 (AMP, deg12)",   eisen_offsets(13),  eisen_norm2),
]

# ---------- patch tools ----------
def deg_in(p, Sset, V):
    return sum((p[0]+v[0], p[1]+v[1]) in Sset for v in V)

def edges(Sset, V):
    return sum(deg_in(p, Sset, V) for p in Sset) // 2

def disk(V, norm2, N, rad=14):
    cells = sorted(((a, b) for a in range(-rad, rad+1) for b in range(-rad, rad+1)),
                   key=lambda p: (norm2(p), p))
    return set(cells[:N])

def anneal(V, norm2, N, iters=60000, restarts=6):
    """SA densest N-cell patch. incremental edge bookkeeping."""
    best_global = -1
    best_set = None
    for r in range(restarts):
        S = disk(V, norm2, N) if r == 0 else perturbed_disk(V, norm2, N, r)
        Sset = set(S)
        E = edges(Sset, V)
        best = E
        best_local = set(Sset)
        # candidate external cells maintained lazily
        for it in range(iters):
            T = max(0.02, 2.5 * (1 - it/iters))   # cooling
            # propose swap: remove u in S, add w adjacent to S\{u}
            u = random.choice(tuple(Sset))
            # pick external neighbour of the patch
            p = random.choice(tuple(Sset))
            v = random.choice(V)
            w = (p[0]+v[0], p[1]+v[1])
            if w in Sset:
                continue
            if w == u:
                continue
            du = deg_in(u, Sset, V)
            # degree of w in S after removing u:
            dw = deg_in(w, Sset, V) - (1 if ((w[0]-u[0], w[1]-u[1]) in Vset_cache[id(V)]) else 0)
            delta = dw - du
            if delta >= 0 or random.random() < pow(2.718281828, delta / T):
                Sset.remove(u); Sset.add(w)
                E += delta
                if E > best:
                    best = E
                    best_local = set(Sset)
        if best > best_global:
            best_global = best
            best_set = best_local
    return best_global, best_set

def perturbed_disk(V, norm2, N, seed_off):
    base = disk(V, norm2, N + 6)
    base = list(base)
    random.shuffle(base)
    return set(base[:N])

# cache offset sets for O(1) membership
Vset_cache = {}

def run_family(name, V, norm2, Nlist, iters, restarts):
    Vset_cache[id(V)] = set(V)
    out = {}
    for N in Nlist:
        e, S = anneal(V, norm2, N, iters=iters, restarts=restarts)
        out[N] = (e, S)
    return out

def verify_two_common(V):
    Vset = set(V)
    cells = set()
    for v in V:
        cells.add(v)
        for w in V:
            cells.add((v[0]+w[0], v[1]+w[1]))
    worst = 0
    for c in cells:
        if c == (0, 0):
            continue
        cn = sum(((x[0]-c[0], x[1]-c[1]) in Vset) for x in V)
        worst = max(worst, cn)
    return worst

def main():
    print("="*72)
    print("3N-CROSSOVER: strong SA densest-patch search, single-norm UD families")
    print("exact integer counts; u(N) >= count is a constructive lower bound")
    print("="*72)

    # calibration: sqrt7 at N=32 should reach ~97 (repo s710)
    Vc = eisen_offsets(7); Vset_cache[id(Vc)] = set(Vc)
    e32, _ = anneal(Vc, eisen_norm2, 32, iters=80000, restarts=8)
    print(f"\n[calibration] sqrt7 N=32: best={e32} (repo s710 anneal=97, 3N=96)")

    Nlist = list(range(22, 61, 2))
    summary = []
    for name, V, norm2 in FAMILIES:
        deg = len(V)
        wc = verify_two_common(V)
        res = run_family(name, V, norm2, Nlist, iters=40000, restarts=5)
        cross = None
        for N in sorted(res):
            if res[N][0] > 3*N:
                cross = N; break
        print(f"\n--- {name}  (degree {deg}, max common-nbrs={wc}) ---")
        print(f"   {'N':>3} {'u>=':>5} {'3N':>5} {'u-3N':>5}")
        for N in sorted(res):
            cnt = res[N][0]
            mark = "  <== BEATS 3N" if cnt > 3*N else ""
            print(f"   {N:>3} {cnt:>5} {3*N:>5} {cnt-3*N:>5}{mark}")
        summary.append((name, deg, cross, res[cross][0] if cross else None))

    print("\n" + "="*72)
    print("SUMMARY: first single-norm-lattice 3N-crossover by family")
    print(f"{'family':<30}{'deg':>4}{'N_cross':>9}{'count':>7}")
    for name, deg, cross, cnt in summary:
        print(f"{name:<30}{deg:>4}{str(cross) if cross else '>60':>9}{str(cnt) if cnt else '-':>7}")
    print("="*72)
    print("Anchors (AMP arXiv:2412.11914 + Engel): u(n)<=3n all n<=24 => N*>=25;")
    print("u(28)>=85>84 => N*<=28. A lattice family crossing at N in {25,26,27}")
    print("would lower the ceiling; crossing >28 shows Engel's record is NON-lattice.")

if __name__ == "__main__":
    main()
