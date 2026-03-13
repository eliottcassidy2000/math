#!/usr/bin/env python3
"""
Fast orbit complex computation for larger Paley tournaments.

KEY OPTIMIZATION: Instead of building ALL diff-seqs and then finding orbits,
we build orbit representatives directly by fixing s_1 = 1 (minimum QR element).

Since Z_m acts on diff-seqs by scaling all components, and 1 is the smallest QR
element, the lexicographically smallest member of any orbit starts with s_1 = 1.

This reduces enumeration by factor m.

opus-2026-03-13-S71b
"""
import time
import sys

def get_QR(p):
    return sorted(set(pow(x, 2, p) for x in range(1, p)))

def build_orbit_reps(p, d, report=True):
    """Build canonical orbit representatives for degree-d diff-seqs.
    Fix s_1 = 1, then check if (1, s_2, ..., s_d) is truly the canonical rep."""
    QR = set(pow(x, 2, p) for x in range(1, p))
    QR_list = sorted(QR)
    m = len(QR_list)

    if d == 0:
        return [()]

    results = []

    def is_canonical(seq):
        """Check if seq is the lexicographically smallest in its Z_m orbit."""
        for q in QR_list[1:]:  # skip q=1
            scaled = tuple((q * s) % p for s in seq)
            if scaled < seq:
                return False
        return True

    def backtrack(seq, ps_list, ps_set):
        if len(seq) == d:
            t = tuple(seq)
            if is_canonical(t):
                results.append(t)
            return
        for s in QR_list:
            new_ps = (ps_list[-1] + s) % p
            if new_ps in ps_set:
                continue
            seq.append(s)
            ps_list.append(new_ps)
            ps_set.add(new_ps)
            backtrack(seq, ps_list, ps_set)
            seq.pop()
            ps_list.pop()
            ps_set.remove(new_ps)

    # Start with s_1 = 1 (but we need to check all starting values
    # since the canonical rep might not start with 1 for some orbits)
    # Actually: for s_1 = q, scaling by q^{-1} gives s_1' = 1.
    # And (1, q^{-1}*s_2, ...) <= (q, s_2, ...) iff 1 <= q, which is true.
    # So the canonical rep ALWAYS starts with s_1 = 1.
    # ... unless the orbit has a smaller tuple with a different s_1.
    # But QR_list[0] = 1, so s_1 = 1 is the smallest possible first element.
    # After scaling by q^{-1}, the first element becomes 1, and the rest become q^{-1}*s_i.
    # So the lex-smallest rep always has s_1 = 1.

    # Optimization: only start with s_1 = 1
    s1 = 1
    ps1 = (0 + s1) % p
    backtrack([s1], [0, ps1], {0, ps1})

    return sorted(results)


def compute_face(seq, i, p):
    d = len(seq)
    if i == 0:
        return seq[1:]
    elif i == d:
        return seq[:-1]
    else:
        merged = (seq[i-1] + seq[i]) % p
        return seq[:i-1] + (merged,) + seq[i+1:]


def zm_orbit_class(seq, QR_list, p):
    best = seq
    for q in QR_list:
        scaled = tuple((q * s) % p for s in seq)
        if scaled < best:
            best = scaled
    return best


def sparse_rank_mod(cols_sparse, q):
    """Sparse Gaussian elimination mod q. Returns rank."""
    pivot_col = {}
    reduced = []
    for col in cols_sparse:
        c = dict(col)
        while c:
            r = min(c.keys())
            v = c[r] % q
            if v == 0:
                del c[r]
                continue
            if r in pivot_col:
                pv = reduced[pivot_col[r]]
                pv_val = pv[r]
                inv_pv = pow(pv_val, q - 2, q)
                factor = (v * inv_pv) % q
                for pr, pval in pv.items():
                    c[pr] = (c.get(pr, 0) - factor * pval) % q
                c = {k: v for k, v in c.items() if v % q != 0}
            else:
                pivot_col[r] = len(reduced)
                reduced.append(c)
                break
    return len(reduced)


def compute_orbit_betti(p, q=997, max_d=None):
    m = (p - 1) // 2
    QR = get_QR(p)
    if max_d is None:
        max_d = 2 * m
    print(f"{'='*60}")
    print(f"P_{p} (m={m}), QR = {QR}, mod {q}")
    print(f"{'='*60}\n")

    t0 = time.time()

    # Build orbit representatives directly
    orbit_reps = {}
    for d in range(max_d + 1):
        t1 = time.time()
        orbit_reps[d] = build_orbit_reps(p, d)
        dt = time.time() - t1
        flag = f" ({dt:.1f}s)" if dt > 0.5 else ""
        print(f"d={d}: {len(orbit_reps[d])} orbit reps{flag}", flush=True)
        if dt > 120:
            print(f"  Too slow, stopping at d={d}")
            max_d = d
            break

    print(f"\nBuilt in {time.time()-t0:.1f}s\n")

    QR_set = set(QR)

    # Orbit constraint (sparse)
    C_orb_cols = {}
    C_orb_n_rows = {}
    Omega_orb = {}

    for d in range(max_d + 1):
        reps = orbit_reps[d]
        if d <= 1:
            C_orb_cols[d] = []
            C_orb_n_rows[d] = 0
            Omega_orb[d] = len(reps)
            continue

        junk_orbits = {}
        cols = []
        for sigma in reps:
            col = {}
            for i in range(1, d):
                merged = (sigma[i-1] + sigma[i]) % p
                if merged not in QR_set:
                    face = list(sigma)
                    face[i-1] = merged
                    del face[i]
                    ft = zm_orbit_class(tuple(face), QR, p)
                    if ft not in junk_orbits:
                        junk_orbits[ft] = len(junk_orbits)
                    r = junk_orbits[ft]
                    sign = ((-1) ** i) % q
                    col[r] = (col.get(r, 0) + sign) % q
            col = {k: v for k, v in col.items() if v % q != 0}
            cols.append(col)

        C_orb_cols[d] = cols
        C_orb_n_rows[d] = len(junk_orbits)
        rank_C = sparse_rank_mod(cols, q)
        Omega_orb[d] = len(reps) - rank_C

    print(f"Omega_orb = {[Omega_orb[d] for d in range(max_d+1)]}")
    chi_orb = sum((-1)**d * Omega_orb[d] for d in range(max_d+1))
    print(f"chi_orb = {chi_orb}\n")

    # Stacking ranks
    R_orb = {0: 0}
    for d in range(1, max_d + 1):
        t1 = time.time()
        reps_d = orbit_reps[d]
        reps_dm1 = orbit_reps[d-1]
        reps_dm1_idx = {rep: i for i, rep in enumerate(reps_dm1)}

        n_junk = C_orb_n_rows[d]
        stacked_cols = []
        for j, sigma in enumerate(reps_d):
            col = {}
            if j < len(C_orb_cols[d]):
                for r, v in C_orb_cols[d][j].items():
                    col[r] = v
            for fi in range(d + 1):
                face = compute_face(sigma, fi, p)
                face_canon = zm_orbit_class(face, QR, p)
                if face_canon in reps_dm1_idx:
                    r = n_junk + reps_dm1_idx[face_canon]
                    sign = ((-1) ** fi) % q
                    col[r] = (col.get(r, 0) + sign) % q
            col = {k: v for k, v in col.items() if v % q != 0}
            stacked_cols.append(col)

        rank_stacked = sparse_rank_mod(stacked_cols, q)
        rank_C = sparse_rank_mod(C_orb_cols[d], q) if C_orb_cols[d] else 0
        R_orb[d] = rank_stacked - rank_C

        dt = time.time() - t1
        print(f"  d={d}: R_orb={R_orb[d]} ({dt:.1f}s)", flush=True)

    R_orb[max_d + 1] = 0

    betti_orb = [Omega_orb[d] - R_orb[d] - R_orb[d+1] for d in range(max_d + 1)]
    print(f"\nR_orb = {[R_orb[d] for d in range(max_d + 2)]}")
    print(f"β_orb = {betti_orb}")
    print(f"chi_orb = {sum((-1)**d * b for d, b in enumerate(betti_orb))}")

    beta_m_orb = (m - 3) // 2
    if m <= max_d:
        print(f"\nExpected β_m^orb = (m-3)/2 = {beta_m_orb}")
        if m < len(betti_orb):
            print(f"  β_{m}^orb = {betti_orb[m]} {'✓' if betti_orb[m] == beta_m_orb else '✗'}")
        if m+1 < len(betti_orb):
            print(f"  β_{m+1}^orb = {betti_orb[m+1]} {'✓' if betti_orb[m+1] == beta_m_orb else '✗'}")

    print(f"\n{'='*60}\n")
    return betti_orb


if __name__ == "__main__":
    compute_orbit_betti(7)
    compute_orbit_betti(11)
    compute_orbit_betti(19)
    # Try P_23 if P_19 works
    if len(sys.argv) > 1 and sys.argv[1] == '--p23':
        compute_orbit_betti(23)
