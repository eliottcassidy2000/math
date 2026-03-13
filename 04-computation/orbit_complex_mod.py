#!/usr/bin/env python3
"""
Orbit complex computation for Paley tournaments using modular arithmetic.
Uses Gaussian elimination mod a prime q for fast rank computation.

opus-2026-03-13-S71b
"""
import time

def get_QR(p):
    return sorted(set(pow(x, 2, p) for x in range(1, p)))

def build_diffseqs(p, d):
    QR = set(pow(x, 2, p) for x in range(1, p))
    QR_list = sorted(QR)
    results = []
    if d == 0:
        return [()]
    def backtrack(seq, ps_list, ps_set):
        if len(seq) == d:
            results.append(tuple(seq))
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
    backtrack([], [0], {0})
    return results

def zm_orbit_class(seq, QR_list, p):
    best = seq
    for q in QR_list:
        scaled = tuple((q * s) % p for s in seq)
        if scaled < best:
            best = scaled
    return best

def compute_face(seq, i, p):
    d = len(seq)
    if i == 0:
        return seq[1:]
    elif i == d:
        return seq[:-1]
    else:
        merged = (seq[i-1] + seq[i]) % p
        return seq[:i-1] + (merged,) + seq[i+1:]

def sparse_rank_mod(cols_sparse, q, n_rows):
    """
    Gaussian elimination on sparse columns mod q.
    cols_sparse: list of dicts {row: value}
    Returns rank.
    """
    pivot_col = {}  # row -> column index
    reduced = []

    for col in cols_sparse:
        c = dict(col)
        while c:
            # Find leading row
            r = min(c.keys())
            v = c[r] % q
            if v == 0:
                del c[r]
                continue
            if r in pivot_col:
                # Reduce using pivot
                pv = reduced[pivot_col[r]]
                pv_val = pv[r]
                # Subtract (v / pv_val) * pivot from c
                inv_pv = pow(pv_val, q - 2, q)
                factor = (v * inv_pv) % q
                for pr, pval in pv.items():
                    c[pr] = (c.get(pr, 0) - factor * pval) % q
                # Clean zeros
                c = {k: v for k, v in c.items() if v % q != 0}
            else:
                pivot_col[r] = len(reduced)
                reduced.append(c)
                break

    return len(reduced)

def main():
    q = 997  # prime for modular arithmetic

    for p in [7, 11]:
        m = (p - 1) // 2
        QR = get_QR(p)
        print(f"{'='*60}")
        print(f"P_{p} (m={m}), QR = {QR}, working mod {q}")
        print(f"{'='*60}\n")

        t0 = time.time()

        # Build diff-seqs and orbits
        Ad = {}
        orbit_reps = {}
        for d in range(2*m + 1):
            Ad[d] = build_diffseqs(p, d)
            if d == 0:
                orbit_reps[d] = [()]
            else:
                seen = set()
                reps = []
                for seq in Ad[d]:
                    canon = zm_orbit_class(seq, QR, p)
                    if canon not in seen:
                        seen.add(canon)
                        reps.append(canon)
                orbit_reps[d] = sorted(reps)
            print(f"d={d}: |A|={len(Ad[d])}, orbits={len(orbit_reps[d])}", flush=True)

        print(f"\nDiff-seqs built in {time.time()-t0:.1f}s\n")

        QR_set = set(QR)

        # Build orbit constraint columns (sparse)
        C_orb_cols = {}  # d -> list of sparse column dicts
        C_orb_n_rows = {}  # d -> number of junk orbit rows
        Omega_orb = {}

        for d in range(2*m + 1):
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
            rank_C = sparse_rank_mod(cols, q, len(junk_orbits))
            Omega_orb[d] = len(reps) - rank_C

        print(f"Omega_orb = {[Omega_orb[d] for d in range(2*m+1)]}")
        chi_orb = sum((-1)**d * Omega_orb[d] for d in range(2*m+1))
        print(f"chi_orb = {chi_orb}\n")

        # Build orbit boundary sparse columns and compute stacking rank
        R_orb = {0: 0}
        for d in range(1, 2*m + 1):
            t1 = time.time()
            reps_d = orbit_reps[d]
            reps_dm1 = orbit_reps[d-1]
            reps_dm1_idx = {rep: i for i, rep in enumerate(reps_dm1)}

            # Build boundary columns (sparse): B maps d-orbits to (d-1)-orbits
            # For stacking, we need [C_d; B_d] where C_d constrains degree-d orbits
            # and B_d maps degree-d orbits to degree-(d-1) orbits.
            # But these have different row spaces!
            #
            # Stacking trick: rank(∂ | Omega) = rank([C_d; B_d]) - rank(C_d)
            # where C_d has rows indexed by "junk orbits of degree d"
            # and B_d has rows indexed by "orbit reps of degree d-1".
            # We stack them vertically: C_d rows first, then B_d rows (shifted).

            n_junk = C_orb_n_rows[d]
            stacked_cols = []
            for j, sigma in enumerate(reps_d):
                col = {}
                # Constraint rows (from C_orb_cols[d])
                if j < len(C_orb_cols[d]):
                    for r, v in C_orb_cols[d][j].items():
                        col[r] = v
                # Boundary rows (shifted by n_junk)
                for fi in range(d + 1):
                    face = compute_face(sigma, fi, p)
                    face_canon = zm_orbit_class(face, QR, p)
                    if face_canon in reps_dm1_idx:
                        r = n_junk + reps_dm1_idx[face_canon]
                        sign = ((-1) ** fi) % q
                        col[r] = (col.get(r, 0) + sign) % q
                col = {k: v for k, v in col.items() if v % q != 0}
                stacked_cols.append(col)

            rank_stacked = sparse_rank_mod(stacked_cols, q, n_junk + len(reps_dm1))
            rank_C = sparse_rank_mod(C_orb_cols[d], q, n_junk) if C_orb_cols[d] else 0
            R_orb[d] = rank_stacked - rank_C

            dt = time.time() - t1
            print(f"  d={d}: R_orb={R_orb[d]} ({dt:.1f}s)", flush=True)

        R_orb[2*m + 1] = 0

        betti_orb = [Omega_orb[d] - R_orb[d] - R_orb[d+1] for d in range(2*m + 1)]
        print(f"\nR_orb = {[R_orb[d] for d in range(2*m + 2)]}")
        print(f"β_orb = {betti_orb}")
        print(f"chi_orb = {sum((-1)**d * b for d, b in enumerate(betti_orb))}")

        # Expected
        beta_m_orb = (m - 3) // 2
        print(f"\nExpected β_m^orb = (m-3)/2 = {beta_m_orb}")
        if m >= 3:
            print(f"  β_{m}^orb = {betti_orb[m]} {'✓' if betti_orb[m] == beta_m_orb else '✗'}")
            print(f"  β_{m+1}^orb = {betti_orb[m+1]} {'✓' if betti_orb[m+1] == beta_m_orb else '✗'}")

        # Full Betti from orbit: β^{(0)} = m * β^{orb}
        beta_full = [m * b if d > 0 else 1 for d, b in enumerate(betti_orb)]
        # Actually β_0^{(0)} = 1, and for d>0: β_d^{(0)} = m * β_d^{orb}
        # Hmm, β_0^{orb} = 1 and β_0^{(0)} = 1. But m * 1 = m ≠ 1.
        # The issue: d=0 has Omega_0 = 1, not m. The Z_m action is NOT free at d=0.
        beta_full[0] = 1
        print(f"\nFull k=0 Betti (m * β_orb): {beta_full}")
        print(f"  β_m = m*(m-3)/2 = {m * beta_m_orb}")

        print(f"\n{'='*60}\n")


if __name__ == "__main__":
    main()
