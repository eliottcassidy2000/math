#!/usr/bin/env python3
"""
Orbit complex for P_19 (m=9) and P_23 (m=11).
Predictions: β_m^{orb} = (m-3)/2 = 3 (P_19) and 4 (P_23).

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

def compute_orbit_betti(p, q=997):
    m = (p - 1) // 2
    QR = get_QR(p)
    print(f"{'='*60}")
    print(f"P_{p} (m={m}), QR = {QR}, working mod {q}")
    print(f"{'='*60}\n")

    t0 = time.time()

    # Build diff-seqs and orbits
    Ad = {}
    orbit_reps = {}
    total_cells = 0
    for d in range(2*m + 1):
        t1 = time.time()
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
        total_cells += len(orbit_reps[d])
        dt = time.time() - t1
        flag = f" ({dt:.1f}s)" if dt > 1 else ""
        print(f"d={d}: |A|={len(Ad[d])}, orbits={len(orbit_reps[d])}{flag}", flush=True)

    print(f"\nTotal orbit cells: {total_cells}")
    print(f"Diff-seqs built in {time.time()-t0:.1f}s\n")

    QR_set = set(QR)

    # Orbit constraint columns (sparse)
    C_orb_cols = {}
    C_orb_n_rows = {}
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

    # Stacking ranks
    R_orb = {0: 0}
    for d in range(1, 2*m + 1):
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

    beta_m_orb = (m - 3) // 2
    print(f"\nExpected β_m^orb = (m-3)/2 = {beta_m_orb}")
    print(f"  β_{m}^orb = {betti_orb[m]} {'✓' if betti_orb[m] == beta_m_orb else '✗'}")
    print(f"  β_{m+1}^orb = {betti_orb[m+1]} {'✓' if betti_orb[m+1] == beta_m_orb else '✗'}")

    print(f"\nFull k=0 Betti: β_m = m*(m-3)/2 = {m * beta_m_orb}")
    print(f"Total β_m = m*(m-3)/2 = {m * beta_m_orb}")
    print(f"Total β_{{m+1}} = m*(m-3)/2 + (p-1)*1 = {m * beta_m_orb + p - 1} = C(m+1,2) = {m*(m+1)//2}")
    print(f"\n{'='*60}\n")

    return betti_orb

if __name__ == "__main__":
    compute_orbit_betti(7)
    compute_orbit_betti(11)
    compute_orbit_betti(19)
    # P_23 may be too large — try if P_19 works
    compute_orbit_betti(23)
