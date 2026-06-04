"""
verify_t11_betti_s_monad.py  (monad-compute)

INV-143 next step: independently recompute the FULL GLMY path-homology Betti
numbers of the Paley tournament T_11 from scratch and check against the cached
claim KNOWN_BETTI[11] = [1,0,0,0,0,5,15,0,0,0,0] (verify beta_5=5, beta_6=15).

Guards MISTAKE-020: max_degree = n-1 = 10 (full complex).

Output is LINE-BUFFERED and TIMED per phase / per eigenspace so partial
verified results are captured even if the full run is long. Each eigenspace's
boundary ranks are an independently-checkable artifact.
"""
import sys, os, time
sys.stdout.reconfigure(line_buffering=True)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from circulant_homology import PaleyHomology, find_nth_root_of_unity

P = 11
EXPECTED_OMEGA = [1, 5, 20, 70, 205, 460, 700, 690, 450, 180, 30]
EXPECTED_BETTI = [1, 0, 0, 0, 0, 5, 15, 0, 0, 0, 0]
MAXD = P - 1  # 10, full complex

print("=" * 64)
print(f"monad-compute: from-scratch Betti verification, Paley T_{P}")
print("=" * 64)

h = PaleyHomology(p=P)
print(f"n={h.n}  |S|={len(h.S)}  root-field prime={h.prime}  max_degree={MAXD}")
print(f"S (QR mod {P}) = {sorted(h.S)}")
sys.stdout.flush()

# ---- Phase 1: enumerate allowed paths + Omega dims (from scratch) ----
t0 = time.time()
h._ensure_enumerated(MAXD + 1)
raw = [len(h._diff_seqs.get(m, [])) for m in range(MAXD + 1)]
print(f"\n[enumeration] raw |A_m| (directed path counts) = {raw}   [{time.time()-t0:.1f}s]")
sys.stdout.flush()

t0 = time.time()
omega = h.omega_dims(max_degree=MAXD, use_cache=False, verbose=False)
chi = sum((-1) ** m * d for m, d in enumerate(omega))
print(f"[omega] recomputed Omega dims = {omega}   [{time.time()-t0:.1f}s]")
print(f"[omega] expected             = {EXPECTED_OMEGA}")
print(f"[omega] MATCH={omega == EXPECTED_OMEGA}   chi={chi} (expect 1)")
sys.stdout.flush()

# ---- Phase 2: boundary ranks per eigenspace (timed, incremental) ----
omega_p = find_nth_root_of_unity(h.n, h.prime)
boundary_ranks = {}
print(f"\n[boundary] computing rank(d_m^(k)) for k=0..{h.n-1}, m=0..{MAXD+1}")
sys.stdout.flush()
for k in range(h.n):
    tk = time.time()
    omega_k = pow(omega_p, k, h.prime)
    ranks_k = [h._boundary_rank_k(m, omega_k) for m in range(MAXD + 2)]
    boundary_ranks[k] = ranks_k
    print(f"  k={k:2d}: ranks={ranks_k}   [{time.time()-tk:.1f}s]")
    sys.stdout.flush()

# ---- Phase 3: assemble Betti ----
betti = []
for m in range(MAXD + 1):
    tot = 0
    for k in range(h.n):
        ker_m = omega[m] - boundary_ranks[k][m]
        im_next = boundary_ranks[k][m + 1]
        tot += ker_m - im_next
    betti.append(tot)

chi_betti = sum((-1) ** m * b for m, b in enumerate(betti))
print(f"\n[betti] recomputed = {betti}")
print(f"[betti] expected   = {EXPECTED_BETTI}")
print(f"[betti] MATCH={betti == EXPECTED_BETTI}")
print(f"[betti] beta_5={betti[5]} (expect 5)   beta_6={betti[6]} (expect 15)")
print(f"[euler] chi from Betti={chi_betti}  chi from Omega={chi}  MATCH={chi_betti==chi}")

ok = (omega == EXPECTED_OMEGA and betti == EXPECTED_BETTI
      and chi == 1 and chi_betti == 1)
print("\n" + "=" * 64)
print(f"OVERALL: {'ALL CHECKS PASS' if ok else 'MISMATCH — INVESTIGATE'}")
print("=" * 64)
