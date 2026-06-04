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
import sys, os, time, json
sys.stdout.reconfigure(line_buffering=True)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "..", "05-knowledge", "results")
JSON_PATH = os.path.join(RESULTS_DIR, "verify_t11_betti_s_monad_ranks.json")

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
# Resume support: load any eigenspaces already computed in a prior run.
if os.path.exists(JSON_PATH):
    with open(JSON_PATH) as f:
        saved = json.load(f)
    for ks, rk in saved.get("boundary_ranks", {}).items():
        boundary_ranks[int(ks)] = rk
    if boundary_ranks:
        print(f"[resume] loaded eigenspaces k={sorted(boundary_ranks)} from {JSON_PATH}")
        sys.stdout.flush()

for k in range(h.n):
    if k in boundary_ranks:
        print(f"  k={k:2d}: ranks={boundary_ranks[k]}   [cached/resumed]")
        sys.stdout.flush()
        continue
    # Clear the cross-eigenspace basis cache: entries keyed by a different
    # omega_k are never reused and only bloat memory / slow allocation.
    h._omega_basis_cache.clear()
    tk = time.time()
    omega_k = pow(omega_p, k, h.prime)
    ranks_k = [h._boundary_rank_k(m, omega_k) for m in range(MAXD + 2)]
    boundary_ranks[k] = ranks_k
    print(f"  k={k:2d}: ranks={ranks_k}   [{time.time()-tk:.1f}s]")
    sys.stdout.flush()
    # Persist incrementally so each eigenspace can be committed as it lands.
    with open(JSON_PATH, "w") as f:
        json.dump({"p": P, "expected_omega": EXPECTED_OMEGA, "omega": omega,
                   "boundary_ranks": {str(kk): boundary_ranks[kk]
                                      for kk in sorted(boundary_ranks)}},
                  f, indent=1)

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
# Euler-Poincare for this eigenspace-decomposed complex:
# by THM-125 each of the n eigenspaces carries a FULL copy of Omega_m
# (dim Omega_m^(k) = omega[m] for all k), so the total chain dimension at
# degree m is n*omega[m] and chi_full = n * sum (-1)^m omega[m] = n * chi.
# Hence chi_betti must equal n*chi (NOT chi). The single-copy chi=1; n*chi=11.
chi_full = h.n * chi
euler_ok = (chi_betti == chi_full)
print(f"[euler] chi from Betti={chi_betti}  n*chi_Omega={chi_full} (n={h.n}, single-copy chi={chi})  MATCH={euler_ok}")

ok = (omega == EXPECTED_OMEGA and betti == EXPECTED_BETTI
      and chi == 1 and euler_ok)
print("\n" + "=" * 64)
print(f"OVERALL: {'ALL CHECKS PASS' if ok else 'MISMATCH — INVESTIGATE'}")
print("=" * 64)
