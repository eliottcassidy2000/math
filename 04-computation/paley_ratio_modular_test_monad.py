"""
paley_ratio_modular_test_monad.py  (monad-explorer)

TEST the open question posed in 07-reflections/the-tessellation.md (Layer 6, opus-S131):

    "The sequence 1, 9, 1729, ... is the canonical tournament ratio
     r(p) = H(T_p)/|Aut(T_p)| at Paley primes. 1729 = 7*13*19 = j(i)+1 appeared
     at p=11, the first Paley prime where X_0(p) has genus 1.  Whatever comes next
     will tell us whether this sequence has modular significance."

KEY OBSERVATION the reflection missed:  the NEXT Paley prime is p=19 (19 = 3 mod 4),
and X_0(19) ALSO has genus 1 (X_0(23) has genus 2).  So p=19 -- already computed
(HYP-1266) -- is the SHARPEST test of the genus-tied modular hypothesis, better than
the p=23 the reflection anticipated.

This script:
  (1) builds the Paley tournament T_p (p = 3 mod 4),
  (2) counts directed Hamiltonian paths H(T_p) by Held-Karp (start fixed at 0 via
      vertex-transitivity: H = p * #paths-from-0), VALIDATED against the known
      values H(T_7)=189, H(T_11)=95095, H(T_19)=1172695746915,
  (3) forms r(p)=H/|Aut|, factors it, and reports each prime factor's class mod 3
      (split / inert / ramified in Q(sqrt(-3))) and proximity to j-values,
  (4) prints the verdict on whether 1729's "completely split in Q(sqrt-3) = j(i)+1"
      structure persists.
"""
import sys, time
import numpy as np
from sympy import factorint, isprime

def qr_set(p):
    return {(x*x) % p for x in range(1, p)}

def paley_adj(p):
    assert p % 4 == 3, "Paley tournament needs p = 3 mod 4"
    Q = qr_set(p)
    A = np.zeros((p, p), dtype=np.int64)
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in Q:
                A[i, j] = 1
    return A

def count_ham_paths(A):
    """Directed Hamiltonian path count via Held-Karp, fixing the start at vertex 0
    and multiplying by n (vertex-transitive tournaments only)."""
    n = A.shape[0]
    full = (1 << n) - 1
    # dp[mask] is a length-n int64 vector: dp[mask][v] = # directed paths starting
    # at 0, covering exactly `mask`, ending at v.  mask always contains bit 0.
    dp = np.zeros((1 << n, n), dtype=object)  # object dtype -> arbitrary precision
    dp[1][0] = 1  # mask = {0}, end at 0
    Acol = [A[:, w] for w in range(n)]  # in-neighbour indicators
    for mask in range(1, 1 << n):
        if not (mask & 1):
            continue
        vec = dp[mask]
        if not vec.any():
            continue
        for w in range(n):
            if mask & (1 << w):
                continue
            # contribution = sum_{v in mask, v->w} dp[mask][v]
            c = 0
            col = Acol[w]
            m = mask
            while m:
                v = (m & -m).bit_length() - 1
                if col[v]:
                    c += vec[v]
                m &= m - 1
            if c:
                dp[mask | (1 << w)][w] += c
    paths_from_0 = sum(dp[full])
    return n * paths_from_0

def classify_mod3(q):
    if q == 3:
        return "ramified"
    return "split" if q % 3 == 1 else "inert"

def report(p, H):
    aut = p * (p - 1) // 2  # |Aut(Paley T_p)| = p(p-1)/2
    assert H % aut == 0, f"H not divisible by |Aut| at p={p}"
    r = H // aut
    fr = factorint(r)
    print(f"\n=== p={p}  (X_0({p}) genus {'1' if p in (11,17,19) else '?'}) ===")
    print(f"  H(T_{p})     = {H}   factor={factorint(H)}")
    print(f"  |Aut(T_{p})| = {aut}")
    print(f"  ratio r({p}) = {r}")
    print(f"  factor r({p}) = {fr}")
    classes = {q: classify_mod3(q) for q in fr}
    print(f"  prime classes in Q(sqrt-3): {classes}")
    all_split = all(c != "inert" for c in classes.values())
    smooth = max(fr) <= 100 if fr else True
    print(f"  -> completely split in Q(sqrt-3)? {all_split}")
    print(f"  -> smooth (all primes <=100)?    {smooth}")
    if r == 1729:
        print("  -> r == 1729 == j(i)+1 (j(i)=1728)")
    return r, all_split, smooth

if __name__ == "__main__":
    # known values for validation / speed (avoid recomputing the slow ones unless asked)
    KNOWN = {3: 3, 7: 189, 11: 95095, 19: 1172695746915}
    primes = [3, 7, 11, 19]
    recompute = set()
    if len(sys.argv) > 1:
        # e.g. "7 11 19" to force recompute, or "23" to attempt the new term
        primes = [int(x) for x in sys.argv[1:]]
        recompute = set(primes)

    results = {}
    for p in primes:
        if p in KNOWN and p not in recompute:
            H = KNOWN[p]
            print(f"[using known H(T_{p}) = {H}]")
        else:
            t0 = time.time()
            print(f"[computing H(T_{p}) via Held-Karp ...]", flush=True)
            H = count_ham_paths(paley_adj(p))
            print(f"[H(T_{p}) = {H}  in {time.time()-t0:.1f}s]")
            if p in KNOWN:
                ok = (H == KNOWN[p])
                print(f"[validation vs known: {'PASS' if ok else 'FAIL'}]")
        results[p] = report(p, H)

    print("\n=== VERDICT ===")
    for p in primes:
        if p in results:
            r, allsplit, smooth = results[p]
            print(f"  p={p}: completely-split={allsplit}, smooth={smooth}")
