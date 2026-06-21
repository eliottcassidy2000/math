#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_support_alternation_kps.py   (kind-pasteur 2026-06-21, HYP-2724)

The CORRECT weight-enumerator decomposition of corr(E), using the ZERO-PADDED
measure kernel K (length k) -- NOT the bare active-coord Q (the THM-538 trap).

corr(E) = Sum_{0!=n in Lambda(E)} K(n).  K(n) depends only on (multiset of nonzero
coeffs, k) since the product over coords is order-free and zeros give chat(0,T)=(1-|T|/7).
We compute the per-support partial sums W_s(E)=Sum_{support-s relations}K and the
truncated reconstruction corr(<=s_max), comparing to the EXACT corr=measS7-iid.

GOAL: (1) does low-support (<=s_max) reconstruct corr with the CORRECT kernel?
      (2) the SIGN-by-support structure (the signed Erdos-Turan seam).
"""
import itertools, os, importlib.util
from fractions import Fraction as Fr

_d = os.path.dirname(__file__)
def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, os.path.join(_d, path))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
rcm = _load("rcm", "lrc_q108_relation_code_mds_kps.py")     # measS7, iid_k
wev = _load("wev", "lrc_q108_weight_enumerator_validate_kpswf2.py")  # Kk (zero-padded kernel)
measS7, iid_k, Kk = rcm.measS7, rcm.iid_k, wev.Kk

_Kcache = {}
def Kpad(nonzero_coeffs, k):
    """zero-padded K: multiset of nonzero coeffs + (k-len) zeros."""
    key = (tuple(sorted(nonzero_coeffs)), k)
    if key in _Kcache: return _Kcache[key]
    vec = list(nonzero_coeffs) + [0]*(k - len(nonzero_coeffs))
    v = Kk(vec).real
    _Kcache[key] = v; return v

def relations_by_support(E, B=2, max_support=6):
    """yield (support, nonzero_coeffs) for each PRIMITIVE relation, |coef|<=B.
       counts n and -n once (canonical)."""
    from math import gcd
    from functools import reduce
    E = [int(e) for e in E]; k = len(E)
    seen = set()
    for s in range(2, max_support+1):
        for combo in itertools.combinations(range(k), s):
            for coefs in itertools.product(range(-B, B+1), repeat=s):
                if any(c == 0 for c in coefs): continue
                if sum(c*E[i] for c,i in zip(coefs, combo)) != 0: continue
                g = reduce(gcd, [abs(c) for c in coefs]); prim = tuple(c//g for c in coefs)
                if prim[0] < 0: prim = tuple(-c for c in prim)
                key = (combo, prim)
                if key in seen: continue
                seen.add(key)
                yield s, prim

def analyze(E, B=2, max_support=6):
    k = len(E)
    corr_exact = float(measS7(E) - iid_k(k))
    W = {s: 0.0 for s in range(2, max_support+1)}
    cnt = {s: 0 for s in range(2, max_support+1)}
    for s, prim in relations_by_support(E, B, max_support):
        W[s] += 2*Kpad(prim, k)   # n and -n both contribute (K(-n)=conj K(n), real part doubles)
        cnt[s] += 1
    return corr_exact, W, cnt

def main():
    battery = {
        "consec {0..7}":   list(range(8)),
        "consec {1..8}":   list(range(1,9)),
        "AP step2":        [0,2,4,6,8,10,12,14],
        "Sidon":           [0,1,3,7,12,20,30,44],
        "two-block":       [0,1,2,3,40,41,42,43],
        "wide near-consec":[0,1,2,3,4,5,6,40],
    }
    print("="*86)
    print("CORRECT (zero-padded K) weight-enumerator decomposition  corr = Sum_s W_s")
    print("(B=|coef|<=2, support<=6; W_s = signed K-weighted support-s mass)")
    print("="*86)
    print(f"{'set':<18}{'corr':>9} | {'W2':>9}{'W3':>9}{'W4':>9}{'W5':>9}{'W6':>9} | {'sum<=6':>9}{'frac':>7}")
    for name, E in battery.items():
        ce, W, cnt = analyze(E, B=2, max_support=6)
        tot = sum(W.values())
        frac = tot/ce if abs(ce)>1e-9 else float('nan')
        print(f"{name:<18}{ce:>9.4f} | {W[2]:>9.4f}{W[3]:>9.4f}{W[4]:>9.4f}{W[5]:>9.4f}{W[6]:>9.4f} | {tot:>9.4f}{frac:>7.2f}")
    print()
    print("SIGN-BY-SUPPORT (the signed Erdos-Turan seam): is there an alternation W2<0,W3>0,...?")
    print("FRACTION = (low-support<=6 truncation) / (exact corr). If ~1, low-support reconstructs corr")
    print("with the CORRECT kernel (vindicating the weight-enumerator route); if <1, tail (height>2 or")
    print("support>6) matters and the bound needs the conditional-convergence/R6 machinery.")
    print("\nDONE.")

if __name__ == "__main__":
    main()
