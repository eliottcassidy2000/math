#!/usr/bin/env python3
"""
paley_cluster_catalan_monad.py
monad-explorer-2026-06-07 (deep-research lane). Confirms the CATALAN law for the
Paley cluster integrals (sharpening HYP-2307):

      A_{2k} = C_k * p^{k+1} + O(p^{k+1/2}),   C_k = k-th Catalan number,

with leading coefficient SIGN = + (the Mobius/inclusion-exclusion sign exactly
cancels the (-1)^k from the k bigons chi(-1)=-1).

The C_k bigon-TREE patterns are the closed length-2k walks that traverse each edge
of a plane tree with k edges exactly twice (once per direction) = Euler tours of
plane trees = Catalan C_k.

General low-memory A_L: fix x_0=0 (translation, factor p); loop the first `nloop`
free coordinates, meshgrid the rest, enforce full distinctness. Choose nloop so the
meshgrid stays <~ 2e7 entries.
"""
import numpy as np

def legendre_array(p):
    chi = np.zeros(p, dtype=np.int64)
    qr = set((x*x) % p for x in range(1, p))
    for d in range(1, p):
        chi[d] = 1 if d in qr else -1
    return chi

def catalan(k):
    from math import comb
    return comb(2*k, k)//(k+1)

def compute_A(p, L, chi, target=2.0e7):
    """A_L = p * sum_{x1..xL distinct, !=0} chi(x1)chi(x2-x1)...chi(xL-x_{L-1}).
       Loop first nloop coords, meshgrid the remaining (L-nloop)."""
    a = np.arange(p)
    # pick nloop so p^(L-nloop) <= target
    nloop = 0
    while p**(L-nloop) > target:
        nloop += 1
    rem = L - nloop
    total = 0
    # iterate loop coords (values for x1..x_nloop)
    import itertools
    for loopvals in itertools.product(range(p), repeat=nloop):
        # x1..x_nloop = loopvals ; must be distinct & !=0
        used = set(loopvals)
        if 0 in used or len(used) != nloop:
            continue
        if rem == 0:
            xs = list(loopvals)
            prod = chi[xs[0] % p]
            for i in range(1, L):
                prod = prod * chi[(xs[i]-xs[i-1]) % p]
            total += int(prod)  # single term
            continue
        grids = np.meshgrid(*([a]*rem), indexing='ij')
        cols = [g.ravel() for g in grids]  # x_{nloop+1} .. x_L
        keep = np.ones(cols[0].shape, dtype=bool)
        for c in cols:
            keep &= (c != 0)
            for v in loopvals:
                keep &= (c != v)
        for i in range(rem):
            for j in range(i+1, rem):
                keep &= (cols[i] != cols[j])
        xs_rem = [c[keep] for c in cols]
        # full x-sequence: [loopvals..., xs_rem...]
        # prod chi(x1) * chi(x2-x1) * ... building incrementally
        # first part among loopvals (scalars), then transitions into arrays
        prev = loopvals[0] if nloop>=1 else None
        # chi(x1)
        if nloop >= 1:
            scal = chi[loopvals[0] % p]
            for i in range(1, nloop):
                scal *= chi[(loopvals[i]-loopvals[i-1]) % p]
            prod = np.full(xs_rem[0].shape, scal, dtype=np.int64)
            prod = prod * chi[(xs_rem[0]-loopvals[-1]) % p]
        else:
            prod = chi[xs_rem[0] % p].copy()
        for i in range(1, rem):
            prod = prod * chi[(xs_rem[i]-xs_rem[i-1]) % p]
        total += int(prod.sum())
    return p * total

def main():
    print("="*84)
    print("THE CATALAN LAW:  A_{2k} = C_k p^{k+1} + O(p^{k+1/2})   (sharpening HYP-2307)")
    print("="*84)
    print(f"Catalan C_1..C_5 = {[catalan(k) for k in range(1,6)]}\n")

    print("L=4 (k=2): A_4/p^3 -> C_2 = 2")
    print(f"  {'p':>3} | {'A_4':>12} | {'A_4/p^3':>9} | {'(2-A_4/p^3)*sqrt p':>18}")
    for p in [7,11,19,23,31,43,47,59,67]:
        chi=legendre_array(p); A4=compute_A(p,4,chi)
        r=A4/p**3
        print(f"  {p:>3} | {A4:>12} | {r:>9.5f} | {(2-r)*p**0.5:>18.4f}")

    print("\nL=6 (k=3): A_6/p^4 -> C_3 = 5  (subleading O(p^{3.5}) oscillates in sign)")
    print(f"  {'p':>3} | {'A_6':>14} | {'A_6/p^4':>9}")
    for p in [7,11,19,23]:
        chi=legendre_array(p); A6=compute_A(p,6,chi)
        print(f"  {p:>3} | {A6:>14} | {A6/p**4:>9.5f}")

    print("\nL=8 (k=4): A_8/p^5 -> C_4 = 14")
    print(f"  {'p':>3} | {'A_8':>16} | {'A_8/p^5':>9}")
    for p in [7,11]:
        chi=legendre_array(p); A8=compute_A(p,8,chi)
        print(f"  {p:>3} | {A8:>16} | {A8/p**5:>9.5f}")

    print("\n=> A_{2k}/p^{2k} -> 0 for all k>=2 (since k+1<2k); only a_2=1 survives; R(p)->e.")

if __name__ == "__main__":
    main()
