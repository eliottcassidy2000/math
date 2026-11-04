#!/usr/bin/env python3
"""audit_A (e), independent brute force.  Shares NO code with cpq_defect_census.py.
For H = <c_P c_Q> of type (p, q, 1^k) enumerate ALL 2^c H-fixed tournaments and ALL 2^c H-fixed graphs as c-bit
words (one bit per H-orbit of pairs), compute
  * the N-action on words directly (N = closure of the normaliser generators; its order is checked against
    pq k! phi(pq) = |C_{S_n}(h)| * #{generators of H}, which is the order of the full normaliser),
  * N-orbit representatives as min over N of g.w (vectorised with numpy),
  * N-evenness of a graph word directly from the definition: no g in N with g.w = w reverses an odd number of
    edges (reversal parity of g on w = popcount(w & rho_g) with rho_g the per-orbit reversal-parity mask),
  * isomorphism types by nauty certificates (directed for tournaments) and sign-evenness from nauty generators,
and then tau_T, tau_E, W, E_T, E_even, E_wit, R, RG from the definitions of THM-4728, and the identity.
Usage: python3 audit_e1_brute.py p q k
"""
import sys, math, time
import numpy as np
import pynauty

def units(m): return [a for a in range(1, m) if math.gcd(a, m) == 1]

def closure(gens, n):
    ident = tuple(range(n)); seen = {ident}; frontier = [ident]
    while frontier:
        nxt = []
        for g in frontier:
            for s in gens:
                x = tuple(s[g[i]] for i in range(n))
                if x not in seen: seen.add(x); nxt.append(x)
        frontier = nxt
    return list(seen)

def main(p, q, k):
    t0 = time.time()
    n = p + q + k
    h = list(range(n))
    for i in range(p): h[i] = (i + 1) % p
    for j in range(q): h[p + j] = p + (j + 1) % q
    h = tuple(h)
    gens = [h]
    for (off, m) in ((0, p), (p, q)):
        for a in units(m):
            g = list(range(n))
            for i in range(m): g[off + i] = off + (a * i) % m
            gens.append(tuple(g))
    if k >= 2:
        g = list(range(n)); g[p + q], g[p + q + 1] = p + q + 1, p + q; gens.append(tuple(g))
        g = list(range(n))
        for i in range(k): g[p + q + i] = p + q + (i + 1) % k
        gens.append(tuple(g))
    N = closure(gens, n)
    phi = lambda m: sum(1 for x in range(1, m + 1) if math.gcd(x, m) == 1)
    assert len(N) == p * q * math.factorial(k) * phi(p * q), (len(N), p * q * math.factorial(k) * phi(p * q))
    # H-orbits on pairs, with a canonical orientation per orbit (first pair u<v oriented u->v, propagated by h)
    pairs = [(u, v) for u in range(n) for v in range(u + 1, n)]
    orb_of = {}; orbits = []; orient = {}   # orient[(u,v) ordered] = +1 if the canonical orientation is u->v
    for pr in pairs:
        if pr in orb_of: continue
        idx = len(orbits); orb = []
        u, v = pr
        while (min(u, v), max(u, v)) not in orb_of:
            orb_of[(min(u, v), max(u, v))] = idx; orb.append((min(u, v), max(u, v)))
            orient[(u, v)] = 1; orient[(v, u)] = -1
            u, v = h[u], h[v]
        orbits.append(orb)
    c = len(orbits)
    assert c == (p + q) // 2 + 2 * k + k * (k - 1) // 2
    W_ALL = np.arange(1 << c, dtype=np.uint32)
    # per g: sigma_g (orbit permutation), eps_g (orientation flip per orbit), rho_g (reversal parity per orbit)
    gdata = []
    for g in N:
        sigma = [None] * c; eps = [0] * c; rho = [0] * c
        for i, orb in enumerate(orbits):
            u, v = orb[0]                      # canonical orientation u->v
            gu, gv = g[u], g[v]
            j = orb_of[(min(gu, gv), max(gu, gv))]
            sigma[i] = j
            eps[i] = 0 if orient[(gu, gv)] == 1 else 1   # image of canonical arc is canonical in orbit j?
            rho[i] = sum(1 for (a, b) in orb if g[a] > g[b]) & 1
        # consistency: sigma is a permutation and the orientation flip is well defined over the whole orbit
        assert sorted(sigma) == list(range(c))
        for i, orb in enumerate(orbits):
            for (a, b) in orb:
                ca, cb = (a, b) if orient[(a, b)] == 1 else (b, a)   # the canonical arc of this pair
                ga, gb = g[ca], g[cb]
                assert orb_of[(min(ga, gb), max(ga, gb))] == sigma[i]
                assert (0 if orient[(ga, gb)] == 1 else 1) == eps[i]   # flip status is constant along the orbit
        gdata.append((sigma, sum(1 << i for i in range(c) if eps[i]), sum(1 << i for i in range(c) if rho[i])))

    def act(words, sigma):
        out = np.zeros_like(words)
        for i in range(c):
            out |= ((words >> np.uint32(i)) & np.uint32(1)) << np.uint32(sigma[i])
        return out

    def popparity(x):
        x = x.copy()
        par = np.zeros_like(x)
        for i in range(c): par ^= (x >> np.uint32(i)) & np.uint32(1)
        return par

    # N-orbit reps (min over N) and N-evenness, for graphs and tournaments
    minG = W_ALL.copy(); minT = W_ALL.copy(); Nodd = np.zeros(1 << c, dtype=bool)
    for (sigma, epsmask, rhomask) in gdata:
        gw = act(W_ALL, sigma)
        minG = np.minimum(minG, gw)
        fixed = (gw == W_ALL)
        Nodd |= fixed & (popparity(W_ALL & np.uint32(rhomask)) == 1)
        gwT = gw ^ np.uint32(sum(1 << sigma[i] for i in range(c) if (epsmask >> i) & 1))
        minT = np.minimum(minT, gwT)
    Neven = ~Nodd
    # N-evenness must be constant on N-orbits
    repsG = {}
    for w, r in zip(W_ALL.tolist(), minG.tolist()):
        repsG.setdefault(r, set()).add(bool(Neven[w]))
    assert all(len(s) == 1 for s in repsG.values()), "N-evenness not constant on an N-orbit"
    R_T = len(set(minT.tolist())); R_G = len(repsG); RG_even = sum(1 for r in repsG if Neven[r])
    print("(p,q,k)=(%d,%d,%d) n=%d c=%d |N|=%d: N-orbits on tournaments R=%d, on graphs %d, N-even graph orbits=%d  [%.1fs]" % (
        p, q, k, n, c, len(N), R_T, R_G, RG_even, time.time() - t0), flush=True)

    # types and signs
    certT = {}; certG = {}; evenG = {}
    orbit_masks = [sum(1 << 0 for _ in ()) for _ in orbits]
    for w in range(1 << c):
        adjT = {v: [] for v in range(n)}; adjG = {v: [] for v in range(n)}; edges = []
        for i, orb in enumerate(orbits):
            bit = (w >> i) & 1
            for (a, b) in orb:
                # tournament: canonical orientation is a->b iff orient[(a,b)] == 1 ; bit 1 reverses it
                src, dst = (a, b) if (orient[(a, b)] == 1) != bool(bit) else (b, a)
                adjT[src].append(dst)
                if bit:
                    adjG[a].append(b); adjG[b].append(a); edges.append((a, b))
        certT[w] = pynauty.certificate(pynauty.Graph(n, directed=True, adjacency_dict=adjT))
        g = pynauty.Graph(n, directed=False, adjacency_dict=adjG)
        certG[w] = pynauty.certificate(g)
        gens_ = pynauty.autgrp(g)[0]
        evenG[w] = all((sum(1 for (a, b) in edges if gg[a] > gg[b]) & 1) == 0 for gg in gens_)
    print("  certificates done [%.1fs]" % (time.time() - t0), flush=True)
    # tournaments: m_T = number of N-orbits per type
    from collections import defaultdict
    orbT = defaultdict(set)
    for w in range(1 << c): orbT[certT[w]].add(int(minT[w]))
    tau_T = len(orbT); E_T = sum(len(s) - 1 for s in orbT.values())
    orbG = defaultdict(set); orbGev = defaultdict(set); signtype = {}
    for w in range(1 << c):
        ct = certG[w]; r = int(minG[w])
        orbG[ct].add(r)
        if Neven[w]: orbGev[ct].add(r)
        signtype.setdefault(ct, set()).add(evenG[w])
    assert all(len(s) == 1 for s in signtype.values()), "sign-evenness not constant on a type"
    sign_even = {ct: next(iter(s)) for ct, s in signtype.items()}
    # sign-even => N-even
    assert all(Neven[w] for w in range(1 << c) if evenG[w])
    tau_E = sum(1 for ct in orbG if sign_even[ct]); E_even = sum(len(orbG[ct]) - 1 for ct in orbG if sign_even[ct])
    witnesses = [ct for ct in orbG if not sign_even[ct] and len(orbGev[ct]) >= 1]
    W = len(witnesses); E_wit = sum(len(orbGev[ct]) - 1 for ct in witnesses)
    defect = tau_T - tau_E; ident = W + E_even + E_wit - E_T
    print("  tau_T=%d tau_E=%d defect=%d W=%d E_T=%d E_even=%d E_wit=%d identity %s; Theorem R (R = N-even graph orbits): %s  [%.1fs]" % (
        tau_T, tau_E, defect, W, E_T, E_even, E_wit, "OK" if defect == ident else "FAILS", "OK" if R_T == RG_even else "FAILS", time.time() - t0), flush=True)
    return dict(p=p, q=q, k=k, n=n, c=c, N=len(N), R=R_T, RG=RG_even, tau_T=tau_T, tau_E=tau_E, defect=defect, W=W, E_T=E_T, E_even=E_even, E_wit=E_wit)

if __name__ == "__main__":
    p, q, k = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    main(p, q, k)
