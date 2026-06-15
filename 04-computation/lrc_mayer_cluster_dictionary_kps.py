"""
lrc_mayer_cluster_dictionary_kps.py

Verifies the Mayer/Ursell cluster-expansion reading of the LRC(14) singular series
L(S) = (6/7)^13 + sum_{exact relations T} (6/7)^{13-|T|} (-1)^{|T|} prod_{v in T} s(t_v),
s(t) = sin(pi*t/7)/(pi*t).

Goals:
 1. Reproduce the ideal/free term (6/7)^13 and confirm L as ideal * [1 + signed cluster corrections].
 2. Confirm the 7-vanishing s(7k)=0 acts as a cluster WEIGHT ZERO (gate of the apex prime 7).
 3. Confirm |T|=2 (pairwise) is ABSOLUTELY convergent and matches the bound |P|<=g^2/(3 v_a v_b).
 4. Confirm |T|=3 SIGNED sum stabilizes (conditional convergence) while ABSOLUTE sum grows (THM-504 B).

This is a STRUCTURAL/numerical sanity pass for the cluster dictionary, not a new theorem.
"""
import math, itertools

def s(t):
    if t == 0:
        return 0.0  # not used; relations have t_v != 0
    return math.sin(math.pi * t / 7.0) / (math.pi * t)

IDEAL = (6.0/7.0)**13

# ---------- 1. 7-vanishing as cluster weight zero ----------
print("=== 1. 7-vanishing (apex-prime-7 gate = cluster weight zero) ===")
for t in [1,2,3,6,7,8,13,14,21]:
    print(f"  s({t:2d}) = {s(t):+.6f}   {'<-- WEIGHT ZERO (7|t)' if t%7==0 else ''}")
print(f"  ideal/free term (6/7)^13 = {IDEAL:.6f}\n")

# ---------- 2. pairwise (|T|=2) lattice factor P(a',b') and the absolute bound ----------
def P_pair(va, vb, K=200000):
    g = math.gcd(va, vb)
    ap, bp = vb//g, va//g           # coprime
    tot = 0.0; abstot = 0.0
    for k in range(1, K+1):
        term = s(k*ap)*s(k*bp)
        tot += term; abstot += abs(term)
    P = 2*tot                        # s odd => product even in k
    bound = g*g/(3.0*va*vb)
    return P, 2*abstot, bound, ap, bp

print("=== 2. |T|=2 pairwise cluster: ABSOLUTE convergence + bound g^2/(3 va vb) ===")
for (va,vb) in [(14,21),(7,14),(98,105),(14,17),(2,3),(3,5)]:
    P, Pabs, bnd, ap, bp = P_pair(va,vb)
    flag = "(7|a' or 7|b' => P=0)" if (ap%7==0 or bp%7==0) else ""
    print(f"  {{{va},{vb}}}: P={P:+.6f}  |P|sum={Pabs:.6f}  bound={bnd:.6f}  ok={abs(P)<=bnd+1e-9} {flag}")
print()

# ---------- 3. |T|=3 signed vs absolute (conditional convergence, THM-504 B) ----------
# Sigma_3 over 3-term exact relations t1 v1 + t2 v2 + t3 v3 = 0 on a generic-coprime triple.
# We enumerate by coefficient cutoff Tmax and report signed and absolute partial sums.
def sigma3_triple(v, Tmax):
    v1,v2,v3 = v
    signed = 0.0; absol = 0.0; count = 0
    for t1 in range(-Tmax, Tmax+1):
        if t1 == 0: continue
        for t2 in range(-Tmax, Tmax+1):
            if t2 == 0: continue
            num = -(t1*v1 + t2*v2)
            if num % v3 != 0:
                continue
            t3 = num // v3
            if t3 == 0 or abs(t3) > Tmax:
                continue
            w = s(t1)*s(t2)*s(t3)
            signed += w; absol += abs(w); count += 1
    return signed, absol, count

print("=== 3. |T|=3 cluster: SIGNED stabilizes, ABSOLUTE grows (conditional convergence) ===")
triple = (3,4,5)   # generic-coprime; reproduces THM-504's qualitative behavior
for Tmax in [13,20,27,40]:
    sg, ab, c = sigma3_triple(triple, Tmax)
    ratio = abs(sg)/ab if ab>0 else 0
    print(f"  triple{triple} Tmax={Tmax:3d}: signed Sigma3={sg:+.5f}  abs A3={ab:.5f}  |sig|/A3={ratio:.3f}  (N3={c})")
print()
print(f"NOTE: signed stable + absolute growing + ratio falling = conditional-convergence signature (THM-504 B).")
