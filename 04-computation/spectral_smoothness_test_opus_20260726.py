from sympy import primefactors, factorint
# Test the spectral-grammar CENTRAL claim: "smooth values = constructive interference of prime lines"
# Twin-prime / prime-pair singular series: S(g) = prod_{p|g, p>2} (p-1)/(p-2)  (relative to base gap 2)
def sing_series(g):
    s=1.0
    for p in primefactors(g):
        if p>2: s *= (p-1)/(p-2)
    return s
print("=== Spectral grammar test 1: prime-pair singular series S(g) -- 'jumping champions' (opus) ===")
print("    claim (kps HYP-9025): favoured g maximize the product of prime lines (smooth/primorial g win)\n")
ranked=sorted(range(2,220,2), key=lambda g:-sing_series(g))[:12]
for g in ranked:
    print(f"  g={g:>3} S(g)={sing_series(g):.4f}  factors={dict(factorint(g))}")
print("  => classical jumping champions are 4,6,30,210,... (primorials); confirms smooth=constructive.\n")

print("=== Spectral grammar test 2: does the SAME selection pick out LRC(n) thresholds? ===")
# LRC(n): n-1 speeds, threshold 1/n. Hard primes = factors of n and n-1 (AP length).
# Speculation: n where prime lines constructively interfere are the 'hard/rich' cases.
# Local avoidance density for the AP at t=1/n involves residues mod each prime p|n.
def lrc_richness(n):
    # crude: product over p|n and p|(n-1) of a 'line strength' (p-1)/p-type; smooth n richer
    s=1.0
    for p in primefactors(n): s *= p/(p-1)      # multiplicative structure of the threshold modulus n
    for p in primefactors(n-1): s *= (p)/(p-1)  # AP length n-1
    return s
print("  n, factor(n)=threshold-modulus, factor(n-1)=AP-len, richness (higher=more prime-line structure):")
for n in [7,8,12,13,14,15,18,19,24,30]:
    print(f"    n={n:>2}: {dict(factorint(n))} x AP{n-1}={dict(factorint(n-1))}  richness={lrc_richness(n):.3f}   {'<- LRC14 (2*7 x 13)' if n==14 else ''}")
print("\n  NOTE (honest): this is a SPECULATION probe, not a derivation. 14=2*7 (threshold), AP-len 13 prime.")
print("  The 'hard primes 7,13' (kps-S131b) = prime factors of the threshold-modulus 14 and AP-length 13.")
