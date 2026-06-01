#!/usr/bin/env python3
"""
lrc_parity_goldbach_lemoine_s521.py   claudebox-2026-06-01-S521

Three parities and the doubled prime (reflection:
07-reflections/lrc-parity-cycles-numbers-doubled-primes-s521.md).

Even n = p + q (Goldbach, two prime atoms); Odd n = p + 2q (Lemoine, one prime +
one DOUBLED prime). The doubled prime 2q is the parity carrier (odd primes sum even;
the x2 supplies the odd adjustment). The factor 2 is the odd<->even bridge across
cycles (2-cover / OCF odd cycles), numbers (2-adic, first-even bridge), and LRC
speeds (the BENIGN doubling 2v vs the OBSTRUCTING additive sum-relation v_i+v_j=v_k).
"""
def sieve(N):
    s = bytearray([1])*(N+1); s[0] = s[1] = 0
    for i in range(2, int(N**.5)+1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return s
N = 4000; s = sieve(N); primes = [i for i in range(2, N) if s[i]]; pset = set(primes)
def goldbach(n): return [(p, n-p) for p in primes if p <= n-p and (n-p) in pset]
def lemoine(n):  return [(n-2*q, q) for q in primes if 2*q < n and (n-2*q) in pset]

def main():
    print("Even n = p+q (Goldbach); Odd n = p + 2q (Lemoine; 2q = DOUBLED prime, the parity carrier):")
    for n in range(6, 21):
        if n % 2 == 0:
            d = goldbach(n); print(f"  n={n:2} EVEN: {len(d)} reps, e.g. {d[0][0]}+{d[0][1]}")
        else:
            d = lemoine(n); print(f"  n={n:2} ODD : {len(d)} reps, e.g. {d[0][0]}+2*{d[0][1]}")
    print("\n  Odd primes sum to an EVEN total; reaching an ODD total needs one EVEN adjustment")
    print("  = exactly one DOUBLED prime 2q (Lemoine). The factor 2 is the odd<->even bridge.")
    print("\n  LRC binding pair v_i+v_j=n (Thm B, tight): n even => same-parity pair (Goldbach-type);")
    print("  n odd => (odd,even) pair, the even member a doubled atom. First-even bridge n=2*odd")
    print("  = one doubling of an odd n. Doubling 2v is BENIGN (2-adic edge); the additive prime")
    print("  PAIR (sum-relation v_i+v_j=v_k) is the obstruction. Doubled primes carry parity, not difficulty.")

if __name__ == "__main__":
    main()
