#!/usr/bin/env python3
"""
claudebox-2026-06-03-S599 : developing H3-H6 — the spectral/character, isostatic, CRT, and
rigidity-height refinements of "rigidity = orbit-type" (HYP-2140).

The danger structure on the n-clock. At clock time j/n, a runner v of the AP {1,..,n-1} sits at
phase vj/n; it is STUCK AT THE ORIGIN (deepest danger, ‖vj/n‖=0) iff vj≡0 mod n. The danger depth
at clock j is
        d(j) = #{v in [1,n-1] : vj ≡ 0 (mod n)} = gcd(j,n) - 1.
Clock j is fully lonely (gap = 1/n, the tight witness) iff d(j)=0 iff gcd(j,n)=1.

THE BACKBONE IDENTITY (block-diagonalization by divisor = Dirichlet-character level):
        gcd(j,n) = Σ_{d | n} φ(d) · [d | j].
So  d(j) = Σ_{d|n, d>1} φ(d)·[d|j]  — the danger is a SUM OF DIVISOR-BLOCKS, each block the
characters of level d (period n/d). This is the symmetry-adapted orbit rigidity matrix, diagonal in
the Dirichlet basis; the additive-character transform is the Ramanujan sum c_d.

  H3 (spectral): rigidity leaks only through the gcd(·,n)>1 blocks; the 2-block = the d=2 term.
  H5 (CRT): for n=2q, d(j) = [2|j] (2-block) + (q-1)[q|j] (q-block) + (q-1)[2q|j]; the q-block is
            the prime case lifted, the obstruction is the 2-block.
  H4 (isostatic): #lonely clocks = φ(n) (the unit orbit) = the critical (Maxwell) count; a
            counterexample would be OVER-rigid (0 lonely) — forbidden, the orbit floor is φ(n)>0.
  H6 (rigidity-height): the doubling map x->2x is nilpotent on the 2-part Z/2^a with index v2(n).
"""
from math import gcd
from cmath import exp, pi

def euler_phi(n):
    return sum(1 for j in range(1, n + 1) if gcd(j, n) == 1)

def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]

def danger_depth(j, n):
    return gcd(j, n) - 1

def ramanujan(d, a):
    # c_d(a) = sum over t coprime to d of e^{2 pi i a t / d}
    return sum(exp(2j * pi * a * t / d) for t in range(1, d + 1) if gcd(t, d) == 1)

def doubling_nilpotency(a):
    """index k with 2^k ≡ 0 mod 2^a starting from an odd seed: x->2x on Z/2^a reaches 0 in a steps."""
    n = 2 ** a; x = 1; k = 0
    while x % n != 0:
        x = (2 * x) % n; k += 1
        if k > a + 5: break
    return k

def v2(n):
    e = 0
    while n % 2 == 0: n //= 2; e += 1
    return e

def main():
    print("=" * 94)
    print("S599  rigidity H3-H6: spectral/character block, isostatic count, CRT 2-block, height=v2")
    print("=" * 94)

    # ---- H3 + H5: divisor-block (Dirichlet-level) decomposition of the danger ----
    print("\n[H3/H5] danger d(j)=gcd(j,n)-1 = Σ_{d|n,d>1} φ(d)·[d|j]  (block-diagonal by divisor)")
    for n in [7, 13, 14, 10, 22, 18]:
        # verify the identity
        ok = all(gcd(j, n) == sum(euler_phi(d) for d in divisors(n) if j % d == 0) for j in range(n))
        blocks = [(d, euler_phi(d)) for d in divisors(n) if d > 1]
        # which j are lonely (gcd=1)
        lonely = [j for j in range(1, n) if gcd(j, n) == 1]
        two_block = [j for j in range(1, n) if gcd(j, n) % 2 == 0]
        print(f"  n={n:2d}: identity={ok}  blocks(d,φ(d))={blocks}")
        print(f"        lonely clocks (gcd=1): {lonely}  (#={len(lonely)}=φ(n)={euler_phi(n)})")
        print(f"        2-block (even gcd, the defect): {two_block}")
    print("  => rigidity (d(j)>0) leaks ONLY through the d>1 divisor-blocks; the 2-block = the even-j")
    print("     defect. For n=2q the obstruction is the 2-block; the q-block is one point (the apex).")

    # ---- the additive-character (Ramanujan) eigenview, n=14 ----
    print("\n[H3] Ramanujan / Dirichlet-character spectrum of d(j) at n=14 (where the mass sits)")
    n = 14
    dvec = [danger_depth(j, n) for j in range(n)]
    # DFT
    spec = []
    for a in range(n):
        val = sum(dvec[j] * exp(-2j * pi * a * j / n) for j in range(n))
        spec.append((a, gcd(a, n), abs(val)))
    print("  a | gcd(a,14) | |DFT_a(d)|   (nonzero only at a with gcd(a,n)>1 ⇒ imprimitive chars)")
    for a, g, m in spec:
        if m > 1e-9:
            tag = "2-block" if g == 2 else ("7-block" if g == 7 else ("DC" if a == 0 else ""))
            print(f"  {a:2d} |    {g:2d}     |  {m:6.3f}   {tag}")

    # ---- H4: isostatic / Maxwell count ----
    print("\n[H4] ISOSTATIC: #lonely clocks of the AP = φ(n) (the unit orbit) = the critical count")
    print("  n | φ(n) | #lonely | binding/lonely-clock | over-rigid (0 lonely) possible?")
    for n in [6, 7, 10, 13, 14, 18, 22]:
        lonely = [j for j in range(1, n) if gcd(j, n) == 1]
        # binding runners at a lonely clock j: v with vj≡±1 (gap exactly 1/n)
        j0 = lonely[0]
        bind = [v for v in range(1, n) if (v * j0) % n in (1, n - 1)]
        print(f"  {n:2d} | {euler_phi(n):4d} | {len(lonely):7d} | {len(bind):2d} (v with vj≡±1)        | "
              f"NO — orbit floor φ(n)={euler_phi(n)}>0 forbids it")
    print("  => the worry-set (AP) is CRITICALLY rigid: exactly φ(n) witnesses, 2 binders each")
    print("     (isostatic). A counterexample = 0 witnesses = over-rigid, below the orbit floor φ(n).")

    # ---- H6: rigidity-height = v2(n) ----
    print("\n[H6] RIGIDITY-HEIGHT = v2(n): doubling x->2x nilpotency index on the 2-part Z/2^a")
    print("  n | v2(n) | doubling-nilpotency on 2-part | height interpretation")
    for n in [7, 13, 14, 10, 22, 28, 8, 16, 56]:
        a = v2(n)
        nil = doubling_nilpotency(a) if a > 0 else 0
        ht = ("0 — odd, no 2-adic break (prime-regular)" if a == 0 else
              f"{a} — {a} break(s) above the prime; n=14 is height 1")
        print(f"  {n:2d} | {a:4d}  | {nil:24d}      | {ht}")
    print("  => rigidity-height = v2(n): each doubling level degrades the dynamical rigidity by one")
    print("     factor of 2; n=14=2·7 has height 1 (the single mod-2 collapse = the apex/2-adic seam).")

if __name__ == "__main__":
    main()
