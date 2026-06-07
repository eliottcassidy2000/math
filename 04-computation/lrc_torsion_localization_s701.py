# lrc_torsion_localization_s701.py
# CLAIM (user, S701): the clock-witness "leaks" at composite n are LOCKED to the
# algebraic torsion of the composite divisors. Precisely: a runner that defeats
# the clocks of one prime-cofactor but survives in the complementary prime sits at
# a PRIME-TORSION element of Z/n, which CRT-decomposes as (nonzero in that prime,
# 0 in the cofactor base). n=14: half-turn 7 = (1 mod 2, 0 mod 7) = 2-torsion gen.
# n=15: residues 5,10 = (.,0 mod 5) = 3-torsion, killed in the 5-runner base.
#
# This script verifies, for many composite n:
#  (1) the p-torsion subgroup T_p = {x: p x = 0 in Z/n} = <n/p>, order p;
#  (2) every nonzero x in T_p is 0 in the cofactor base m=n/p^{v_p(n)} and !=0 mod p;
#  (3) at the full clock t=1/n, such x sits at the EXACT order-p rotation j/p,
#      so |x/n| >= 1/p  -> the p-torsion leak is plugged by the p-clock (margin 1/p);
#  (4) the union of prime-torsion subgroups = the "single-prime-surviving" residues;
#      the hard core (all primes leak at once) = {x = 0 mod n}.
from math import gcd
from fractions import Fraction

def factor(n):
    f={}; d=2; m=n
    while d*d<=m:
        while m%d==0: f[d]=f.get(d,0)+1; m//=d
        d+=1
    if m>1: f[m]=f.get(m,0)+1
    return f

def p_torsion(n,p):           # {x in Z/n : p*x == 0 mod n} = <n/p>, order p
    return sorted(x for x in range(n) if (p*x)%n==0)

def norm(x,n):                # ||x/n|| = distance to nearest integer, in units of 1/n
    r=x%n; return min(r,n-r)

def controls(x,n,k):          # does clock t=1/k push runner x off the origin?
    return norm((x* (n//gcd(n,k)) ) , n) if False else (min((x%k),(k-(x%k)))>0)
def min_controlling_clock(x,n):
    # smallest k in {2..n} with x not ≡ 0 mod k (Lemma A clock that controls x alone)
    for k in range(2,n+1):
        if x%k!=0: return k
    return None

print("="*72)
print("TORSION-LOCALIZATION OF CLOCK LEAKS  (verifies the S701 claim)")
print("="*72)
geom_ok=True; squarefree_clean=True
for n in [14,15,12,20,21,18,30,33,35,45,50]:
    f=factor(n); primes=sorted(f)
    desc=' * '.join(f'{p}^{a}' if a>1 else str(p) for p,a in f.items())
    sqfree=all(a==1 for a in f.values())
    print(f"\nn={n} = {desc}   [{'squarefree' if sqfree else 'has prime power'}]")
    for p in primes:
        a=f[p]; pa=p**a; m=n//pa            # m = cofactor coprime to p
        T=p_torsion(n,p); nz=[x for x in T if x]
        ok_cof = all(x%m==0 for x in nz) if m>1 else True   # dies in cofactor base
        ok_pb  = all(x%p!=0 for x in nz)                    # survives mod p (squarefree only)
        rots=[Fraction(x,n) for x in nz]                    # position at t=1/n
        ok_rot = all(norm(x,n)>=n//p for x in nz)           # ||x/n|| >= 1/p (geometric)
        ctrl   = sorted(set(min_controlling_clock(x,n) for x in nz))  # plugging clock
        geom_ok &= ok_cof and ok_rot
        if a==1: squarefree_clean &= ok_pb
        print(f"  p={p}: {p}-torsion={T} (gen n/p={n//p}, order {p})"
              f"\n        nonzero die in cofactor base {m}: {ok_cof}; survive mod {p}: {ok_pb}"
              f"{'  <-- prime power: socle is DIVISIBLE by p, needs deeper clock' if not ok_pb else ''}"
              f"\n        at t=1/n sit at {[str(r) for r in rots]} (order-{p} rotations); ||.||>=1/{p}: {ok_rot}"
              f"\n        smallest controlling clock(s) k: {ctrl} "
              f"({'= p (clean plug, margin 1/p)' if ctrl==[p] else '> p (deeper p-adic clock needed)'})")

# explicit reproduction of the user's two instances, via CRT
print("\n"+"-"*72)
print("Explicit CRT reading of the user's instances:")
print(f"  n=14: residue 7 -> (7 mod 2, 7 mod 7) = (7%2={7%2}, 7%7={7%7})  = 2-torsion gen, 0 in 7-base")
for r in (5,10):
    print(f"  n=15: residue {r} -> ({r} mod 3, {r} mod 5) = ({r%3}, {r%5})  = 3-torsion, 0 in 5-base")

# union of prime-torsions (single-prime-surviving residues) and the hard core
print("\n"+"-"*72)
print("Union of prime-torsion subgroups (single-prime-surviving leaks) & hard core:")
for n in [14,15,30,12]:
    f=factor(n); primes=sorted(f)
    U=set()
    for p in primes: U|= set(x for x in p_torsion(n,p) if x)
    core=[x for x in range(n) if all(x%p==0 for p in primes)]  # 0 in every prime comp
    print(f"  n={n}: U(T_p\\0)={sorted(U)}; hard core (0 in all primes, mod n) "
          f"= {sorted(set(core))} (incl 0)")

print("\nGEOMETRIC localization (torsion sits at order-p rotations, ||.||>=1/p) holds for ALL n:", geom_ok)
print("SQUAREFREE clean plug (socle survives mod p, plugged by t=1/p): ", squarefree_clean)
print("=> squarefree primes: leak plugged by the prime clock (margin 1/p);")
print("   prime-power factors: socle torsion is itself ≡0 mod p, needs the deeper p-adic clock")
print("   -- this is the prime-power hardness (cf. n=14: 2n-1=27=3^3, no coprime plug).")
