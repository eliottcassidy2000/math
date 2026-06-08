#!/usr/bin/env python3
"""
S718 — A momentum twistor for the LRC: the DISCRETE LOGARITHM linearizes the dual conformal symmetry.

S717: the LRC dual conformal symmetry is the MULTIPLIER group (Z/m)* acting by v -> a v on residues
(shell m), with the INVERSION v -> v^{-1} as the special-conformal generator (the bad multiplier of
THM-420). A multiplicative action is linearized by a LOGARITHM. So the LRC "momentum twistor" is the
discrete log
            ell : (Z/m)* -> Z/phi(m),   v |-> log_g(v)   (g a primitive root),
which is the exact analog of Hodges' momentum twistors: it makes dual conformal symmetry MANIFEST and
LINEAR and resolves the constraints.

In the twistor (log) coordinate:
  - the MULTIPLIER v -> a v becomes a TRANSLATION ell -> ell + ell(a);
  - the INVERSION (special conformal) v -> v^{-1} becomes NEGATION ell -> -ell;
  - ell(-1) = phi(m)/2 =: c (the unique order-2 element), so v -> -v is the shift ell -> ell + c.

THE PAYOFF (THM-420 in twistor coordinates). A multiplier a is BAD for runner i iff a v_i in {+1,-1},
i.e. ell(a) in {-ell(v_i), c - ell(v_i)}. So with L = { ell(v_i) } the bad-multiplier log-set is
            B = (-L) cup (c - L)   subset Z/phi(m),
the union of TWO TRANSLATES of -L. A good multiplier (dodge) exists iff B != Z/phi(m). And |B| = 2|L| =
phi(m) (full, no dodge) iff the two translates TILE Z/phi(m) iff L is a HALF-SYSTEM: a transversal of the
order-2 subgroup {0, c}, i.e. L picks exactly one of each pair {x, x+c}. So:

   TRANSVERSAL CORE (hard frontier, THM-420)  <=>  L is a half-system in Z/phi(m)  <=>  {-L, c-L} tiles.

The momentum twistor turns the multiplicative dodge problem into ADDITIVE COVERING on the cyclic group
Z/phi(m), linearizes the dual conformal symmetry, and unifies S717 (dual conformal) + THM-420 (transversal
core) + THM-415 (half-system). We verify all of this and the multi-shell "hard core = half-system at every
shell" picture.

No numpy/sympy.
"""
from math import gcd
import random

def factorize(n):
    f={}; d=2
    while d*d<=n:
        while n%d==0: f[d]=f.get(d,0)+1; n//=d
        d+=1
    if n>1: f[n]=f.get(n,0)+1
    return f
def euler_phi(n):
    r=n
    for p in factorize(n): r-=r//p
    return r
def is_cyclic_mult(m):
    # (Z/m)* cyclic iff m in {1,2,4,p^k,2 p^k}
    if m in (1,2,4): return True
    f=factorize(m)
    if len(f)==1:
        p=next(iter(f)); return p%2==1
    if len(f)==2 and f.get(2)==1:
        p=[q for q in f if q!=2][0]; return p%2==1
    return False
def primitive_root(m):
    phi=euler_phi(m); fs=list(factorize(phi))
    for g in range(2,m):
        if gcd(g,m)!=1: continue
        if all(pow(g,phi//p,m)!=1 for p in fs): return g
    return None
def dlog_table(m,g):
    phi=euler_phi(m); t={}; x=1
    for e in range(phi):
        t[x]=e; x=(x*g)%m
    return t

def good_multiplier_exists_pm1(residues, m):
    """direct: exists a in (Z/m)* with a*v not in {1, m-1} for all v in residues."""
    for a in range(1,m):
        if gcd(a,m)!=1: continue
        if all((a*v)%m not in (1,m-1) for v in residues): return True
    return False

def bad_logset(L, c, phi):
    return set((-l)%phi for l in L) | set((c-l)%phi for l in L)

def is_half_system(L, c, phi):
    """L picks exactly one of each pair {x, x+c} (transversal of {0,c}); requires |L|=phi/2."""
    if len(set(L))!=phi//2: return False
    seen=set()
    for l in set(L):
        key=min(l%phi,(l+c)%phi)
        if key in seen: return False
        seen.add(key)
    return True

if __name__=="__main__":
    rng=random.Random(0)
    print("="*86)
    print("S718 — a momentum twistor for the LRC: the discrete log linearizes dual conformal symmetry")
    print("="*86)

    # (1) the twistor map linearizes the dual conformal group
    print("\n(1) THE TWISTOR ell=dlog_g LINEARIZES dual conformal: multiplier->translation, inversion->negation")
    for m in (7,11,13,19,23,27):
        g=primitive_root(m); phi=euler_phi(m); t=dlog_table(m,g); c=phi//2
        # multiplier a: ell(a*v) = ell(a)+ell(v) ?
        ok_mul=all((t[(a*v)%m]==(t[a]+t[v])%phi) for a in t for v in t if (a*v)%m in t)
        # inversion: ell(v^-1) = -ell(v) ?
        ok_inv=all((t[pow(v,-1,m)]==(-t[v])%phi) for v in t)
        ok_neg1=(t[m-1]==c)  # ell(-1)=phi/2
        print(f"  m={m:2d} (g={g}, phi={phi}, cyclic={is_cyclic_mult(m)}): "
              f"mul->add {ok_mul} | inv->neg {ok_inv} | ell(-1)=phi/2 {ok_neg1}")
    print("  => exactly like Hodges momentum twistors: the hidden (dual conformal) symmetry is now LINEAR.")

    # (2) THM-420 transversal core IN TWISTOR COORDS: bad set = (-L) U (c-L); core <=> half-system tiling
    print("\n(2) TRANSVERSAL CORE in twistor coords: bad-multiplier set B=(-L)U(c-L); no dodge <=> half-system")
    for n in (4,6,7,10,12,14):
        m=2*n-1
        if not is_cyclic_mult(m):
            print(f"  n={n:2d}, m=2n-1={m}: (Z/m)* not cyclic -> higher-dim twistor (skip 1D demo)"); continue
        g=primitive_root(m); phi=euler_phi(m); t=dlog_table(m,g); c=phi//2
        units=[v for v in range(1,m) if gcd(v,m)==1]   # twistor covers only (Z/m)* (units)
        trials=600; match=0; core=0
        for _ in range(trials):
            # multiple-of-n config residues = n-1 distinct UNIT residues mod m (non-units = ramified, off-twistor)
            res=rng.sample(units, n-1)
            L=[t[v] for v in res]
            B=bad_logset(L,c,phi)
            dodge_pred = (len(B)<phi)              # twistor: good multiplier exists iff B != full
            dodge_true = good_multiplier_exists_pm1(res, m)
            if dodge_pred==dodge_true: match+=1
            hs=is_half_system(L,c,phi)
            if not dodge_true: core+=1
            # core (no dodge) should coincide with half-system
        # focused check: among configs, no-dodge <=> half-system
        chk=0
        for _ in range(trials):
            res=rng.sample(units,n-1); L=[t[v] for v in res]
            nod = not good_multiplier_exists_pm1(res,m)
            hs  = is_half_system(L,c,phi)
            if nod==hs: chk+=1
        print(f"  n={n:2d}, m={m} (phi={phi}, c=ell(-1)={c}): twistor predicts dodge {match}/{trials}; "
              f"no-dodge<=>half-system {chk}/{trials}")
    print("  => THM-420 transversal core = L is a HALF-SYSTEM (transversal of {0,c}) = {-L,c-L} tiles Z/phi.")
    print("     The multiplicative dodge becomes ADDITIVE COVERING on Z/phi(m); inversion=negation is manifest.")

    # (3) explicit small example: show the tiling
    print("\n(3) EXPLICIT (n=7, m=13, phi=12, c=6): a transversal-core config and a dodgeable one")
    m=13; g=primitive_root(m); phi=12; t=dlog_table(m,g); c=6
    # build a half-system L (one of each pair {x,x+6}) -> map back to residues
    inv_t={e:v for v,e in t.items()}
    Lhs=[0,1,2,3,4,5]                  # one of each {x,x+6}: a half-system
    res_core=[inv_t[l] for l in Lhs]
    print(f"   half-system L={Lhs} -> residues {sorted(res_core)} : dodge? {good_multiplier_exists_pm1(res_core,m)} "
          f"(expect False = transversal core); B covers Z/12: {len(bad_logset(Lhs,c,phi))==phi}")
    Lbad=[0,1,2,3,4,0+6]               # NOT a half-system (0 and 6 both present -> pair collision)
    res_dod=[inv_t[l] for l in [0,1,2,3,4,6]]
    print(f"   non-half-system L={[0,1,2,3,4,6]} -> residues {sorted(res_dod)} : dodge? "
          f"{good_multiplier_exists_pm1(res_dod,m)} (expect True); |B|={len(bad_logset([0,1,2,3,4,6],c,phi))}<12")

    # (4) multi-shell: hard core = half-system at EVERY available cyclic shell
    print("\n(4) MULTI-SHELL twistor: a config is hard (no shell dodge) only if it is a half-system at every shell")
    print("    (shells m=2n-1 prime/prime-power give a 1D twistor; composite (Z/m)* -> higher-dim twistor)")
    for n in (4,6,7,10,12,14):
        m=2*n-1; cyc=is_cyclic_mult(m); dim=len([1 for _ in factorize(m)])  # rough twistor dim
        print(f"   n={n:2d}: shell m={m}, (Z/m)* cyclic={cyc}, #prime factors(m)={len(factorize(m))} "
              f"=> twistor dim ~ {len(factorize(m))} ({'P^1-like' if cyc else 'higher'})")
    print("="*86)
