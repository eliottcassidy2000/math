# lrc_antipodal_shellpartner_torsion_s702.py
# POKE STEERING TASK 1: connect the signed-floor binding shell-partner q=a+b (HYP-2296)
# to the prime-torsion / fiber picture (THM-421/427/428).
#
# THESIS (S702): the antipodal involution  sigma: x |-> -x  on Z/q unifies all three moduli.
#  - A shell-partner {a,b} (a+b ≡ 0 mod q) is a sigma-ORBIT on Z/q; synchronization
#    (THM-425 L0: ||a k/q||=||b k/q||) is exactly the sigma-INVARIANCE of ||.|| (||-x||=||x||).
#  - The SELF-partners (a ≡ b) are the sigma-FIXED POINTS = 2-torsion T_2(Z/q)={0, q/2};
#    a nonzero one exists iff q is EVEN, and it is the half-turn q/2.
#  - FACE SPLIT: the shell modulus 2n-1 is ALWAYS ODD  => sigma has NO nonzero fixed point
#    => the shell is "antipodally free" (every residue in a genuine 2-orbit). The clock
#    modulus n can be EVEN => the half-turn n/2 is a sigma-fixed 2-torsion leak (poke's n=14
#    residue 7 = 14/2). So the prime 2 distinguishes the two faces antipodally.
#  - sigma commutes with CRT: a shell-partner decomposes into per-prime sigma-orbits
#    (the "fiber alignment" poke asks for).
#
# This script: (1) exact Gstar(S) minimizers for small n; (2) sigma/torsion analysis of every
# binding shell-partner q=a+b: orbit vs fixed point, parity of q, factorization, per-prime
# CRT split, gcd(q,n), gcd(q,2n-1); (3) the clock half-turn vs shell antipodal-freeness table.
from math import gcd
from fractions import Fraction
from itertools import combinations
from functools import reduce

def ndist(x,q):                      # ||x/q|| in units of 1/q : min(x mod q, q-(x mod q))
    r=x%q; return min(r,q-r)

def M_and_witness(W):
    # exact maximin of min_{w in W} ||w t||; optimum denom = a pair-sum <= 2*maxW
    W=sorted(set(W));
    if not W: return Fraction(1,2), None, []
    Q=2*max(W)
    best=Fraction(-1); bq=bk=None
    for q in range(2,Q+1):
        for k in range(1,q):
            g=min(ndist(w*k,q) for w in W)
            val=Fraction(g,q)
            if val>best:
                best=val; bq,bk=q,k
    binding=[w for w in W if Fraction(ndist(w*bk,bq),bq)==best]
    return best,(bq,bk),binding

def Gstar(S):
    S=list(S); r=len(S); best=Fraction(-1); info=None
    # cuts: A ranges over subsets containing S[0]  (dedup complementary cuts)
    for mask in range(1<<(r-1)):
        A=[S[0]]+[S[1+i] for i in range(r-1) if (mask>>i)&1]
        B=[v for v in S if v not in A]
        W=[]
        for x,y in combinations(A,2): W.append(abs(x-y))
        for x,y in combinations(B,2): W.append(abs(x-y))
        for x in A:
            for y in B: W.append(x+y)
        W=[w for w in W if w>0]
        if not W: continue
        m,wit,binding=M_and_witness(W)
        if m>best:
            best=m; info=(A,B,W,wit,binding)
    return best,info

def factor(n):
    f={};d=2;m=n
    while d*d<=m:
        while m%d==0: f[d]=f.get(d,0)+1;m//=d
        d+=1
    if m>1: f[m]=f.get(m,0)+1
    return f

def analyze_binding(a,b,q,n):
    # a+b should ≡ 0 mod q (shell-partner). sigma-orbit analysis.
    sigma_ok = (a+b)%q==0
    selfp = (a-b)%q==0                      # a ≡ b: fixed point (self-partner)
    fixed = [x for x in range(q) if (2*x)%q==0]   # 2-torsion of Z/q = sigma fixed set
    f=factor(q)
    # per-prime CRT split of a (mod each p^e) and whether the pair is self-partner in that fiber
    perprime={}
    for p,e in f.items():
        pe=p**e
        perprime[pe]=((a%pe),(b%pe), (a+b)%pe==0, (a-b)%pe==0)
    return dict(sigma_orbit=sigma_ok, self_partner=selfp, q_parity=('even' if q%2==0 else 'odd'),
                fixedpts=fixed, factor_q=f, gcd_q_n=gcd(q,n), gcd_q_2n1=gcd(q,2*n-1),
                perprime=perprime)

print("="*78)
print("POKE TASK 1 — antipodal (sigma) structure of the signed-floor shell-partner q=a+b")
print("="*78)

# (1) exact minimizers by small search + the HYP-2296 published minimizers as seeds
searches=[(5,14),(6,11),(7,10)]
found={}
for n,B in searches:
    r=n-1; best=Fraction(2); bestS=[]
    for S in combinations(range(1,B+1),r):
        if reduce(gcd,S)!=1: continue
        g,info=Gstar(S)
        if g<best:
            best=g; bestS=[(S,info)]
        elif g==best:
            bestS.append((S,info))
    found[n]=(best,bestS)
    print(f"\n--- n={n} (B<={B}): inf Gstar = {best} = {float(best):.4f}, n*inf={float(n*best):.3f}, "
          f"{len(bestS)} minimizer(s) ---")
    for S,info in bestS[:4]:
        A,Bs,W,wit,binding=info
        q,k=wit
        # find the shell-partner pair within binding: two binding rel-speeds summing to q
        pair=None
        for x in binding:
            for y in binding:
                if x<=y and (x+y)==q: pair=(x,y)
        print(f"   S={S}  cut A={A}|B={Bs}  t*={k}/{q}  M={W and ''}{found[n][0]}  binding={binding}  pair={pair}")
        if pair:
            an=analyze_binding(pair[0],pair[1],q,n)
            print(f"      q={q} [{an['q_parity']}]  factor={an['factor_q']}  sigma-orbit:{an['sigma_orbit']} "
                  f"self-partner:{an['self_partner']}  2-torsion(fixed)={an['fixedpts']}")
            print(f"      gcd(q,n)={an['gcd_q_n']}  gcd(q,2n-1)={an['gcd_q_2n1']}  per-prime(a,b,sum0?,self?)={an['perprime']}")

# (2) the HYP-2296 published minimizers, analyzed directly
print("\n"+"="*78)
print("Published HYP-2296 minimizers (mover-only & observer-inclusive) — sigma/torsion table")
print("="*78)
published=[
  ("mover (2,3,4,6,8)",   9,10, 6),   # 3/19  q=19
  ("mover (4,5,8,10,15)",  4,23, 6),  # 4/27  q=27
  ("mover (2,4,7,10,11,12)",19,23, 7),# 5/42  q=42
  ("obs   (5,6,7,8,9)",    6,13, 6),  # 2/19  q=19
  ("obs   (5,7,8,9,11)",  11,18, 6),  # 3/29  q=29
]
for tag,a,b,n in published:
    q=a+b; an=analyze_binding(a,b,q,n)
    print(f"\n {tag}: binding {{{a},{b}}}  q={q} [{an['q_parity']}] factor={an['factor_q']}")
    print(f"    sigma-orbit:{an['sigma_orbit']}  self-partner:{an['self_partner']}  2-torsion={an['fixedpts']}")
    print(f"    gcd(q,n)={an['gcd_q_n']} gcd(q,2n-1)={an['gcd_q_2n1']}  per-prime(a%pe,b%pe,sum≡0,self):{an['perprime']}")

# (3) the FACE SPLIT table: clock n (half-turn n/2 when even) vs shell 2n-1 (always odd, no fixed pt)
print("\n"+"="*78)
print("FACE SPLIT: clock Z/n half-turn (2-torsion) vs shell Z/(2n-1) antipodal-freeness")
print("="*78)
print(" n | n even? | clock half-turn n/2 (sigma-fixed 2-torsion) | shell 2n-1 | (2n-1) odd => no fixed pt | block q=4n-5 [odd]")
for n in range(5,16):
    ht = n//2 if n%2==0 else None
    print(f" {n:2d} | {str(n%2==0):5s}   | {('r='+str(ht)+'  (= n/2, the half-turn leak)') if ht else 'none (n odd: only 0 fixed)':40s} "
          f"| {2*n-1:3d}      | {str((2*n-1)%2==1):4s} (always)            | {4*n-5}")
print("\n=> shell 2n-1 and block q=4n-5 are ALWAYS ODD: antipodally free (genuine 2-orbits, no self-partner).")
print("   The clock n carries the half-turn self-partner n/2 exactly when n is EVEN (n=14: r=7).")
