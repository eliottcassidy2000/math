from fractions import Fraction as F
from itertools import combinations
from math import comb, gcd

def forbidden(v, delta):
    r=delta/v; out=[]
    for k in range(v):
        c=F(k,v); L=(c-r)%1; H=(c+r)%1
        if L<H: out.append((L,H))
        else: out+=[(L,F(1)),(F(0),H)]
    return merge(out)
def merge(iv):
    iv=sorted(iv); res=[]
    for l,h in iv:
        if res and l<=res[-1][1]: res[-1]=(res[-1][0],max(res[-1][1],h))
        else: res.append((l,h))
    return res
def depth_dist(V, delta):
    ev={}
    for v in V:
        for (lo,hi) in forbidden(v,delta):
            ev[lo]=ev.get(lo,0)+1; ev[hi]=ev.get(hi,0)-1
    pts=sorted(set(list(ev)+[F(0),F(1)])); pk={}; d=0; prev=F(0)
    for p in pts:
        if p>prev: pk[d]=pk.get(d,F(0))+(p-prev)
        d+=ev.get(p,0); prev=p
    return pk
def Krav(k,n,x):  # K_k(n,x) = sum_j (-1)^j C(x,j) C(n-x,k-j)
    return sum((-1)**j * comb(x,j) * comb(n-x,k-j) for j in range(k+1))

# n=14: 13 runners, delta=1/14
n=13; delta=F(1,14)
configs={
 "WALL {1..11,13,14}":[1,2,3,4,5,6,7,8,9,10,11,13,14],
 "AP {1..13} (no apex)":list(range(1,14)),
 "random {1,2,4,7,8,11,13,16,17,19,22,23,25}":[1,2,4,7,8,11,13,16,17,19,22,23,25],
}
print(f"n={n} delta={delta}  2delta=1/7  baseline (1-4d)=5/7  indep p0=(6/7)^13={float((F(6,7))**13):.5f}")
print(f"{'config':30} {'p0_direct':>10} {'(1/2^n)Sum rho_j':>16} {'match':>6}")
data={}
for name,V in configs.items():
    pk=depth_dist(V,delta)
    p0=pk.get(0,F(0))
    # Krawtchouk transform rho_j = sum_k K_j(n,k) p_k
    rho=[sum(Krav(j,n,k)*pk.get(k,F(0)) for k in range(n+1)) for j in range(n+1)]
    recon=sum(rho)/(2**n)
    data[name]=(pk,rho,p0)
    print(f"{name:30} {float(p0):>10.5f} {float(recon):>16.5f} {str(p0==recon):>6}")

print("\n=== resonance levels rho_j vs independent baseline C(13,j)(5/7)^j (deviation = resonance) ===")
base=[comb(n,j)*float(F(5,7))**j for j in range(n+1)]
for name,(pk,rho,p0) in data.items():
    print(f"\n{name}:")
    print("  j:       "+" ".join(f"{j:>8}" for j in range(7)))
    print("  rho_j:   "+" ".join(f"{float(rho[j]):>8.3f}" for j in range(7)))
    print("  baseline:"+" ".join(f"{base[j]:>8.3f}" for j in range(7)))
    print("  excess:  "+" ".join(f"{float(rho[j])-base[j]:>+8.3f}" for j in range(7)))

# Normalized lonely measure: p0/indep_p0  (how much the resonance helps/hurts vs independent)
print("\n=== Krawtchouk-normalized lonely ratio  p0 / (6/7)^13 ===")
ip0=float(F(6,7)**13)
for name,(pk,rho,p0) in data.items():
    print(f"  {name:30} ratio = {float(p0)/ip0:>7.3f}  ({'lonelier' if float(p0)>ip0 else 'less lonely'} than independent)")
