#!/usr/bin/env python3
"""Is the function-field INV unconditional?  Test over F_5.  (boxeph-S91)"""
from itertools import product, combinations
P=5
def norm(a):
    a=list(a)
    while a and a[-1]%P==0: a.pop()
    return tuple(x%P for x in a)
def deg(a):
    a=norm(a); return len(a)-1
def pmul(a,b):
    a=norm(a); b=norm(b)
    if not a or not b: return ()
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): r[i+j]+=x*y
    return norm(r)
def _degl(a):
    i=len(a)-1
    while i>=0 and a[i]%P==0: i-=1
    return i
def pmod(a,b):
    a=list(norm(a)); b=norm(b); lb=len(b); inv=pow(b[-1],P-2,P)
    while _degl(a)>=lb-1:
        da=_degl(a); d=da-(lb-1); c=(a[da]*inv)%P
        for i,y in enumerate(b): a[d+i]=(a[d+i]-c*y)%P
    return norm(a)
def pgcd(a,b):
    a=norm(a); b=norm(b)
    while b: a,b=b,pmod(a,b)
    return a
def monic(a):
    a=norm(a)
    if not a: return a
    inv=pow(a[-1],P-2,P); return norm(tuple((x*inv)%P for x in a))
def polys_upto_deg(d):
    seen=set(); res=[]
    for L in range(1,d+2):
        for c in product(range(P),repeat=L):
            a=norm(c)
            if a and a not in seen: seen.add(a); res.append(a)
    return res
def M_logdist(V,maxdegQ=2):
    best=None; info=None
    for Q in [q for q in polys_upto_deg(maxdegQ) if deg(q)>=1 and q==monic(q)]:
        dQ=deg(Q)
        for a in polys_upto_deg(dQ-1):
            if deg(a)>=dQ or pgcd(a,Q)!=(1,): continue
            worst=0; ok=True; res={}
            for v in V:
                r=pmod(pmul(v,a),Q); res[v]=r
                if not r: ok=False; break
                worst=min(worst,deg(r)-dQ)
            if not ok: continue
            if best is None or worst>best: best=worst; info=(a,Q,dict(res))
    return best,info
def vanishing_poly():
    r=[0]*(P+1); r[P]=1; r[1]=(r[1]-1)%P; return norm(r)
def evalp(a,c):
    s=0
    for x in reversed(norm(a)): s=(s*c+x)%P
    return s
def roots(a): return [c for c in range(P) if evalp(a,c)==0]
def is_covering(V):
    cov=set()
    for v in V: cov|=set(roots(v))
    return cov==set(range(P))
def topc(r,dQ): return r[dQ-1] if len(r)>=dQ else 0

print("THE F_5 DEEP WELL: F_5^* (constants) + killer t^5-t")
V=[(c,) for c in range(1,5)]+[vanishing_poly()]
v1,_=M_logdist(V,1); v2,i2=M_logdist(V,2)
print(f"  covering={is_covering(V)}  level1={v1} (None=no level1=covering)  level2={v2}")
if i2:
    a,Q,res=i2; dQ=deg(Q); tops=[topc(res[v],dQ) for v in V[:-1]]
    print(f"  a={a} Q={Q}; core top-coeffs={sorted(tops)}==F_5^*? {sorted(tops)==[1,2,3,4]}; killer top={topc(res[V[-1]],dQ)}")

print("\nHUNT: covering + level-2-only-lonely 5-families; is the 4-core's top-coeff set = F_5^* (distinct nonzero)?")
tested=0; inv_ok=0; bad=[]
lowpolys=[p for p in polys_upto_deg(1)]  # deg<=1
for core in combinations(lowpolys,4):
    V=list(core)+[vanishing_poly()]
    if len(set(V))!=5 or not is_covering(V): continue
    v1,_=M_logdist(V,1)
    if v1==-1: continue
    v2,i2=M_logdist(V,2)
    if v2!=-1 or i2 is None: continue
    tested+=1
    a,Q,res=i2; dQ=deg(Q); tops=[topc(res[v],dQ) for v in V[:-1]]
    if len(set(tops))==4 and 0 not in tops: inv_ok+=1
    else: bad.append(([c for c in core],sorted(tops)))
print(f"  tested {tested}; 4-core top-coeffs = F_5^* (4 distinct nonzero): {inv_ok}/{tested}")
for core,tops in bad[:8]: print(f"    core={core} tops={tops}")
print("\nVERDICT:", "UNCONDITIONAL (FF-INV holds by packing)" if bad==[] else "has residual (packing needs the difference-closed hypothesis)")
