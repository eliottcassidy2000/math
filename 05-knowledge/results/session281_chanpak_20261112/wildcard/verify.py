#!/usr/bin/env python3
"""Exact physical observation certificate. No assertion is disabled by -O."""
from pathlib import Path
import importlib.util,json
from fractions import Fraction
from math import gcd
from functools import lru_cache
ROOT=Path(__file__).resolve().parents[4]
p=ROOT/'05-knowledge/results/session275_exact_selection_20261111/switches/carrier.py'
spec=importlib.util.spec_from_file_location('inherited_carrier',p)
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)

def check(ok,label):
    if not ok:raise RuntimeError(label)

def hp(T):
    """Backward predecessor recurrence, independent of carrier's forward DP."""
    n=len(T);incoming=[sum(1<<u for u in range(n) if T[u]>>v&1) for v in range(n)]
    @lru_cache(None)
    def end(S,v):
        R=S^(1<<v)
        if not R:return 1
        pred=R&incoming[v];answer=0
        while pred:
            z=pred&-pred;pred-=z;answer+=end(R,z.bit_length()-1)
        return answer
    return sum(end((1<<n)-1,v) for v in range(n)) if n else 1

def delete(T,v):
    ids=[i for i in range(len(T)) if i!=v]
    return tuple(sum(1<<j for j in range(len(ids)) if T[ids[i]]>>ids[j]&1) for i in range(len(ids)))

def det(M):
    A=[[Fraction(x) for x in row] for row in M];ans=Fraction(1)
    for i in range(len(A)):
        j=next((j for j in range(i,len(A)) if A[j][i]),None)
        if j is None:return 0
        if j!=i:A[i],A[j]=A[j],A[i];ans=-ans
        z=A[i][i];ans*=z
        for j in range(i+1,len(A)):
            q=A[j][i]/z
            for k in range(i,len(A)):A[j][k]-=q*A[i][k]
    return int(ans) if ans.denominator==1 else str(ans)

I=(0,);S=(2,0);E=(0,1);EE=(0,1,3);SS=(6,4,0);C=(4,1,2)
PROBES=[(I,I),(I,E),(I,EE),(I,C),(S,I),(S,E),(C,I)]
KEYS=('H','h','Q','U','V','d','e')

def observe(T):return [hp(c.fuse(c.fuse(L,T),R)) for L,R in PROBES]
def recover(p):
    H=p[0]
    check(all(v%2==1 for v in p),'physical parity')
    U=Fraction(p[3]-p[0],2);V=Fraction(p[6]-p[0],2)
    d=p[0]+U-p[1];e=p[0]+V-p[4]
    h=Fraction(p[1]+p[3]-p[2]-p[0],2);Q=p[5]+p[0]-p[1]-p[4]
    vals=(H,h,Q,U,V,d,e)
    check(all(Fraction(v).denominator==1 for v in vals),'integral recovery')
    return dict(zip(KEYS,map(int,vals)))

def left_observe(T):return [hp(c.fuse(L,T)) for L in (I,S,SS,C)]
def left_recover(q):
    H=q[0];V=Fraction(q[3]-q[0],2);e=q[0]+V-q[1]
    h=Fraction(q[1]+q[3]-q[2]-q[0],2)
    return [H,h,V,e]

count=0
for n in range(2,6):
    for T in c.tournaments(n):
        p=c.profile(T);out=observe(T);got=recover(out)
        check(got=={k:p[k] for k in KEYS},f'profile recovery n={n},T={T}')
        check(gcd(*out)==gcd(*(p[k] for k in KEYS)),'full orbit content certificate')
        check(left_recover(left_observe(T))==[p[k] for k in ('H','h','V','e')],'left recovery')
        c.hp.cache_clear();count+=1

A=(0,41,3,5,15,29);B=(0,33,3,7,15,29)
pa=c.profile(A);pb=c.profile(B)
check(pa==pb,'full packet collision')
check(recover(observe(A))==recover(observe(B)),'literal probe collision')
ma=sorted(hp(delete(A,v)) for v in range(6));mb=sorted(hp(delete(B,v)) for v in range(6))
check(ma==[1,3,3,5,5,9] and mb==[1,1,5,5,5,9],'reroot menu hostile')

coef=[[1,0,0,0,0,0,0],[1,0,0,1,0,-1,0],[1,-2,0,3,0,-1,0],
 [1,0,0,2,0,0,0],[1,0,0,0,1,0,-1],[1,0,1,1,1,-1,-1],[1,0,0,0,2,0,0]]
check(abs(det(coef))==8,'seven-observation determinant')
Ts=[(0,1),(2,0),(0,1,3),(2,0,3),(4,1,2),(6,0,2),(2,0,3,7)]
transition=[ [c.profile(T)[k] for k in KEYS] for T in Ts]
check(det(transition)==16,'physical transition span')
Hankel=[[hp(c.fuse(L,R)) for R in (I,E,EE,C)] for L in (I,S,SS,C)]
check(det(Hankel)==-16,'physical Hankel rank four')
# k=1 boundary must not be silently reconstructed with k=0.
iobs=observe(I);iq=recover(iobs)
check(iq['Q']==2 and c.profile(I)['Q']==1,'singleton exposes k correction')
# Add eighth count and recover k= p_extra - p0 -3u-v +dd+ee -3Q.
# Here use full eight-coordinate observation matrix determinant instead.
coords=[(0,0),(1,1),(2,2),(0,1),(1,3),(0,2),(2,3),(0,3)]
basis=[]
for i,j in coords:
    Z=[[0]*4 for _ in range(4)];Z[i][j]=1
    if (i,j)==(0,0):Z[3][3]=1
    basis.append(Z)
probes8=PROBES[:6]+[(S,EE)]+PROBES[6:]
co8=[]
for L,R in probes8:
    ml=c.matrix(c.profile(L));mr=c.matrix(c.profile(R));row=[]
    for Z in basis:
        m=c.mul(c.mul(ml,Z),mr);row.append(m[0][0]+m[0][3])
    co8.append(row)
check(det(co8)==-16,'eight-entry determinant')
# Old central pair is a positive control, not new discovery.
C41=(56,13,25,48,34,6);C43=(52,13,24,49,34,6)
check([b-a for a,b in zip(observe(C41),observe(C43))]==[2]*7,'central delta experiment')
# Repaired bank from S280; recreate C2 directly through its one-apex word100.
# Root beats first spine vertex; spine is 1->2->3.
C2=(2,12,9,1)
R273=(386,508,489,497,485,193,257,320,32)
bank=[c.fuse(C2,C41),c.fuse(C2,C43),R273]
profiles=[c.profile(T) for T in bank]
check([p['H'] for p in profiles]==[269,271,273],'old bank counts')
check(all((p['h'],p['V'],p['e'])==(15,441,15) for p in profiles),'renewable receiver')
# Literal probes at order11 are affordable and independently test the four-count gate.
lobs=[left_observe(T) for T in bank]
check([[y-x for x,y in zip(lobs[0],lobs[i])] for i in (1,2)]==[[2]*4,[4]*4],'repair iff positive control')
# Six-bank upper receivers differ by -480; the C left probe amplifies that by 2.
# Exact formula rather than large literal order17 enumeration, inherited S280 inputs.
hs=[24875,24877,24879,24881,24883,24885]
vs=[28119]*3+[27639]*3
q3=[h+2*v for h,v in zip(hs,vs)]
check(q3[3]-q3[0]==6-960,'six bank receiver hostile')
result=dict(status='PASS',literal_nontrivial_inputs=count,probe_added_vertices_max=2,
 observation_det=det(coef),transition_span_det=det(transition),hankel=Hankel,hankel_det=det(Hankel),
 full_eight_entry_det=det(co8),collision_packet=pa,collision_root_menus=[ma,mb],
 singleton_observations=iobs,renewable_bank_left_observations=lobs,six_bank_probe_delta=q3[3]-q3[0])
print(json.dumps(result,indent=2,sort_keys=True))
