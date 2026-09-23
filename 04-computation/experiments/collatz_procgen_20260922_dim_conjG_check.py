"""Conjecture G check beyond the certified census: for each j in J_LIST, take the candidates x=-p/3^j with
1 <= |x| < 3/2 lying in classes of Bad_M (thread dump), find the nearest (2-adically) point h of smaller
exponent among the census (CENSUS_FILE) and the candidates of smaller j, and test whether x is the canonical
Theorem-P step x = h - 2^i/3^alpha_h(i) with condition (A).  When every parent is certified hostile, such an x
is PROVED hostile by Theorem P (this is how the one prover-undecided point at j=27 is settled).
usage: python3 ..._dim_conjG_check.py DUMPFILE M CENSUS_FILE ALPHA_BIN j1,j2,...
"""
import subprocess, sys
from fractions import Fraction as F
dump,M,cfile,alpha_bin,jl=sys.argv[1],int(sys.argv[2]),sys.argv[3],sys.argv[4],[int(t) for t in sys.argv[5].split(',')]
MOD=1<<M
bad=[int(l.split()[0]) for l in open(dump)]
cen=[F(l.strip()) for l in open(cfile) if l.strip()]
def j3(x):
    q=x.denominator; e=0
    while q%3==0: q//=3; e+=1
    return e
def v2(n): n=abs(n); return (n & -n).bit_length()-1
newpts={}
for j in jl:
    q=3**j; L=[]
    for x in bad:
        p=(-x*q)%MOD
        if q<=p<q*3//2 and p%2==1 and p%3: L.append(F(-p,q))
    newpts[j]=sorted(set(L))
pool=list(cen); rows=[]
for j in jl:
    for x in newpts[j]:
        best=None
        for h in pool:
            if j3(h)>=j: continue
            i=v2((x-h).numerator)
            if best is None or i>best[0]: best=(i,h)
        rows.append((x,j,best[1],best[0]))
    pool+=newpts[j]
need={}
for x,j,h,i in rows: need[h]=max(need.get(h,0),i)
hs=list(need)
out=subprocess.run([alpha_bin,'26',str(max(need.values()))]+[t for h in hs for t in (str(h.numerator),str(h.denominator))],capture_output=True,text=True,check=True).stdout.strip().splitlines()
al={h:[int(t) for t in line.split('alpha:')[1].split()] for h,line in zip(hs,out)}
for j in jl:
    R=[r for r in rows if r[1]==j]
    canon=[r for r in R if r[0]==r[2]-F(2**r[3],3**al[r[2]][r[3]-1])]
    A=[r for r in canon if F(1,2)-F(1,2*3**j) > (-r[2]-1)+F(2**r[3],3**j)]
    print(f"j={j}: candidates with 1<=|x|<3/2 in Bad_{M}: {len(R)}; canonical steps from the nearest lower-exponent point: {len(canon)}; satisfying (A): {len(A)}; deficit 1.585j-i: {sorted(set(round(1.585*j-r[3]) for r in R))}")
    nc=[str(r[0]) for r in R if r not in canon]
    if nc: print("   non-canonical:",nc)
