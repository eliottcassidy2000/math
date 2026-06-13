"""C2b ANGLE 1b + minimizer structure (fast, no M_exact in hot loop). (1) discrete 2n-1 witness
works <=> shell-collision present? (2) the M=2/(2n-1) MINIMIZERS: do they always contain n and
its shell-partner n-1 (which collide forbidden sets, since (n-1)≡-n mod 2n-1)? opus-2026-06-03-S599k."""
from itertools import combinations
from fractions import Fraction as F
from math import gcd
def nrm(a,m): a%=m; return min(a,m-a)
def discrete_works(V,m): return any(min(nrm(v*j,m) for v in V)>=2 for j in range(1,m))
def has_shell(V,m): return any((V[i]+V[k])%m==0 for i in range(len(V)) for k in range(i+1,len(V)))
def M_exact(V):
    ds=set(V)
    for a,b in combinations(V,2): ds.add(a+b); ds.add(abs(a-b))
    ds.discard(0); best=F(-1)
    for d in ds:
        for mm in range(d):
            t=F(mm,d); val=min(min((v*t)%1,1-(v*t)%1) for v in V)
            if val>best: best=val
    return best
def main():
    print("ANGLE 1b — discrete 2n-1 witness works  <=>  shell-collision present?")
    print(" n | #C2b | works&shell | works&noshell | fails&shell | fails&noshell")
    for n in range(3,9):
        m=2*n-1; B=2*n; a=b=c=d=0
        for V in combinations(range(1,B+1),n-1):
            if not any(v%n==0 for v in V): continue
            g=0
            for v in V: g=gcd(g,v)
            if g!=1: continue
            w=discrete_works(V,m); s=has_shell(V,m)
            if w and s: a+=1
            elif w and not s: b+=1
            elif (not w) and s: c+=1
            else: d+=1
        print(f" {n:2d} | {a+b+c+d:5d} | {a:4d} | {b:4d} | {c:4d} | {d:4d}")
    print("\nMINIMIZER STRUCTURE — do all M=2/(2n-1) minimizers contain BOTH n and n-1 (shell-partners,")
    print(" since (n-1)+n=2n-1≡0; their forbidden j's ±n^{-1}=±2 and ±(n-1)^{-1}=∓2 COLLIDE)?")
    print(" n | minimizers (M=2/(2n-1)) | all contain {n, n-1}?")
    for n in range(3,8):
        m=2*n-1; B=3*n; thr=F(2,m); mins=[]
        for V in combinations(range(1,B+1),n-1):
            if not any(v%n==0 for v in V): continue
            g=0
            for v in V: g=gcd(g,v)
            if g!=1: continue
            if M_exact(V)==thr: mins.append(V)
        allcontain=all((n in V and (n-1) in V) for V in mins)
        print(f" {n:2d} | {len(mins)} mins, e.g. {mins[:3]} | {allcontain}")
if __name__=='__main__': main()
