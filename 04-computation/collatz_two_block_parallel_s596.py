#!/usr/bin/env python3
"""Collatz <-> LRC parallel. Collatz cycle on odd a_1->...->a_k (a->(3a+1)/2^e):
a_1(2^E - 3^k) = S, S=sum_{i} 3^{k-i} 2^{sigma_{i-1}}, E=sum e_i. Two-block: 2^E (2-adic)
vs 3^k (odd); obstruction = 2^E-3^k divides the bounded sum S -- the SAME shape as the
LRC S595 rank-1 two-block (det = w n u u' l). opus-2026-06-03-S596."""
from itertools import product
from math import gcd, log2
def cycle_solution(es):
    k=len(es); E=sum(es)
    den=2**E-3**k
    if den<=0: return None
    # S = sum_{i=1..k} 3^{k-i} 2^{sigma_{i-1}}, sigma_0=0
    S=0; sig=0
    for i in range(k):
        S+=3**(k-1-i)*2**sig; sig+=es[i]
    if S%den!=0: return None
    a1=S//den
    if a1<=0 or a1%2==0: return None
    # verify it's a genuine cycle of odd numbers
    a=a1; seq=[a]
    for e in es:
        x=3*a+1
        if x % (2**e)!=0: return None
        a=x//(2**e)
        if a%2==0 and a!=a1: pass
        seq.append(a)
    if a!=a1: return None
    return a1,E,den,S,seq[:-1]
def main():
    print("Collatz cycle two-block: a_1*(2^E - 3^k) = S; valid cycle needs (2^E-3^k) | S, a_1 odd>0")
    found=[]
    for k in range(1,8):
        for es in product(range(1,6),repeat=k):
            sol=cycle_solution(es)
            if sol:
                a1,E,den,S,seq=sol
                if True:
                    found.append((k,E,den,a1,tuple(sorted(set(seq)))))
    # report distinct cycles
    seen=set(); out=[]
    for k,E,den,a1,verts in found:
        if verts in seen: continue
        seen.add(verts); out.append((k,E,den,a1,verts))
    print(f"  cycles found (e_i in 1..5, k<=7): {len(out)}")
    for k,E,den,a1,verts in out[:6]:
        print(f"   k={k} E={E}: 2^E-3^k={den}; a_1={a1}; cycle odds={verts}")
    print()
    print("The two-block gap 2^E - 3^k (the obstruction denominator) near 0 = where cycles could hide:")
    for k in range(1,15):
        E=round(k*log2(3))  # E ~ k log2(3) makes 2^E ~ 3^k
        for Ee in (E,E+1):
            d=2**Ee-3**k
            print(f"   k={k:2d}: E={Ee}, 2^E-3^k={d}  (|.|/2^E={abs(d)/2**Ee:.4f})")
if __name__=='__main__': main()
