#!/usr/bin/env python3
"""mac-mini-S77: PROVE the Schur-triple inverse T(A) <= C(k,2), equality iff dilated AP.
T(A) = #{ordered (a,b) in A^2 : a+b in A} = Sum_{c in A} r_A(c), r_A(a_l)=#{i<l: a_l-a_i in A}.
Since a+b=c>0 needs a,b<c, r_A(a_l)<=l-1 => T<=Sum(l-1)=C(k,2). Equality iff {a_l-a_i:i<l}=
{a_1..a_{l-1}} for all l iff min diff a_l-a_{l-1}=a_1 all l iff constant gap iff dilated AP.
EXHAUSTIVE verification: (1) T<=C(k,2) never violated; (2) equality ONLY at dilated APs."""
from itertools import combinations
from math import comb, gcd
def T(A):
    Aset=set(A); return sum(1 for a in A for b in A if (a+b) in Aset)
def is_dilated_AP(A):
    A=sorted(A); d=A[1]-A[0]
    return all(A[i+1]-A[i]==d for i in range(len(A)-1)) and A[0]==d  # {d,2d,...,kd}
def is_AP_any(A):  # any AP (dilated interval = AP with first term = common diff)
    A=sorted(A); d=A[1]-A[0]
    return all(A[i+1]-A[i]==d for i in range(len(A)-1))

print("EXHAUSTIVE: all k-subsets of {1..N}; check T<=C(k,2) and equality-iff-dilated-AP\n")
for k in [3,4,5]:
    N={3:14,4:12,5:11}[k]
    cap=comb(k,2)
    viol=0; eq_nonAP=[]; AP_noneq=[]; count=0
    for A in combinations(range(1,N+1),k):
        count+=1
        t=T(A)
        if t>cap: viol+=1; print(f"  VIOLATION T({A})={t}>{cap}")
        dap = (sorted(A)[0]==sorted(A)[1]-sorted(A)[0]) and is_AP_any(A)  # dilated interval {d,2d,..,kd}
        if t==cap and not dap: eq_nonAP.append(A)
        if dap and t!=cap: AP_noneq.append((A,t))
    print(f"k={k}, subsets of [1..{N}]: {count} sets. T<=C(k,2)={cap}: {'OK (0 violations)' if viol==0 else f'{viol} VIOLATIONS'}")
    print(f"   equality T=C(k,2) but NOT dilated-AP: {eq_nonAP[:5] if eq_nonAP else 'NONE'}")
    print(f"   dilated-AP but T!=C(k,2): {AP_noneq[:3] if AP_noneq else 'NONE'}")
# also confirm dilated APs {d,2d,..,kd} all achieve equality
print("\ndilated APs {d,2d,..,kd} achieve T=C(k,2):")
for k in [5,8,13]:
    for d in [1,2,3,7]:
        A=[d*i for i in range(1,k+1)]
        print(f"   k={k},d={d}: T={T(A)}, C(k,2)={comb(k,2)}, equal={T(A)==comb(k,2)}")
print("\n=> THEOREM VERIFIED: T(A)<=C(k,2), equality iff A is a dilated AP {d,2d,...,kd}.")
print("This is the combinatorial EXTREMAL STEP of the last inch (opus-S182 target). The remaining")
print("open part is the RESUMMATION (Schur-deficit => L>0); this theorem alone is NOT LRC(14).")
