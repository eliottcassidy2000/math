"""Exact all-cycle transposition controls, retaining one distinguished vertex."""
from fractions import Fraction
from itertools import combinations, permutations
from math import factorial, perm
import sympy as s

gates=0


def check(good,label):
    global gates
    if not good:
        raise RuntimeError(label)
    gates+=1


def identity(left,right,label):
    check(s.cancel(left-right)==0,label)


def formula(n,k,r):
    q=r*(n-2-r)
    e=factorial(n-2)//factorial(n-k)
    if k<=5:
        return 2*q*e//(n-2)
    B=n*n-7*n+14-2*q
    bracket=(n-k+2)*(n-4)*(n-5)+(k-5)*B
    return 2*q*factorial(n-6)//factorial(n-k)*bracket


cycles_total=0
for n in range(5,10):
    edge_index={edge:j for j,edge in enumerate(combinations(range(n),2))}
    star=[]
    for r in range(n-1):
        star.append(sum((1<<edge_index[(u,x)]) for u in (0,1) for x in range(2,r+2)))
    counts=[[0]*(n-1) for _ in range(n+1)]
    one=[[0]*(n-1) for _ in range(n+1)]
    both=[[0]*(n-1) for _ in range(n+1)]
    total=0
    for k in range(3,n+1):
        for chosen in combinations(range(n),k):
            root=chosen[0]
            distinguished=(0 in chosen)+(1 in chosen)
            for p in permutations(chosen[1:]):
                if p[0]>p[-1]:
                    continue
                word=(root,)+p
                mask=sum(1<<edge_index[tuple(sorted((word[j],word[(j+1)%k])))]
                         for j in range(k))
                total+=1
                for r in range(n-1):
                    negative=(mask&star[r]).bit_count()%2
                    counts[k][r]+=negative
                    if distinguished==1:
                        one[k][r]+=negative
                    elif distinguished==2:
                        both[k][r]+=negative
        for r in range(n-1):
            q=r*(n-2-r)
            check(one[k][r]==2*q*perm(n-4,k-3),"exactly one distinguished vertex")
            if k==3:
                expect_both=0
            elif k==4:
                expect_both=2*q
            elif k==5:
                expect_both=4*q*perm(n-4,k-4)
            else:
                B=n*n-7*n+14-2*q
                expect_both=4*q*perm(n-4,k-4)+2*q*B*(k-5)*perm(n-6,k-6)
            check(both[k][r]==expect_both,"both distinguished vertices with near/far distinction")
            check(counts[k][r]==formula(n,k,r),"full literal transposition weight")
    check(all(value==0 for value in one[n]),"Hamilton one-vertex term vanishes")
    check(one[3][2]>0 and both[3][2]==0,"both-only counting hostile in triangle layer")
    cycles_total+=total
    print(f"n={n} all_cycle_lengths=3..{n} cycles={total} every_split_and_distinguished_stratum=PASS")

n,k,q=s.symbols("n k q")
D=(n-2)*(n-3)*(n-4)*(n-5)
B=n*n-7*n+14-2*q
R=q*((n-k+2)*(n-4)*(n-5)+(k-5)*B)/D
H=q*((n-2)*(n-3)-2*q)/((n-2)*(n-3)*(n-4))
identity(R.subs(k,n),H,"Hamilton boundary")
identity(R.subs(k,5),q/(n-2),"length-five boundary")
identity(R-H,2*q*(n-k)*(q-n+3)/D,"normalized Hamilton dominance k>=5")
identity(q/(n-2)-H,2*q*(q-n+3)/((n-2)*(n-3)*(n-4)),
         "normalized short-layer dominance")
identity(R.subs(q,2*(n-4)),2*(n*n-7*n+22-2*k)/((n-2)*(n-3)),
         "two-disagreement lower endpoint")
for order in range(9,151):
    for length in range(3,order+1):
        e=factorial(order-2)//factorial(order-length)
        for r in (2,(order-2)//2,order-4):
            norm=Fraction(formula(order,length,r),2*e)
            ham=Fraction(formula(order,order,r),2*factorial(order-2))
            check(norm>=ham,"finite auxiliary uniform dominance")
print(f"directly_enumerated_cycles={cycles_total} explicit_failure_gates={gates}")
print("uniform_layer_scope=all_3<=k<=n_at_n>=9")
print("missing_term=2*r*(n-2-r)*falling(n-4,k-3)")
print("RESULT=PASS")
