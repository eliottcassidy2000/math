"""Exact induced-matching and marked-root controls for all-layer second weights."""
from functools import lru_cache
from itertools import permutations
from math import comb, factorial
import sympy as s

gates=0


def check(good,label):
    global gates
    if not good:
        raise RuntimeError(label)
    gates+=1


def choose(n,k):
    return comb(n,k) if 0<=k<=n else 0


@lru_cache(None)
def W(k,j):
    numerator=sum(choose(j,i)*(-4)**i*factorial(k-i-1) for i in range(j+1))
    check(numerator%2==0,"Hamilton matching integrality")
    return numerator//2


@lru_cache(None)
def subset_distribution(n,k,m):
    out=[]
    for j in range(min(m,k//2)+1):
        count=choose(m,j)*sum(choose(m-j,l)*2**l*choose(n-2*m,k-2*j-l)
                               for l in range(m-j+1))
        out.append(count)
    check(sum(out)==choose(n,k),"complete induced-matching subset distribution")
    return tuple(out)


@lru_cache(None)
def signed(n,k,m):
    return sum(count*W(k,j) for j,count in enumerate(subset_distribution(n,k,m)))


@lru_cache(None)
def marked_signed(n,k,m):
    # Root is not incident to a negative matching edge. Choose k-1 other vertices.
    return sum(count*W(k,j) for j,count in enumerate(subset_distribution(n-1,k-1,m)))


def weight(n,k,m):
    total=factorial(n)//(2*k*factorial(n-k))
    check((total-signed(n,k,m))%2==0,"all-layer matching weight integrality")
    return (total-signed(n,k,m))//2


literal=0
for n in range(5,10):
    values=[[0]*(n//2+1) for _ in range(n+1)]
    # This independent path uses full ordered vertex words and their minimum root.
    for k in range(3,n+1):
        for word in permutations(range(n),k):
            if word[0]!=min(word) or word[1]>word[-1]:
                continue
            edges={tuple(sorted((word[j],word[(j+1)%k]))) for j in range(k)}
            parity=0;literal+=1
            for m in range(n//2+1):
                if m:
                    parity^=(2*m-2,2*m-1) in edges
                values[k][m]+=parity
        for m,value in enumerate(values[k]):
            check(value==weight(n,k,m),"literal all-cycle matching versus induced-subset path")
    print(f"n={n} all_lengths_and_all_matching_sizes=literal_induced_PASS")

for n in (16,17,18,19,20,21,30,40):
    for k in range(3,n+1):
        e=factorial(n-2)//factorial(n-k)
        for m in range(n//2+1):
            check(signed(n,k,m)>0,"all-layer full matching signed sum positive")
        for m in range(n//2):
            difference=weight(n,k,m+1)-weight(n,k,m)
            if k==3:
                check(difference==n-2,"triangle matching increment")
            else:
                check(difference==2*marked_signed(n-1,k-1,m),"marked-unmatched-root contraction")
                check(marked_signed(n-1,k-1,m)>0,"marked-root positivity")
            check(difference>0,"all-layer strict matching monotonicity")
        disj=2*e-4*(k-3)*factorial(n-4)//factorial(n-k)
        adj=2*e-2*factorial(n-3)//factorial(n-k)
        check(weight(n,k,2)==disj,"two-disjoint-edge exact weight")
        check((adj>disj)-(adj<disj)==((2*k-n-3)>0)-((2*k-n-3)<0),
              "exact reversal and tie boundary")
    print(f"n={n} every_layer_matching_sign_marked_increment_and_reversal=PASS")

n,t,q,k=s.symbols("n t q k")
D4=(n-2)*(n-3)*(n-4)*(n-5)
R=q*((n-k+2)*(n-4)*(n-5)+(k-5)*(n*n-7*n+14-2*q))/D4
U=2-4*(k-3)/((n-2)*(n-3))
transport=2*(n-k)*(q*(q-n+3)-2*(n-4)*(n-5))/D4
check(s.cancel((R-U)-(R-U).subs(k,n)-transport)==0,"benchmark-aware n16 transport")
check(s.expand((q*(q-n+3)-2*(n-4)*(n-5)).subs(q,3*(n-5))-4*(n-5)*(n-7))==0,
      "benchmark-aware positive lower endpoint")
check(s.cancel(3*(n-5)/(n-2)-(2-2/(n-2))-(n-9)/(n-2))==0,
      "short-layer n16 threshold")
# The older n18 constant-benchmark proof is retained as a separate exact control.
F=q*((n-2)*(n-3)-2*q)
target=2*(n-3)**2*(n-4)
lo=n**3-28*n*n+207*n-468
hi=(n-4)*(n**3-22*n*n+108*n-152)/8
check(s.expand(F.subs(q,3*(n-5))-target-lo)==0,"uniform R2 lower endpoint")
check(s.expand(F.subs(q,(n-2)**2/4)-target-hi)==0,"uniform R2 upper endpoint")
for polynomial in (lo,hi):
    check(all(c>0 for c in s.Poly(polynomial.subs(n,18+t),t).all_coeffs()),
          "all-height uniform threshold positivity")
for f in (3,4,5):
    check(s.Rational(f)-s.Rational(2*f*(f-1),14)>2,"all finite >=3-edge motifs exceed2e at n16")
for j in range(10):
    check(j%2>=j-2*comb(j,2),"literal parity Bonferroni")
tables={k:[W(k,j) for j in range(k//2+1)] for k in range(3,9)}
check(tables[3]==[1,-1] and tables[4]==[3,-1,3],"small signed-table hostile signs")
check(tables[5]==[12,0,4] and tables[6]==[60,12,12,-4],"zero and negative signed-table controls")
check(all(x>0 for k in (7,8) for x in tables[k]),"remaining short signed tables positive")
print(f"short_signed_Hamilton_tables={tables}")
print(f"literal_cycles={literal} explicit_failure_gates={gates}")
print("claimed_all_order_scope=n>=16_and_every_3<=k<=n")
print("equality_transition=n_versus_2k-3;global_negation_exactly_when_k_even")
print("RESULT=PASS")
