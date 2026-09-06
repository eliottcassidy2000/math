"""Exact structural-carrier controls for cumulative second cycle weights."""
from functools import lru_cache
from itertools import combinations
from math import comb, factorial
from fractions import Fraction
import sympy as s

gates=0


def check(good,label):
    global gates
    if not good:
        raise RuntimeError(label)
    gates+=1


@lru_cache(None)
def fact(n):
    return factorial(n)


def matching(n,k,m):
    # Exact cycle edge-subset expansion, independent of induced-subset averaging.
    return sum((-4)**(j-1)*comb(m,j)*comb(n-2*j,k-2*j)*fact(k-j-1)
               for j in range(1,min(m,k//2)+1))


def motif(n,k,edges):
    total=0
    for j in range(1,len(edges)+1):
        for selected in combinations(edges,j):
            vertices={v for edge in selected for v in edge}
            adj={v:set() for v in vertices}
            for u,v in selected:
                adj[u].add(v);adj[v].add(u)
            if max(map(len,adj.values()))>2:
                continue
            unseen=set(vertices);components=[]
            while unseen:
                start=unseen.pop();component={start};todo=[start]
                while todo:
                    u=todo.pop()
                    for v in adj[u]&unseen:
                        unseen.remove(v);component.add(v);todo.append(v)
                components.append(component)
            cyclic=any(all(len(adj[v])==2 for v in component) for component in components)
            if cyclic:
                count=int(len(components)==1 and k==j)
            else:
                h=len(components)
                count=(2**(h-1)*comb(n-j-h,k-j-h)*fact(k-j-1)
                       if k>=j+h else 0)
            total+=(-2)**(j-1)*count
    return total


motifs={"P3":((0,1),(1,2)),"P4":((0,1),(1,2),(2,3)),
        "C3":((0,1),(1,2),(0,2)),"C4":((0,1),(1,2),(2,3),(0,3)),
        "C5":((0,1),(1,2),(2,3),(3,4),(0,4))}
anti_wins=0;edge_ties=0;tested_pairs=0;central_winners=0
for n in range(16,61):
    labels=[f"M{m}" for m in range(n//2+1)]+list(motifs)
    positive={label:0 for label in labels};negative={label:0 for label in labels}
    E=anti=disj=0
    for k in range(3,n+1):
        e=fact(n-2)//fact(n-k)
        N=fact(n)//(2*k*fact(n-k))
        E+=e
        if k%2:
            anti+=N
        disj+=2*e-4*(k-3)*fact(n-4)//fact(n-k)
        for label in labels:
            w=matching(n,k,int(label[1:])) if label.startswith("M") else motif(n,k,motifs[label])
            positive[label]+=w;negative[label]+=N-w if k%2 else w
        if k<4:
            continue
        tested_pairs+=1
        adj=Fraction(2*(n-3),n-2)*E
        expected=min(adj,disj,anti)
        spectrum={label:value for label,value in positive.items() if label not in ("M0","M1")}
        spectrum.update({"minus_"+label:value for label,value in negative.items()})
        actual=min(spectrum.values())
        check(actual==expected,"complete structural-carrier cumulative minimum")
        expected_ties=set()
        if adj==expected:expected_ties.add("P3")
        if disj==expected:expected_ties.add("M2")
        if anti==expected:expected_ties.add("minus_M0")
        check({label for label,value in spectrum.items() if value==actual}==expected_ties,
              "complete cumulative equality types")
        check(actual>E,"first-weight separation")
        check(Fraction(actual,E)>=Fraction(2*n*(n-1),(n+1)**2),"uniform cumulative second isolation")
        if k%2:
            check(anti>2*E,"odd-top antibalanced exclusion")
        if expected_ties=={"minus_M0"}:anti_wins+=1
        if len(expected_ties)>1:edge_ties+=1
        if k%2==0 and abs(2*(k-1)-(n+1))<=2:
            upper=Fraction(2*n*(n*n-1),(n+3)*(n*n-13))
            check(Fraction(anti,E)<=upper,"sharp constant-six central witness bound")
            check(expected_ties=={"minus_M0"},"all-order central antibalanced winner control")
            central_winners+=1
        if n==16 and k==10:
            check((E,anti,adj,disj)==(140807044,234801776,261498796,260444884),
                  "antibalanced defeats both two-edge competitors")
            print(f"hostile_n16_top10 E={E} anti={anti} adjacent={adj} disjoint={disj}")

n,x,t=s.symbols("n x t")
check(s.expand(n*(n-1)*(n-2)-(3*n-2)**2-(n**3-12*n*n+14*n-4))==0,
      "negative-matching paired polynomial")
check(all(c>0 for c in s.Poly((n**3-12*n*n+14*n-4).subs(n,16+t),t).all_coeffs()),
      "all-height paired-matching positivity")
check(s.expand((3*n-2)**2-8*x*(3*n-2-2*x)-(4*x-3*n+2)**2)==0,
      "quadratic pairing bound")
check(s.expand(n*(n-1)-(n+3)**2/2-(n-9)*(n+1)/2)==0,
      "negative-single-edge paired positivity")
check(s.cancel((2-4/(n-2))-2*n*(n-1)/(n+1)**2-
               2*(n*n-9*n-4)/((n-2)*(n+1)**2))==0,
      "antibalanced ratio is common lower envelope")
check(all(c>0 for c in s.Poly((n*n-9*n-4).subs(n,16+t),t).all_coeffs()),
      "common lower envelope all-height positivity")
gamma=2*n*(n-1)/(n+1)**2
upper=2*n*(n*n-1)/((n+3)*(n*n-13))
eta=2-2*(n-1)/((n-2)*(n-3))
check(s.cancel(gamma/((1-4/(n+1)**2)*(1-12/(n*n-1)))-upper)==0,
      "sharp worst-scale upper bound simplification")
numerator=2*(2*n**4-29*n**3+55*n*n+149*n-273)
check(s.cancel(eta-upper-numerator/((n-3)*(n-2)*(n+3)*(n*n-13)))==0,
      "central antibalanced strict margin identity")
check(all(c>0 for c in s.Poly(numerator.subs(n,16+t),t).all_coeffs()),
      "central antibalanced winner at every order")
check(s.limit(n*(2-gamma),n,s.oo)==s.limit(n*(2-upper),n,s.oo)==6,
      "sharp asymptotic worst-scale isolation constant")
print(f"structural_universe=n16..60_all_top_lengths4..n cases={tested_pairs}")
print(f"strict_antibalanced_winners={anti_wins} multiple_type_ties={edge_ties}")
print(f"central_even_top_winners={central_winners} sharp_worst_scale_constant=6")
print(f"explicit_failure_gates={gates}")
print("RESULT=PASS")
