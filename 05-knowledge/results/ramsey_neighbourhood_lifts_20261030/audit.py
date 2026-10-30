"""Separate exact audits; does not import the producer or its chart generator."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
GATES=0
def need(ok,label):
    global GATES
    GATES+=1
    if not ok:raise RuntimeError(label)

def main():
    data=json.loads((HERE/'verification.json').read_text())
    # Native cycle test recovers the entire maximal-support family.
    for n in range(3,12,2):
        arc=lambda i,j:0<(j-i)%n<=n//2
        cycles=[sum(1<<i for i in s) for s in combinations(range(n),3)
                if all(sum(arc(i,j) for j in s if j!=i)==1 for i in s)]
        trans=[m for m in range(1<<n) if not any(m&t==t for t in cycles)]
        good=set(trans)
        maximal={m for m in trans if not any(m|(1<<i) in good for i in range(n) if not m>>i&1)}
        intervals={sum(1<<((s+i)%n) for i in range(n//2+1)) for s in range(n)}
        need(maximal==intervals,'literal carousel supports are exactly half-circle intervals')
    # Direct summatory evaluation, not the binary recurrence implementation.
    g=[sum(2**i.bit_count() for i in range(k)) for k in range(257)]
    for a in range(65):
        for b in range(65):
            need(g[a]+g[b]<=g[a+b],'G superadditivity')
            for c in range(33):need(g[a]+g[b]+g[c]<=g[max(a+b,a+c,b+c)],'three-block integer law')
    # Classical maximum recurrence is an independent way to compute w.
    w=[F(0),F(1)]
    for k in range(2,257):
        w.append(max(F(3,2)*w[min(j,k-j)]+w[max(j,k-j)] for j in range(1,k)))
        h=k//2
        need(w[k]==(F(5,2)*w[h] if k%2==0 else F(3,2)*w[h]+w[h+1]),'balanced maximum recurrence')
    need([str(x) for x in w[1:33]]==data['graphs']['slopes_1_to_32'],'slope receipt')
    # Integer scaling constructions, checked without materializing huge graphs.
    @lru_cache(None)
    def witness(a,b):
        if a==1:return b,1,b
        need(b%2==0,'realization divisibility guard')
        h=a//2;aa=[h]*5
        if a%2:aa[0]=aa[2]=h+1
        kids=[witness(x,b//2) for x in aa]
        return (sum(x[0] for x in kids),max(kids[i][1]+kids[(i+1)%5][1] for i in range(5)),
                max(kids[i][2]+kids[(i+2)%5][2] for i in range(5)))
    for a in range(1,257):
        b=1<<(a-1).bit_length();n,om,al=witness(a,b)
        need((om,al)==(a,b) and n==w[a]*b,'all-slope recursive realization')
    # Residue closure needed for disjoint unions in the exact first four rows.
    h=lambda b:25*(b//4)+(0,4,10,16)[b%4]
    f=lambda b:4*b-int(b%4!=0)
    for x in range(1,17):
        for y in range(1,17):
            need(h(x)+h(y)<=h(x+y) and f(x)+f(y)<=f(x+y),'residue superadditivity')
    for b in range(3,103):need(11*b//2<=h(b),'nonuniform root bound')
    # Exhaust all 243 possible positive omega patterns at budget four.
    for a in product(range(1,4),repeat=5):
        if all(a[i]+a[(i+1)%5]<=4 for i in range(5)) and a!=(2,)*5:
            need(sum(w[x] for x in a)<=11,'all nonuniform root patterns')
    # Handshake defect cost also has an exact integer form, no floating point.
    for d in range(35):
        h0=34-d;p=6 if h0==34 else 3*(h0%2)
        value=3*h0*h0+(30+5-3*87-3*6)*h0+(87-30-1)*(87-1)+p
        need(value==data['neighbourhoods']['F_h0_to_34'][h0],'independent neighbourhood polynomial')
        need(value+6>=40*d,'defect potential')
    # Deliberately remove a needed guard: isolated vertices allow unbounded x
    # if only pair inequalities are retained. Singletons must remain in G1.
    need(w[2]>w[1],'missing singleton guard hostile')
    need(data['graphs']['capacities_A1_to_4_B1_to_8'][2][1]==7<4*2,'integral bound eight is unrealizable')
    print(json.dumps({'status':'PASS','independent_gates':GATES,'largest_recursive_slope_budget':256}))
if __name__=='__main__':main()
