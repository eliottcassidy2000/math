"""Exact checks for six-state carriers; run from any working directory.
No theorem is inferred from the finite bounds below. Uses the standard library.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import product, permutations
from math import prod
from pathlib import Path
import json

OUT = Path(__file__).with_name('receipts.json')

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def six(a,b):
    a,b=F(a),F(b)
    return (a,b,b/a,1/a,1/b,a/b)

def orbit(pair):
    out=[]
    while pair not in out:
        out.append(pair)
        pair=(pair[1],pair[1]/pair[0])
    require(pair==out[0], 'pair==out[0]')
    return out

def factor(n):
    out={}
    p=2
    while p*p<=n:
        while n%p==0:
            out[p]=out.get(p,0)+1
            n//=p
        p+=1
    if n>1: out[n]=out.get(n,0)+1
    return out

def divisors(n):
    out=[1]
    for p,c in factor(n).items():
        out=[d*p**j for d in out for j in range(c+1)]
    return sorted(out)

def integral_states(C):
    # Independent search in h0,h1; opposite-pair product forces both to divide C^2.
    states=[]
    for h0,h1 in product(divisors(C*C),repeat=2):
        vals=tuple(C*x for x in six(F(h0,C),F(h1,C)))
        if all(v.denominator==1 for v in vals):
            vals=tuple(map(int,vals))
            require(all(vals[(i-1)%6]*vals[(i+1)%6]==C*vals[i] for i in range(6)), 'all(vals[(i-1)%6]*vals[(i+1)%6]==C*vals[i] for i in range(6))')
            states.append(vals)
    return states

def hexagon(m):
    return {(r,s) for r in range(2*m+1) for s in range(2*m+1) if m<=r+s<=3*m}

def step(p,m):
    r,s=p
    return 2*m-s,r+s-m

def marked_order(p,m):
    r,s=p
    return tuple([('A',i) for i in range(r)]+[('B',i) for i in range(s)]+[('A',i) for i in range(r,2*m)]+[('B',i) for i in range(s,2*m)]),r+s

def cycle_from(f,start):
    out=[]
    n=start
    while n not in out:
        out.append(n)
        n=f(n)
    require(n==start, 'n==start')
    return out

def raw(n): return 3*n-1 if n%2 else n//2

def terras(n): return (3*n-1)//2 if n%2 else n//2

def odd_step(n):
    n=3*n-1
    while n%2==0: n//=2
    return n

def adjacency_cyclic_book(k):
    # C3[1,1,k], ordered blocks X,Y,Z with X->Y->Z->X; Z transitive.
    n=k+2
    a=[[False]*n for _ in range(n)]
    a[0][1]=True
    for j in range(2,n):
        a[1][j]=True
        a[j][0]=True
        for t in range(j+1,n): a[j][t]=True
    for i in range(n):
        for j in range(i+1,n):
            if not a[i][j] and not a[j][i]: a[j][i]=True
    return a

def ordinal(blocks):
    n=sum(map(len,blocks))
    a=[[False]*n for _ in range(n)]
    start=0
    for block in blocks:
        m=len(block)
        for i in range(m):
            for j in range(m): a[start+i][start+j]=block[i][j]
            for j in range(start+m,n): a[start+i][j]=True
        start+=m
    return a

def hamilton(a):
    # End-state subset DP; does not use SCC/ordinal multiplicativity.
    n=len(a)
    if n==0: return 1
    incoming=[sum(1<<i for i in range(n) if a[i][j]) for j in range(n)]
    dp=[[0]*n for _ in range(1<<n)]
    for j in range(n): dp[1<<j][j]=1
    for mask in range(1,1<<n):
        rest=mask
        while rest:
            bit=rest&-rest; rest-=bit; j=bit.bit_length()-1
            before=mask^bit
            if not before: continue
            opts=before&incoming[j]
            val=0
            while opts:
                kbit=opts&-opts; opts-=kbit
                val+=dp[before][kbit.bit_length()-1]
            dp[mask][j]=val
    return sum(dp[-1])

def circulant(q,S):
    return [[i!=j and (j-i)%q in S for j in range(q)] for i in range(q)]

def parity_cycle(q,S):
    # Forward affine numerator in T^q(x)=(3^k*x+c)/2^q.
    c,k=0,0
    for j in range(q):
        if j in S:
            c=3*c+(1<<j)
            k+=1
    D=(1<<q)-3**k
    x=F(c,D)
    n=x
    vals=[]
    for j in range(q):
        vals.append(n)
        require(n.denominator%2==1, 'n.denominator%2==1')
        require((n.numerator%2)==int(j in S), '(n.numerator%2)==int(j in S)')
        n=(3*n+1)/2 if n.numerator%2 else n/2
    require(n==x, 'n==x')
    return x,D,vals

def mm(a,b): return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(len(b))) for j in range(len(b[0]))) for i in range(len(a)))

def mpow(a,n):
    out=tuple(tuple(int(i==j) for j in range(len(a))) for i in range(len(a)))
    for _ in range(n): out=mm(out,a)
    return out

results={}
vals=sorted({F(a,b) for a in range(-8,9) if a for b in range(1,9)})
period_counts=Counter()
for a,b in product(vals,repeat=2):
    p=len(orbit((a,b)))
    expected=1 if (a,b)==(1,1) else 3 if a*a==b*b==1 else 6
    require(p==expected, 'p==expected')
    period_counts[p]+=1
results['rational_period_grid']={'numerators':[-8,8],'denominators':[1,8],'distinct_values':len(vals),'pairs':len(vals)**2,'periods':dict(period_counts)}

integral=[]
for C in range(1,201):
    states=integral_states(C)
    expected=prod(3*c*(c+1)+1 for c in factor(C).values())
    require(len(states)==expected, 'len(states)==expected')
    require(sum(len(set(tuple(s[i:]+s[:i]) for i in range(6)))==1 for s in states)==1, 'sum(len(set(tuple(s[i:]+s[:i]) for i in range(6)))==1 for s in states)==1')
    require(all(all(v%2 for v in s) for s in states)==bool(C%2), 'all(all(v%2 for v in s) for s in states)==bool(C%2)')
    integral.append({'C':C,'states':len(states),'nontrivial_cycles':(len(states)-1)//6})
results['integral_classification']={'range':[1,200],'all_passed':True,'selected':[r for r in integral if r['C'] in (1,5,17,25,85,125,200)]}
results['C85_cycles']=sorted({min(tuple(s[i:]+s[:i]) for i in range(6)) for s in integral_states(85)})

hex_receipts=[]
for m in range(1,13):
    states=hexagon(m)
    require(len(states)==3*m*(m+1)+1, 'len(states)==3*m*(m+1)+1')
    for p in states:
        seen=[]; q=p
        while q not in seen:
            require(q in states, 'q in states')
            order,cut=marked_order(q,m)
            require(len(order)==len(set(order))==4*m, 'len(order)==len(set(order))==4*m')
            seen.append(q); q=step(q,m)
        require(q==p, 'q==p')
        require(len(seen)==(1 if p==(m,m) else 6), 'len(seen)==(1 if p==(m,m) else 6)')
    hex_receipts.append({'radius':m,'states':len(states),'nontrivial_cycles':(len(states)-1)//6})
results['autonomous_marked_hexagons']=hex_receipts

A=((0,1),(-1,1))
I=((1,0),(0,1)); minusI=((-1,0),(0,-1))
require(mpow(A,3)==minusI and mpow(A,6)==I, 'mpow(A,3)==minusI and mpow(A,6)==I')
v=((1,),(0,)); columns=[]
for i in range(6):
    columns.append(v); v=mm(A,v)
B=tuple(tuple(columns[j][i][0] for j in range(6)) for i in range(2))
P=tuple(tuple(int(i==(j+1)%6) for j in range(6)) for i in range(6))
require(mm(B,P)==mm(A,B), 'mm(B,P)==mm(A,B)')
results['sharp_positive_six_lift']={'A':A,'quotient_B':B,'permutation_P':P,'BP_equals_AB':True}

c=F(-29,16)
f=lambda x:x*x+c
preperiodic=[F(j,4) for j in (-7,-5,-3,-1,1,3,5,7)]
require(all(f(x) in preperiodic for x in preperiodic), 'all(f(x) in preperiodic for x in preperiodic)')
quadratic=cycle_from(f,F(-7,4))
require(quadratic==[F(-7,4),F(5,4),F(-1,4)], 'quadratic==[F(-7,4),F(5,4),F(-1,4)]')
results['quadratic_parameter']={'c':str(c),'cycle':list(map(str,quadratic)),'full_rational_preperiodic_graph_proved_in_note':{str(x):str(f(x)) for x in preperiodic}}
clocks=[]
for start in (1,5,17):
    r,t,o=cycle_from(raw,start),cycle_from(terras,start),cycle_from(odd_step,start)
    require([n for n in r if n%2]==o, '[n for n in r if n%2]==o')
    require(len(t)==len(r)-len(o), 'len(t)==len(r)-len(o)')
    clocks.append({'minimum':start,'raw':r,'terras':t,'odd_only':o,'periods':[len(r),len(t),len(o)]})
require([x['periods'] for x in clocks]==[[2,1,1],[5,3,2],[18,11,7]], "[x['periods'] for x in clocks]==[[2,1,1],[5,3,2],[18,11,7]]")
require(sum(n%2 for n in clocks[-1]['raw'])==7, "sum(n%2 for n in clocks[-1]['raw'])==7")
results['minus_collatz_clocks']=clocks

small=[]
for q in (3,5,7):
    k=(q-1)//2
    for signs in product((0,1),repeat=k):
        S={j if signs[j-1] else q-j for j in range(1,k+1)}
        x,D,points=parity_cycle(q,S)
        require(D>0 and all(t>0 for t in points), 'D>0 and all(t>0 for t in points)')
        small.append({'q':q,'S':sorted(S),'H':hamilton(circulant(q,S)),'D':D,'x':str(x)})
results['circulant_collatz_dictionary']=small
results['circulant_H_by_q']={str(q):sorted({r['H'] for r in small if r['q']==q}) for q in (3,5,7)}
require(len(results['circulant_H_by_q']['7'])>1, "len(results['circulant_H_by_q']['7'])>1")

bookA=adjacency_cyclic_book(2); bookB=adjacency_cyclic_book(4)
require(hamilton(bookA)==5 and hamilton(bookB)==17, 'hamilton(bookA)==5 and hamilton(bookB)==17')
phase=(2,1); actual=[]
for i in range(6):
    order,cut=marked_order(phase,1)
    blocks=[bookA if kind=='A' else bookB for kind,index in order]
    prefix=hamilton(ordinal(blocks[:cut]))
    whole=hamilton(ordinal(blocks))
    require(prefix==5**phase[0]*17**phase[1], 'prefix==5**phase[0]*17**phase[1]')
    require(whole==5**2*17**2==7225, 'whole==5**2*17**2==7225')
    actual.append({'phase_coordinates':phase,'order':[kind+str(index+1) for kind,index in order],'cut':cut,'prefix_vertices':sum(map(len,blocks[:cut])),'prefix_H':prefix,'full_H':whole})
    phase=step(phase,1)
require(phase==(2,1), 'phase==(2,1)')
require([a['prefix_H'] for a in actual]==[425,1445,289,17,5,25], "[a['prefix_H'] for a in actual]==[425,1445,289,17,5,25]")
results['actual_tournament_carrier']={'A':{'name':'C3[1,1,2]','n':4,'H':5},'B':{'name':'C3[1,1,4]','n':6,'H':17},'whole_vertices':20,'phases':actual,'independent_full_subset_DP':True}

OUT.write_text(json.dumps(results,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':'all checks passed','receipt':str(OUT),'rational_pairs':len(vals)**2,'integral_C_checked':200,'marked_radii_checked':12,'circulant_H':results['circulant_H_by_q'],'target':7225,'prefixes':[x['prefix_H'] for x in actual]},indent=2))
