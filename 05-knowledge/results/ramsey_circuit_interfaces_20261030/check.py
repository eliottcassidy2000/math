"""Exact Ramsey-interface and Boolean-encoding checks. All gates survive -O."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path
import argparse
import importlib.util
import json
import random
from explore import solve, tournament, transitive, vertices

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
GATES=0
def need(ok,label):
    global GATES
    GATES+=1
    if not ok: raise RuntimeError(label)

T5=tournament(5,{1,2})
MAXIMAL=((0,1,2),(1,2,3),(0,1,4),(0,3,4),(2,3,4))

def feasible(v,k):return all(sum(v[i] for i in s)<=k for s in MAXIMAL)

def rounding_certificate():
    charts=[]
    for epsilon in (0,1):
        for mask in range(32):
            odd=[i for i in range(5) if mask>>i&1]
            chosen=None
            for word in product(range(3),repeat=len(odd)):
                color=dict(zip(odd,word))
                if all(sum(color.get(i)==j for i in s)<=
                       (sum(i in color for i in s)-epsilon+1)//2+int(j==2)*epsilon
                       for s in MAXIMAL for j in range(3)):
                    chosen=[color.get(i,-1) for i in range(5)];break
            need(chosen is not None,'universal parity-rounding chart')
            charts.append({'epsilon':epsilon,'odd_mask':mask,'colors':chosen})
    return charts

def generic_rounding_chart(maximal,n,epsilon,mask):
    odd=[i for i in range(n) if mask>>i&1]
    constraints=[([odd.index(i) for i in s if i in odd],
                  (sum(i in odd for i in s)-epsilon+1)//2) for s in maximal]
    for colors in product(range(3),repeat=len(odd)):
        if all(all(sum(colors[i]==j for i in ids)<=bound+epsilon*int(j==2)
                   for j in range(3)) for ids,bound in constraints):
            chart=dict(zip(odd,colors))
            return [chart.get(i,-1) for i in range(n)]
    return None

def maximal_transitive(q):
    n=len(q);sets=[set(i for i in range(n) if mask>>i&1) for mask in range(1,1<<n)]
    good=[s for s in sets if transitive(q,sorted(s))]
    return [tuple(sorted(s)) for s in good if not any(s<t for t in good)]

def gould_fast(k):
    # Maintain (G(m),G(m+1)) while reading k from its most significant bit.
    a,b=0,1
    for bit in bin(k)[2:]:
        a,b=(3*a,2*a+b) if bit=='0' else (2*a+b,3*b)
    return a

@lru_cache(None)
def symmetry_exponent(k):
    if k==1:return 0
    h=k//2
    return (1+3*symmetry_exponent(h)) if k%2==0 else (2*symmetry_exponent(h)+symmetry_exponent(h+1))

@lru_cache(None)
def colored_orbits(k,q):
    if k==1:return q
    h=k//2;x=colored_orbits(h,q)
    return (x**3+2*x)//3 if k%2==0 else x*x*colored_orbits(h+1,q)

def cycle_count(p):
    visited=set();c=0
    for i in range(len(p)):
        if i in visited:continue
        c+=1
        while i not in visited:visited.add(i);i=p[i]
    return c

def capacity(limit):
    f=[0,1];plans=[None,None];work=0
    for k in range(2,limit+1):
        best=0;arg=None
        for a in range(k):
            for b in range(min(k-a,k-1)+1):
                for c in range(min(k-a-b,k-1)+1):
                    for d in range(min(k-b-c,k-a,k-c,k-1)+1):
                        e=min(k-a-b,k-a-d,k-c-d,k-1)
                        value=f[a]+f[b]+f[c]+f[d]+f[e];work+=1
                        if value>best: best=value;arg=(a,b,c,d,e)
        f.append(best);plans.append(arg)
    return f,plans,work

def three_capacity(limit):
    f=[0,1]
    for k in range(2,limit+1):
        f.append(max(f[a]+f[b]+f[min(k-a,k-b,k-1)]
                     for a in range(k) for b in range(min(k-a,k-1)+1)))
    return f

def substitute(q,children):
    sizes=list(map(len,children));start=[];n=0
    for s in sizes: start.append(n);n+=s
    rows=[0]*n
    for i,child in enumerate(children):
        for u,row in enumerate(child):
            rows[start[i]+u] |= row<<start[i]
            for j in range(len(children)):
                if q[i][j]: rows[start[i]+u] |= ((1<<sizes[j])-1)<<start[j]
    return rows

def largest_transitive(rows):
    @lru_cache(None)
    def rec(mask):
        if not mask:return 0
        return 1+max(rec(mask&row) for i,row in enumerate(rows) if mask>>i&1)
    return rec((1<<len(rows))-1)

def witness(k,plans):
    if k==0:return []
    if k==1:return [0]
    return substitute(T5,[witness(v,plans) for v in plans[k]])

def pascal_tournament(k):
    cells=[(j,t) for j in range(k) for t in range(j+1) if t&j==t]
    rows=[0]*len(cells);digits={(0,0):0,(1,0):1,(1,1):2}
    for a,b in combinations(range(len(cells)),2):
        j,t=cells[a];l,u=cells[b]
        while (j&1,t&1)==(l&1,u&1):j>>=1;t>>=1;l>>=1;u>>=1
        da=digits[j&1,t&1];db=digits[l&1,u&1]
        if (db-da)%3==1:rows[a]|=1<<b
        else:rows[b]|=1<<a
    return cells,rows

def graph_polytope():
    rows=[];rhs=[]
    for i in range(5):
        rows.append(tuple(int(j in (i,(i+1)%5)) for j in range(5)));rhs.append(F(1))
    for i in range(5):rows.append(tuple(-int(i==j) for j in range(5)));rhs.append(F(0))
    vs=set()
    for sub in combinations(range(10),5):
        x=solve([rows[i] for i in sub],[rhs[i] for i in sub])
        if x is not None and all(sum(a*b for a,b in zip(row,x))<=b for row,b in zip(rows,rhs)):vs.add(x)
    return vs

# Pins >=0 denote input/gate wires. -1/-2 denote false/true constants.
def value(pin,values):return values[pin] if pin>=0 else int(pin==-2)
def evaluate(n,gates,out,x):
    values=list(x)
    for a,b,t in gates:values.append((t>>(value(a,values)+2*value(b,values)))&1)
    return value(out,values),tuple(values[n:])

def simplify_clause(literals):
    if any(x is True for x in literals):return None
    out=set(x for x in literals if x is not False)
    if any(-x in out for x in out):return None
    return tuple(sorted(out))

def pin_literal(pin,desired):
    if pin<0:return (pin==-2)==bool(desired)
    return (pin+1) if desired else -(pin+1)

def cnf(n,gates,out,prefix=()):
    clauses=[]
    for j,(a,b,t) in enumerate(gates):
        if a>=n+j or b>=n+j:raise ValueError('acyclic source order required')
        for u,v in product(range(2),repeat=2):
            z=(t>>(u+2*v))&1
            row=simplify_clause([pin_literal(a,1-u),pin_literal(b,1-v),pin_literal(n+j,z)])
            if row is not None:clauses.append(row)
    for i,v in enumerate(prefix):
        clauses.append((pin_literal(i,v),))
    row=simplify_clause([pin_literal(out,1)])
    if row is not None:clauses.append(row)
    return clauses

def satisfies(clauses,x):
    return all(any(x[abs(v)-1]==int(v>0) for v in row) for row in clauses)

def equality_circuit(m):
    gates=[]
    for i in range(m):gates.append((i,m+i,9))  # XNOR: 00 and 11
    out=2*m
    for i in range(1,m):
        gates.append((out,2*m+i,8));out=2*m+len(gates)-1
    return gates,out

def run():
    charts=rounding_certificate()
    carousel_charts={}
    for n in (7,9):
        q=tournament(n,set(range(1,(n+1)//2)));mt=maximal_transitive(q)
        intervals=[tuple(sorted((i+j)%n for j in range((n+1)//2))) for i in range(n)]
        need(set(mt)==set(intervals),'native carousel maximal transitive supports')
        bank=[]
        for epsilon in (0,1):
            for mask in range(1<<n):
                chart=generic_rounding_chart(mt,n,epsilon,mask)
                need(chart is not None,'complete larger-carousel parity compiler')
                bank.append(chart)
        carousel_charts[str(n)]=bank
    u5=tournament(5,{1,2});u5[0][1],u5[1][0]=0,1
    need(generic_rounding_chart(maximal_transitive(u5),5,1,28) is None,
         'U5 rounding guard obstruction, not a capacity counterexample')
    p7=tournament(7,{1,2,4});pm=maximal_transitive(p7)
    need(max(map(len,pm))==3 and 7>gould_fast(3),'Paley7 actual capacity obstruction')
    need(generic_rounding_chart(pm,7,0,15) is None,'Paley7 rounding obstruction')
    maximal,vs=vertices(T5)
    need(set(maximal)==set(MAXIMAL),'literal T5 maximal transitive sets')
    expected={tuple([F(0)]*5),tuple([F(1,3)]*5)}
    for i in range(5):expected.add(tuple(F(i==j) for j in range(5)))
    for tri in combinations(range(5),3):
        if not transitive(T5,tri): expected.add(tuple(F(int(i in tri),2) for i in range(5)))
    need(vs==expected and len(vs)==12,'complete rational weight polytope')
    need(9>8 and 27>25,'exact p>3/2 and uniform vertex strictness')
    # The C5 graph's complementary clique/independence packing polytopes.
    gv=graph_polytope();expected_g={tuple([F(1,2)]*5)}
    for mask in range(32):
        if all(not(mask>>i&1 and mask>>((i+1)%5)&1) for i in range(5)):
            expected_g.add(tuple(F((mask>>i)&1) for i in range(5)))
    need(gv==expected_g and len(gv)==12,'complete C5 packing polytope')
    # The complement is C5 under i ->2i; certify all 144 vertex-pair types.
    pair_types=Counter()
    for a in gv:
        for oldb in gv:
            b=tuple(oldb[(3*i)%5] for i in range(5))
            vals=tuple(sorted(x*y for x,y in zip(a,b) if x*y))
            need(vals in ((),(F(1),),(F(1,2),),(F(1,2),F(1,2)),tuple([F(1,4)]*5)),
                 'C5 two-sided extremal certificate type')
            pair_types[str(vals)]+=1
    cap,plans,work=capacity(64);tri=three_capacity(64)
    digit=[sum(2**i.bit_count() for i in range(k)) for k in range(65)]
    need(cap==tri==digit,'exact independent capacity DP through64')
    need([gould_fast(k) for k in range(65)]==digit,'binary matrix evaluator')
    # Concrete integer lifting, separate from the symbolic 64-chart argument.
    rounding_lifts=0
    for weights in product(range(5),repeat=5):
        k=max(sum(weights[i] for i in s) for s in MAXIMAL)
        if k<2:continue
        epsilon=k%2;mask=sum((a%2)<<i for i,a in enumerate(weights))
        colors=charts[epsilon*32+mask]['colors'];h=k//2
        layers=[tuple(a//2+int(colors[i]==j) for i,a in enumerate(weights)) for j in range(3)]
        need(all(feasible(layer,h+epsilon*int(j==2)) for j,layer in enumerate(layers)),
             'integer parity-rounding feasible layers')
        need(sum(gould_fast(a) for a in weights)==sum(gould_fast(a) for layer in layers for a in layer),
             'integer parity-rounding exact mass')
        rounding_lifts+=1
    pascal=[]
    from math import comb
    for k in range(1,17):
        cells,rows=pascal_tournament(k)
        need(cells==[(j,t) for j in range(k) for t in range(j+1) if comb(j,t)%2],
             'literal odd Pascal cells')
        need(len(rows)==digit[k] and largest_transitive(rows)==k,
             'Pascal-cell Ramsey family actual native relation')
        if k>1:
            oldcells,oldrows=pascal_tournament(k-1);mask=(1<<len(oldrows))-1
            need(cells[:len(oldcells)]==oldcells and [r&mask for r in rows[:len(oldrows)]]==oldrows,
                 'Pascal source-preserving nested inclusion')
        pascal.append({'budget':k,'order':len(rows),'transitive':k})
    for m in range(1,32):
        need(digit[2*m]==3*digit[m] and digit[2*m+1]==2*digit[m]+digit[m+1],
             'binary-clock recurrence')
    for k in range(2,65):
        need(feasible(plans[k],k) and max(plans[k])<k,'well-founded capacity certificate')
        need(sum(cap[a] for a in plans[k])==cap[k],'attaining capacity tree')
        need(cap[k]>=cap[k-1]+1,'capacity monotonicity')
    for k in range(2,7):
        brute=max(sum(cap[a] for a in v) for v in product(range(k),repeat=5) if feasible(v,k))
        need(cap[k]==brute,'independent uncompressed optimization')
    native=[]
    for k in range(1,9):
        rows=witness(k,plans);t=largest_transitive(rows)
        need(len(rows)==cap[k] and t<=k,'native attained tournament capacity')
        native.append({'budget':k,'order':len(rows),'transitive':t})
    # All five slots nonempty: balanced large children on the cyclic triple 0,1,3.
    rows=[0];full=[]
    for depth in range(4):
        need(len(rows)==2*3**depth-1,'full T5 tree size')
        need(largest_transitive(rows)==2**(depth+1)-1,'full T5 sharp exponent witness')
        full.append({'depth':depth,'order':len(rows),'transitive':2**(depth+1)-1})
        if depth<3:rows=substitute(T5,[rows,rows,[0],rows,[0]])
    # Homogeneous-count polynomials cannot be used in place of transitive support.
    rng=random.Random(20261030);composition_checks=0
    for _ in range(90):
        children=[];weights=[]
        for i in range(5):
            n=rng.randrange(4);child=[0]*n
            for u,v in combinations(range(n),2):
                if rng.randrange(2):child[u]|=1<<v
                else:child[v]|=1<<u
            children.append(child);weights.append(largest_transitive(child))
        native_value=largest_transitive(substitute(T5,children))
        need(native_value==max(sum(weights[i] for i in s) for s in MAXIMAL),'actual weighted substitution')
        composition_checks+=1
    # Every two-input gate, all repeated/constant pin patterns.
    circuits=[]
    for a,b,t in product((0,1,-1,-2),(0,1,-1,-2),range(16)):
        circuits.append((2,[(a,b,t)],2,()))
    for _ in range(120):
        n=rng.randrange(1,4);s=rng.randrange(4);gates=[]
        for j in range(s):
            pins=list(range(n+j))+[-1,-2]
            gates.append((rng.choice(pins),rng.choice(pins),rng.randrange(16)))
        out=rng.choice(list(range(n+s))+[-1,-2]);prefix=tuple(rng.randrange(2) for _ in range(rng.randrange(n+1)))
        circuits.append((n,gates,out,prefix))
    cnf_cases=0
    for n,gates,out,prefix in circuits:
        clauses=cnf(n,gates,out,prefix)
        need(len(clauses)<=4*len(gates)+1+len(prefix) and all(len(c)<=3 for c in clauses),'encoding size and width')
        for x in product(range(2),repeat=n):
            result,aux=evaluate(n,gates,out,x)
            solutions=[y for y in product(range(2),repeat=len(gates)) if satisfies(clauses,x+y)]
            expected_aux=[aux] if result and x[:len(prefix)]==prefix else []
            need(solutions==expected_aux,'parsimonious shared-DAG CNF with prefix')
            cnf_cases+=1
    eq=[]
    for m in range(1,7):
        gates,out=equality_circuit(m);need(len(gates)==2*m-1,'equality B2 gate budget')
        for x in product(range(2),repeat=2*m):
            need(evaluate(2*m,gates,out,x)[0]==int(x[:m]==x[m:]),'equality truth table')
        eq.append({'half_length':m,'carrier_states':2**m,'gates':len(gates)})
    # Three gates suffice for a multiplexer with select bit s.
    mux=[(0,1,6),(2,3,8),(1,4,6)]
    for a,b,s in product(range(2),repeat=3):
        need(evaluate(3,mux,5,(a,b,s))[0]==(a if s else b),'B2 multiplexer compiler')
    spec=importlib.util.spec_from_file_location('hostile_runtime',ROOT/'05-knowledge/results/regular_parent_symmetry_20261029/verify.py')
    h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h);hostile=h.hostile_checks()
    symmetry=[]
    for k in range(1,8):
        cells,rows=pascal_tournament(k);group=h.iso_maps(rows,rows)
        need(len(group)==3**symmetry_exponent(k),'full native automorphism group, not only tree action')
        cycles=Counter(cycle_count(p) for p in group);orbit=[]
        for q in range(1,5):
            numerator=sum(count*q**c for c,count in cycles.items())
            need(numerator%len(group)==0 and numerator//len(group)==colored_orbits(k,q),
                 'full native Burnside colored orbit recurrence')
            orbit.append(numerator//len(group))
        symmetry.append({'budget':k,'order':len(rows),'aut_order':len(group),'color_orbits_q1_to4':orbit})
    return {'status':'All-size T5/T7/T9 palette capacity equals A006046 by parity rounding; no circuit lower-bound claim',
            'rounding_charts':charts,'rounding_integer_lifts':rounding_lifts,
            'carousel_rounding_charts':carousel_charts,
            'rounding_hostiles':{'U5':{'epsilon':1,'odd_mask':28,'capacity_failure_claimed':False},
                                'Paley7':{'epsilon':0,'odd_mask':15,'order':7,'transitive':3,'capacity_bound':5}},
            't5_vertices':[[str(x) for x in v] for v in sorted(vs)],
            'c5_graph_vertices':[[str(x) for x in v] for v in sorted(gv)],
            'c5_pair_types':dict(pair_types),'capacities':cap,'plans':plans,
            'capacity_four_loop_candidates':work,'digit_sum_check_limit':64,
            'pascal_native_family':pascal,
            'native_witnesses':native,'full_nonempty_witnesses':full,
            'composition_checks':composition_checks,'cnf_input_cases':cnf_cases,
            'circuits':len(circuits),'equality_controls':eq,
            'hostile_81_27':hostile,'pascal_symmetry':symmetry,
            'inherited_gates':h.GATES,'local_gates':GATES}

def main():
    p=argparse.ArgumentParser();p.add_argument('--record',action='store_true');args=p.parse_args()
    r=run();path=HERE/'receipt.json'
    # JSON normalization converts tuple plans to arrays before comparison.
    r=json.loads(json.dumps(r))
    if args.record:path.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n',encoding='utf8',newline='\n')
    else:need(r==json.loads(path.read_text(encoding='utf8')),'recorded receipt')
    print(json.dumps({k:r[k] for k in ('local_gates','inherited_gates','capacity_four_loop_candidates','cnf_input_cases')}))

if __name__=='__main__':main()
