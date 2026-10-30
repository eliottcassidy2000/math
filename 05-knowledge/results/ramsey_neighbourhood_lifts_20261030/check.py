"""Exact finite controls for the all-size proofs; explicit gates survive -O."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations, product
from pathlib import Path
import argparse
import importlib.util
import json
from math import gcd
import random
from carousel import chart, check

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
GATES=0
def need(ok,label):
    global GATES
    GATES+=1
    if not ok: raise RuntimeError(label)

@lru_cache(None)
def g(k):
    if k<2:return k
    h=k//2
    return 3*g(h) if k%2==0 else 2*g(h)+g(h+1)

@lru_cache(None)
def w(k):
    if k<2:return F(k)
    h=k//2
    return F(5,2)*w(h) if k%2==0 else F(3,2)*w(h)+w(h+1)

def rows_from_bits(n,mask):
    rows=[0]*n
    for t,(i,j) in enumerate(combinations(range(n),2)):
        if mask>>t&1:rows[i]|=1<<j
        else:rows[j]|=1<<i
    return rows

def induced_bits(rows,vs):
    return sum(((rows[i]>>j)&1)<<t for t,(i,j) in enumerate(combinations(vs,2)))

def tau(rows):
    @lru_cache(None)
    def rec(mask):
        if not mask:return 0
        return 1+max(rec(mask&rows[i]) for i in range(len(rows)) if mask>>i&1)
    return rec((1<<len(rows))-1)

def carousel_controls():
    hist={};forced_tests=0
    for n in range(3,16,2):
        for mask in range(1<<n):
            for eps in (0,1):check(n,mask,eps,chart(n,eps,mask))
            # Separate geometric test uses actual intervals, not the construction.
            marked=[i for i in range(n) if mask>>i&1]
            if len(marked)>=4:
                singles=set()
                for s in range(n):
                    inside=[i for i in marked if (i-s)%n<=n//2]
                    if len(inside)==1:singles.add(marked.index(inside[0]))
                need(len(singles)<=2,'at most two singleton-visible marks')
                if len(singles)==2:
                    a,b=sorted(singles);need((a-b)%len(marked) in (1,len(marked)-1),'forced marks adjacent')
                forced_tests+=1
        hist[str(n)]=2**(n+1)
    rng=random.Random(20261030);weighted=0
    for n in (17,19,31,65,127):
        for _ in range(250):
            v=[rng.randrange(41) for i in range(n)]
            k=max(sum(v[(s+i)%n] for i in range(n//2+1)) for s in range(n))
            odd=sum((a%2)<<i for i,a in enumerate(v));eps=k%2
            c=chart(n,eps,odd);check(n,odd,eps,c)
            need(sum(g(a) for a in v)<=g(k),'all-weight carousel inequality')
            layers=[[a//2+int(a%2==1 and c[i]==j) for i,a in enumerate(v)] for j in range(3)]
            need(sum(g(a) for a in v)==sum(g(a) for layer in layers for a in layer),'exact mass split')
            weighted+=1
    return {'complete_charts':hist,'forced_geometry_cases':forced_tests,'larger_weighted_queries':weighted}

def tournament_controls():
    u=[sum(1<<j for j in range(5) if 0<(j-i)%5<=2) for i in range(5)]
    u[0]^=1<<1;u[1]^=1
    bad={induced_bits(u,p) for p in permutations(range(5))}
    histogram={};total=0
    for n in range(1,7):
        counts=Counter();subsets=list(combinations(range(n),5))
        for bits in range(1<<(n*(n-1)//2)):
            rows=rows_from_bits(n,bits)
            if any(induced_bits(rows,s) in bad for s in subsets):continue
            t=tau(rows);need(n<=g(t),'native U5-free capacity');counts[t]+=1;total+=1
        histogram[str(n)]={str(k):v for k,v in sorted(counts.items())}
    paley=[sum(1<<j for j in range(7) if (j-i)%7 in (1,2,4)) for i in range(7)]
    witnesses=[s for s in combinations(range(7),5) if induced_bits(paley,s) in bad]
    need(tau(paley)==3 and 7>g(3) and witnesses,'Paley7 outside U5-free scope')
    for k in range(1,17):
        cells=[(j,t) for j in range(k) for t in range(j+1) if t&j==t];rows=[0]*len(cells)
        for i,j in combinations(range(len(cells)),2):
            a,b=cells[i];c,d=cells[j]
            while (a&1,b&1)==(c&1,d&1):a>>=1;b>>=1;c>>=1;d>>=1
            digit={(0,0):0,(1,0):1,(1,1):2}
            if (digit[c&1,d&1]-digit[a&1,b&1])%3==1:rows[i]|=1<<j
            else:rows[j]|=1<<i
        need(len(rows)==g(k) and tau(rows)==k,'literal Pascal attainer')
    spec=importlib.util.spec_from_file_location('hostile_native',ROOT/'05-knowledge/results/regular_parent_symmetry_20261029/verify.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    hostile=module.hostile_checks()
    return {'labelled_U5_free_through_6':total,'tau_histograms':histogram,'Paley7_U5_witnesses':witnesses,'Pascal_budgets':16,'hostile_81_27':hostile}

def formula(a,b):
    if a==1:return b
    if a==2:return 5*b//2
    if a==3:return 4*b-int(b%4!=0)
    if a==4:return 25*(b//4)+(0,4,10,16)[b%4]
    raise ValueError(a)

@lru_cache(None)
def budgets(k,step):
    return tuple(v for v in product(range(1,k),repeat=5) if all(v[i]+v[(i+step)%5]<=k for i in range(5)))

def graph_controls():
    # Independent exact decomposition DP: no use of the proposed row formulas.
    table={};plans={};candidates=0
    for a in range(1,5):
        for b in range(1,9):
            best=1;plan=('vertex',)
            for x in range(1,a):
                val=table[x,b]+table[a-x,b]
                if val>best:best=val;plan=('join',x)
            for y in range(1,b):
                val=table[a,y]+table[a,b-y]
                if val>best:best=val;plan=('union',y)
            for aa in budgets(a,1):
                for bb in budgets(b,2):
                    val=sum(table[x,y] for x,y in zip(aa,bb));candidates+=1
                    if val>best:best=val;plan=('C5',aa,bb)
            table[a,b]=best;plans[a,b]=plan
            need(best==formula(a,b),'four exact graph capacity rows')
            need(best<=min(w(a)*b,w(b)*a),'dual recursive bound')
    # Native graph witnesses from the independent optimal DP plans.
    def sub(q,children):
        out=[];offset=0;starts=[]
        for child in children:starts.append(offset);offset+=len(child)
        for i,child in enumerate(children):
            for row in child:
                v=row<<starts[i]
                for j,other in enumerate(children):
                    if q[i][j]:v|=((1<<len(other))-1)<<starts[j]
                out.append(v)
        return out
    @lru_cache(None)
    def build(a,b):
        plan=plans[a,b]
        if plan[0]=='vertex':return (0,)
        if plan[0]=='join':children=[build(plan[1],b),build(a-plan[1],b)];q=[[0,1],[1,0]]
        elif plan[0]=='union':children=[build(a,plan[1]),build(a,b-plan[1])];q=[[0,0],[0,0]]
        else:children=[build(x,y) for x,y in zip(plan[1],plan[2])];q=[[int((i-j)%5 in (1,4)) for j in range(5)] for i in range(5)]
        return tuple(sub(q,children))
    def clique(rows):
        @lru_cache(None)
        def rec(mask):
            if not mask:return 0
            bit=mask&-mask;i=bit.bit_length()-1
            return max(rec(mask^bit),1+rec(mask&rows[i]))
        return rec((1<<len(rows))-1)
    for (a,b),n in table.items():
        rows=build(a,b);full=(1<<n)-1
        need(len(rows)==n and clique(rows)<=a and clique([full^row^(1<<i) for i,row in enumerate(rows)])<=b,'native graph realization')
    pack=0
    for k in range(1,13):
        for v in product(range(k+1),repeat=5):
            if all(v[i]+v[(i+1)%5]<=k for i in range(5)):
                need(sum(w(x) for x in v)<=2*w(k),'C5 fractional-color packing');pack+=1
    for a in range(1,65):
        for b in range(1,65):need(w(a)+w(b)<=w(a+b),'K2 superadditivity')
    return {'capacities_A1_to_4_B1_to_8':[[table[a,b] for b in range(1,9)] for a in range(1,5)],'exact_DP_candidates':candidates,'native_graphs':len(table),'packing_vectors':pack,'slopes_1_to_32':[str(w(i)) for i in range(1,33)]}

def neighbourhood_controls():
    bank=json.loads((ROOT/'05-knowledge/results/ramsey_saturated_ports_20261030/verification.json').read_text())['banks']['r36_17.g6']['entries']
    edges=[row['edges'] for row in bank];need(edges==[41,41,40,41,40,42,42],'critical17 short-proof hostile')
    vals=[]
    for h in range(35):
        d=34-h;p=6 if h==34 else 3*(h%2)
        value=-12+40*d+3*d*d+p
        need(value+6>=40*d,'deficit stability inequality');vals.append(value)
    return {'critical17_edges':edges,'F_h0_to_34':vals,'hypothetical_n87_total_deficits':[2,4,6,8,10,12],'minimum_saturated_vertices':75,'bound_improvement_claimed':False}

def incoming_spectral_controls():
    # Literal two witnesses from Session120's complete N=231 census.
    n=231;u=[sum(1<<j for j in range(5) if 0<(j-i)%5<=2) for i in range(5)]
    u[0]^=2;u[1]^=1;bad={induced_bits(u,p) for p in permutations(range(5))};out=[]
    for nonlinear,expected in ((False,17225670),(True,29082438)):
        symbol=[0]
        for x in range(1,n):
            m=n//gcd(n,x);t=x//(n//m)
            s={p:(1 if pow(t%p,(p-1)//2,p)==1 else -1) for p in (3,7,11) if m%p==0}
            if m in (3,7,11):v=-s[m]
            elif m==21:v=s[3] if nonlinear else s[7]
            elif m==33:v=-s[11] if nonlinear else s[3]
            elif m==77:v=-s[11] if nonlinear else s[11]
            elif nonlinear:v=(s[3]+s[7]-s[11]+s[3]*s[7]*s[11])//2
            else:v=-s[3]*s[7]*s[11]
            symbol.append(v)
        need(all(symbol[-x]==-symbol[x] for x in range(1,n)),'actual odd symbol')
        moment=n*sum(sum(symbol[x]*symbol[(x+h)%n] for x in range(n))**2 for h in range(n))
        need(moment==expected,'incoming quartic optimum replay')
        rows=[sum(1<<j for j in range(n) if symbol[(j-i)%n]==1) for i in range(n)]
        witness=next((s for s in combinations(range(12),5) if induced_bits(rows,s) in bad),None)
        need(witness is not None,'spectral minimizer does not pass U5-free guard')
        out.append({'nonlinear':nonlinear,'fourth_moment':moment,'U5_witness':witness})
    return out

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    out={}
    for name,fn in [('carousel',carousel_controls),('tournaments',tournament_controls),('graphs',graph_controls),('neighbourhoods',neighbourhood_controls),('incoming_spectral',incoming_spectral_controls)]:
        out[name]=fn();print(name,'passed',flush=True)
    out['gates']=GATES
    # Canonicalize tuples as JSON lists before comparing a saved receipt.
    out=json.loads(json.dumps(out))
    if args.write:(HERE/'verification.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    elif (HERE/'verification.json').exists():need(out==json.loads((HERE/'verification.json').read_text()),'published receipt replay')
    print(json.dumps({'gates':GATES,'status':'PASS'}))
if __name__=='__main__':main()
