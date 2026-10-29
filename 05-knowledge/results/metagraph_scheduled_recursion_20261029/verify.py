"""Same-agent exact tests with distinct partition and labelled-orbit oracles."""
import ast
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import importlib.util
from itertools import permutations,combinations
import json
from math import factorial,gcd
from pathlib import Path
import sys

from scheduled_counts import count_prefix
from precision_counts import normalized_low_bits,base_valuation,cycle_floor,term,v2

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def need(ok,message):
    if not ok:
        raise ArithmeticError(message)


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parts(n,lower=1):
    if n==0:
        yield ()
    for a in range(lower,n+1,2):
        for rest in parts(n-a,a):
            yield (a,)+rest


def literal_roots(n):
    edges=list(combinations(range(n),2))
    index={e:i for i,e in enumerate(edges)}
    perms=list(permutations(range(n)))
    def move(code,p):
        result=0
        for i,(a,b) in enumerate(edges):
            x,y=p[a],p[b]
            bit=((code>>i)&1) ^ (x>y)
            result|=bit << index[tuple(sorted((x,y)))]
        return result
    def canon(code,rooted=False):
        return min(move(code,p) for p in perms if not rooted or p[0]==0)
    rooted={c:canon(c,True) for c in range(1<<len(edges))}
    reps=sorted(set(rooted.values()))
    ids={r:i for i,r in enumerate(reps)}
    Q=[]
    for c in reps:
        row=Counter()
        for v in range(n):
            p=list(range(n));p[0],p[v]=p[v],p[0]
            row[ids[rooted[move(c,p)]]]+=1
        need(sum(row.values())==n and all(v%2 for v in row.values()),'literal odd root orbit weights')
        Q.append(row)
    classes=len({canon(c) for c in reps})
    for i,row in enumerate(Q):
        product=Counter()
        for j,a in row.items():
            for k,b in Q[j].items():
                product[k]+=a*b
        need(dict(product)=={j:n*a for j,a in row.items()},'literal reset Q squared equals nQ')
    trace=sum(row.get(i,0) for i,row in enumerate(Q))
    need(trace==n*classes,'literal root trace class count')
    if n>=3 and n%2:
        need(v2((trace>>((n-1)//2))-1)==1,'new sharp odd root-trace bit')
    return dict(n=n,classes=classes,pointed=len(reps),trace=trace)


def run():
    baseline=load('baseline87','05-knowledge/results/session105_metagraph_recursion_20261029/enumeration/baseline87.py')
    channels=load('channels89','05-knowledge/results/session89_consecutive_tournament_pairs_20261028/counting/channels.py')
    cs=((1,0),(1,3),(2,0),(2,1),(3,2),(4,1),(1,1))
    values,_=channels.fused_prefix(80,cs)
    checks=0
    for (r,s),prefix in zip(cs,values):
        for n,value in enumerate(prefix):
            base=base_valuation(n,r,s)
            need(value%(1<<base)==0,'inherited integral normalization')
            normalized=value>>base
            for bits in (1,2,3,5,8,12):
                expected=normalized% (1<<bits)
                need(normalized_low_bits(n,bits,r,s,False)[0]==expected,'precision truncation/full response comparison')
                need(normalized_low_bits(n,bits,r,s,True)[0]==expected,'sharp shortcut comparison')
                checks+=1
            if n>=3 and n%2:
                need(v2(n*normalized-1)==2*r+2*s-1,'sharp first correction for every odd channel')
    term_checks=0
    for n in range(1,29):
        sums=[Fraction() for _ in cs]
        for lam in parts(n):
            k=len(lam)
            energy=sum((a-1)//2 for a in lam)+sum(gcd(a,b) for i,a in enumerate(lam) for b in lam[i+1:])
            z=1
            for a,m in Counter(lam).items():z*=a**m*factorial(m)
            for c,(r,s) in enumerate(cs):
                sums[c]+=Fraction(1 << (r*energy+s*k),z)
                if n!=2:
                    exponent,odd=term(lam,r,s,bool(n%2))
                    need(exponent==r*energy+s*k-v2(z)-base_valuation(n,r,s),'normalized term formula')
                    if k>2:
                        need(exponent>=cycle_floor(k,r,s,bool(n%2)),'cycle-number cutoff bound')
                    term_checks+=1
        need(all(total.denominator==1 and total.numerator==values[c][n] for c,total in enumerate(sums)),
             'independent Fraction odd-partition sums')
    scheduler=[]
    for bound in (0,1,2,3,10,30,60,100,150,200):
        actual,stats=count_prefix(bound)
        expected,old=baseline.degree_bucket_dp(bound)
        need(actual==expected,'whole-prefix scheduled/reference equality')
        scheduler.append(dict(bound=bound,transitions=stats['transitions'],reference_transitions=old['transitions'],
                              hash=sha256(','.join(hex(v) for v in actual).encode()).hexdigest()))
    for color in (2,3,4):
        need(count_prefix(40,color)[0]==baseline.degree_bucket_dp(40,color)[0],'literal-color whole-prefix equality')
    literal=[literal_roots(n) for n in range(1,6)]
    bank_path=ROOT/'05-knowledge/results/consecutive_flags_order10_20261028/generation.json'
    bank=json.loads(bank_path.read_text())['results']
    rooted_index=cs.index((1,1))
    for row in bank:
        n=row['order']
        need(row['classes']==values[0][n] and row['pointed_states']==values[rooted_index][n-1],
             'inherited class bank and pointed state counts')
        if n>=4 and n%2==0:
            need(v2((n-1)*(row['pointed_states']>>(n//2))-1)==3,'new root-bank second constraint')
    native=load('native_previous','05-knowledge/results/strong_attachment_tail_compiler_20261028/native.py')
    c3=(2,4,1);A=tuple(x|8 for x in c3)+(0,)
    pair=[native.core.substitute(c3,children)[0] for children in ([A]*3,[native.core.converse(A),A,A])]
    groups=[native.oracle.isomorphisms(t,t,True) for t in pair]
    orders=list(map(len,groups))
    transport=[sum(g[3]==7 for g in group) for group in groups]
    cyclic=[{c for c in combinations(range(12),3) if all(sum(t[a]>>b&1 for b in c if a!=b)==1 for a in c)} for t in pair]
    need(orders==[81,27] and transport==[27,0] and cyclic[0]==cyclic[1] and len(cyclic[0])==67,
         '81/27 keeps its joint transporter obstruction')
    large=[]
    for n,bits,r,s in ((1001,5,1,0),(1001,6,2,1),(1000001,2,1,0),(1000001,6,2,1)):
        residue,stats=normalized_low_bits(n,bits,r,s)
        if n==1001 and (r,s)==(2,1):
            need(residue==normalized_low_bits(n,bits,r,s,False)[0],'large sharp formula/full three-cycle check')
        large.append(dict(n=n,bits=bits,r=r,s=s,residue=residue,stats=stats))
    legacy_path=ROOT/'04-computation/fast_a000568_v2_s90di.py'
    tree=ast.parse(legacy_path.read_text(encoding='utf-8'))
    nodes=[node for node in tree.body if isinstance(node,(ast.Import,ast.ImportFrom,ast.FunctionDef))]
    namespace={}
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(legacy_path),'exec'),namespace)
    old_small=[namespace['T_fast'](n) for n in (2,3)]
    need(old_small==[2,4],'historical even-cycle filter failure witness')
    need([namespace['T_faster'](n) for n in (2,3)]==old_small,
         'both historical routines retain the same wrong permutation universe')
    sources=[ROOT/'05-knowledge/results/session105_metagraph_recursion_20261029/enumeration/baseline87.py',
             ROOT/'05-knowledge/results/session89_consecutive_tournament_pairs_20261028/counting/channels.py',bank_path]
    return dict(precision_comparisons=checks,term_checks=term_checks,channels=cs,
                scheduled=scheduler,literal_root_resets=literal,bank_orders_checked=len(bank),
                hostile=dict(orders=orders,transport=transport,cyclic_triples=67),large_partial_outputs=large,
                historical_wrong_counts=old_small,sources={str(p.relative_to(ROOT)).replace('\\','/'):sha256(p.read_bytes()).hexdigest() for p in sources})


if __name__=='__main__':
    result=json.loads(json.dumps(run()))
    target=HERE/'verification.json'
    if '--write' in sys.argv:
        target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    else:
        need(result==json.loads(target.read_text(encoding='utf-8')),'receipt replay')
    print('PASS precision recurrence, sharp correction, scheduled prefixes, literal roots, inherited bank and81/27')
