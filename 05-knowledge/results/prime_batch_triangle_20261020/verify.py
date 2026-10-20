#!/usr/bin/env python3
"""Exact finite controls. Checks survive python -O; proofs are in README.md."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, permutations, product
from math import comb, factorial, gcd
from pathlib import Path
import hashlib
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
PREV = HERE.parent / 'continuation_width_symmetry_20261020'
def load_module(name, path):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result

old=load_module('inherited_width_verify',PREV/'verify.py')
saved_verify=sys.modules.get('verify')
sys.modules['verify']=old
try:
    inherited=load_module('inherited_width_controls',PREV/'independent_controls.py')
finally:
    if saved_verify is None: del sys.modules['verify']
    else: sys.modules['verify']=saved_verify

CHECKS = Counter()
def need(value, label):
    CHECKS[label] += 1
    if not value:
        raise RuntimeError(label)


def induced(adj, vertices):
    return [sum(((adj[a] >> b) & 1) << j for j, b in enumerate(vertices))
            for a in vertices]


def canonical(adj):
    # Independent of feedback/group routines: degree blocks + full relabelling.
    groups = defaultdict(list)
    for i, row in enumerate(adj):
        groups[row.bit_count()].append(i)
    groups = [groups[d] for d in sorted(groups)]
    return min(tuple(induced(adj, sum((list(t) for t in order), [])))
               for order in product(*(permutations(g) for g in groups)))


def module(adj, vertices):
    chosen = set(vertices)
    return all(len({(adj[a] >> b) & 1 for a in chosen}) == 1
               for b in range(len(adj)) if b not in chosen)


def triangles(adj):
    return [c for c in combinations(range(len(adj)), 3)
            if all(sum((adj[a] >> b) & 1 for b in c) == 1 for a in c)]


def local_tests(adj, weights):
    ts = triangles(adj)
    common = {(a,b) for a in range(len(adj)) for b in range(len(adj))
              if adj[a] >> b & 1 and weights[a]*weights[b] == 1}
    for tri in ts:
        common &= {(a,b) for a,b in common if a in tri and b in tri}
    rigid = not any(all(weights[a] == 1 for a in tri) and module(adj,tri)
                    for tri in ts)
    return bool(ts) and not common, rigid


def weighted_feedback(adj, weights):
    n = len(adj)
    dp = [0]*(1<<n)
    for mask in range(1,1<<n):
        dp[mask] = min(dp[mask^(1<<v)] + weights[v]*sum(
            weights[w] for w in range(n) if mask>>w&1 and adj[v]>>w&1)
                       for v in range(n) if mask>>v&1)
    return dp[-1]


def sync(sign):
    adj = [0]*6
    for i in range(3):
        adj[i] |= 1<<((i+1)%3)
        adj[3+i] |= 1<<(3+(i+sign)%3)
        for j in range(3):
            if i == j: adj[i] |= 1<<(3+j)
            else: adj[3+j] |= 1<<i
    return adj


def census():
    classes = {()}
    rows = []
    sync_types = {canonical(sync(sign)):sign for sign in (1,-1)}
    need(len(sync_types)==2,'two synchronized types nonisomorphic')
    a = [1]
    for n in range(1,8):
        a.append(a[n-1]+(a[n-3] if n>=3 else 0))
    for n in range(8):
        if n:
            new = set()
            for prev in classes:
                for mask in range(1<<(n-1)):
                    adj = list(prev)+[mask]
                    for i in range(n-1):
                        if not (mask>>i&1): adj[i] |= 1<<(n-1)
                    new.add(canonical(adj))
            classes = new
        need(len(classes)==(1,1,1,2,4,12,56,456)[n],
             'unlabelled extension census matches A000568 through seven')
        slacks = Counter()
        synchronized = 0
        for adj in classes:
            f, order = old.feedback(adj)
            group = old.automorphisms(adj)
            orbits = sorted({tuple(sorted({p[v] for p in group}))
                             for v in range(n)})
            support = sum(len(o) for o in orbits if len(o)>1)
            delta = 3*f-support
            slacks[delta] += 1
            for p in (5,7):
                z, power = len(group), 0
                while z%p==0:
                    z//=p; power+=1
                need(delta>=power*(3*(p*p-1)//8-p),
                     'prime-power symmetry slack lower bound')
            if delta != 6: continue
            need(all(len(o) in (1,3) for o in orbits),
                 'slack six only singleton and triangle orbits')
            nonmodules = [o for o in orbits if len(o)==3 and not module(adj,o)]
            if nonmodules:
                synchronized += 1
                need(len(nonmodules)==2,'exactly two nonmodule triangle orbits')
                vertices = sorted(sum((list(o) for o in nonmodules),[]))
                need(module(adj,vertices),'synchronized union is a module')
                need(canonical(induced(adj,vertices)) in sync_types,
                     'central module one of two synchronized types')
                need(all(module(adj,o) for o in orbits if o not in nonmodules),
                     'all other orbit modules quiet')
            else:
                q = induced(adj,[o[0] for o in orbits])
                weights = list(map(len,orbits))
                need(weighted_feedback(q,weights)==2,
                     'nonsynchronized quotient weighted feedback two')
                need(local_tests(q,weights)==(True,True),
                     'actual nonsynchronized quotient passes local grammar')
                need(len(group)==3**weights.count(3),
                     'nonsynchronized full group product')
        expected = 2*sum(a[i]*a[n-6-i] for i in range(n-5)) if n>=6 else 0
        need(synchronized==expected,'synchronized sector exact generating series')
        rows.append({'n':n,'classes':len(classes),
                     'slack_counts':dict(sorted(slacks.items())),
                     'synchronized':synchronized})
    return rows


def synchronized_ports():
    receipts=[]
    for sign in (1,-1):
        center=sync(sign)
        f, order=old.feedback(center)
        # Independent 720-order check, not the inherited subset DP.
        brute=min(sum(bool(center[p[j]]>>p[i]&1)
                      for i in range(6) for j in range(i+1,6))
                  for p in permutations(range(6)))
        need(brute==f==4,'synchronized feedback by two algorithms')
        receipts.append({'sign':sign,'adjacency':center,'optimal_order':order,
                         'feedback':f,'group':old.automorphisms(center)})
        for total in range(6):
            for left_n in range(total+1):
                for left in inherited.atom_words(left_n):
                    for right in inherited.atom_words(total-left_n):
                        adj=inherited.ordered_join([inherited.atom_word(left),center,
                                                    inherited.atom_word(right)])
                        t=(left+right).count(3)
                        u=(left+right).count(1)
                        group=old.automorphisms(adj)
                        need(old.feedback(adj)[0]==t+4,'generated synchronized feedback')
                        need(len(group)==3**(t+1),'generated synchronized group order')
                        middle=set(range(left_n,left_n+6))
                        triangle_ports=[]
                        start=0
                        for size in left+(6,)+right:
                            if size==3: triangle_ports.append(set(range(start,start+3)))
                            start+=size
                        for ports in [()] + [(i,) for i in range(len(adj))] + list(combinations(range(len(adj)),2)):
                            h=sum(bool(o&set(ports)) for o in triangle_ports)
                            b=bool(middle&set(ports))
                            actual=sum(all(p[i]==i for i in ports) for p in group)
                            need(actual==3**(t+1-h-b),'shared central phase and every port pair')
                        fixedfree=sum(all(p[i]!=i for i in range(len(adj))) for p in group)
                        need(fixedfree==(2**(t+1) if not u else 0),
                             'synchronized fixed-free mass')
    return receipts


def compositions(total, slots):
    if slots==1:
        yield (total,)
    else:
        for first in range(total+1):
            for rest in compositions(total-first,slots-1):
                yield (first,)+rest


def grammar():
    counts=Counter()
    # Exhaust all two-edge endpoint templates and gap words to total order eight.
    for m in (3,4):
        for defects in combinations(list(combinations(range(m),2)),2):
            if set(sum((list(e) for e in defects),[]))!=set(range(m)):continue
            for extra in range(9-m):
                for sizes in compositions(extra,m+1):
                    for gaps in product(*(list(inherited.atom_words(s)) for s in sizes)):
                        blocks=[]; endpoints=[]; qweights=[]
                        for i,gap in enumerate(gaps):
                            for w in gap:
                                blocks.append(inherited.atom_word((w,)))
                                qweights.append(w)
                            if i<m:
                                endpoints.append(len(blocks));blocks.append([0]);qweights.append(1)
                        q=old.adjacency(len(blocks),0)
                        for i,j in defects:
                            aa,bb=endpoints[i],endpoints[j]
                            q[aa]^=1<<bb;q[bb]^=1<<aa
                        # Expand the actual quotient while keeping every label.
                        n=sum(qweights); adj=[0]*n; starts=[]; offset=0
                        for block in blocks:
                            starts.append(offset)
                            for i,row in enumerate(block):adj[offset+i]|=row<<offset
                            offset+=len(block)
                        for i in range(len(blocks)):
                            for j in range(len(blocks)):
                                if q[i]>>j&1:
                                    for v in range(starts[i],starts[i]+qweights[i]):
                                        adj[v]|=((1<<qweights[j])-1)<<starts[j]
                        wf=weighted_feedback(q,qweights)
                        accepted=all(local_tests(q,qweights))
                        qgroup=[p for p in old.automorphisms(q)
                                if all(qweights[i]==qweights[p[i]] for i in range(len(q)))]
                        need(accepted==(wf==2 and len(qgroup)==1),
                             'local grammar equals weighted feedback and coloured rigidity')
                        group=old.automorphisms(adj)
                        f=old.feedback(adj)[0]
                        support=sum(any(p[i]!=i for p in group) for i in range(n))
                        need(not accepted or (3*f-support==6 and len(group)==3**qweights.count(3)),
                             'every accepted endpoint grammar has exact slack and group')
                        counts['presentations']+=1
                        counts['accepted' if accepted else 'rejected']+=1
    return dict(counts)


def smooth_and_prime():
    # Exact finite mass-restricted partition counts via two different products.
    partitions=[]
    for p in (3,5,7,11,17):
        for M in (12,24,48,96):
            allowed=[m for m in range(1,M+1) if old.smooth_part(m,p)==m]
            one=[1]+[0]*M
            two=[1]+[0]*M
            for m in allowed:
                for j in range(m,M+1):one[j]+=one[j-m]
                for _ in range(2):
                    for j in range(m,M+1):two[j]+=two[j-m]
            convolution=[sum(one[j]*one[k-j] for j in range(k+1)) for k in range(M+1)]
            need(two==convolution,'two-colour smooth product equals convolution')
            # Exact exponential majorant at x=1/2, without floating analysis.
            product_value=Q(1)
            for m in allowed:product_value/=(1-Q(1,2**m))**2
            need(sum(two)<=2**M*product_value,'finite partition Euler majorant exact')
            partitions.append({'p':p,'M':M,'parts':len(allowed),'count':sum(two)})
    banks=[]
    for L in (5,10,12,20,30,50):
        ps=[p for p in old.primes_to(2*L) if p>L]
        q=next(p for p in old.primes_to(4*L) if p>2*L)
        S=sum(ps); horizon=(q+1)*S
        profiles=set(); cells=defaultdict(list)
        for bits in product((0,1),repeat=len(ps)):
            subset=tuple(p for p,b in zip(ps,bits) if b)
            history=tuple(q*p for p in subset)
            vec=tuple(sum(gcd(a,j) for a in history) for j in (1,)+tuple(ps))
            expected=(len(subset),)+tuple(len(subset)+(p-1)*(p in subset) for p in ps)
            need(vec==expected,'prime-batch lower bank exact responses')
            need(sum(history)+S<=horizon,'lower bank all future words affordable')
            need(all(old.largest_prime(a)==q for a in history),'past lies in one completed prime batch')
            profiles.add(vec);cells[len(subset),sum(subset)].append(subset)
        need(len(profiles)==2**len(ps),'prime-batch lower profiles distinct')
        cell=max(cells.values(),key=len)
        mat=[[Q(2**sum(p-1 for p in set(a)&set(b))) for b in cell] for a in cell]
        for i in range(len(mat)):
            pivot=mat[i][i]
            need(pivot>0,'common-degree continuation matrix positive definite')
            for j in range(i+1,len(mat)):
                for k in range(i+1,len(mat)):
                    mat[j][k]-=mat[j][i]*mat[i][k]/pivot
        banks.append({'L':L,'primes':ps,'q':q,'horizon':horizon,
                      'profiles':len(profiles),'common_cell_rank':len(cell)})
    return {'smooth_partitions':partitions,'prime_banks':banks}


def geometry():
    rows=[]
    for colours in (3,4,5,6):
        triples=list(permutations(range(colours),3))
        neighbours=[{j for j,b in enumerate(triples) if all(x!=y for x,y in zip(a,b))}
                    for a in triples]
        support2=0;minimum=None;maximum=0
        for i,a in enumerate(triples):
            for j,b in enumerate(triples):
                count=len(neighbours[i]&neighbours[j])
                support2+=count>0
                minimum=count if minimum is None else min(minimum,count)
                maximum=max(maximum,count)
                if colours>=4:need(count>0,'four through six-colour two-link extension')
                if colours==6:need(count>=24,'six-colour greedy count lower bound')
                if colours==3:
                    pos={x:k for k,x in enumerate(a)}
                    perm=[pos[x] for x in b]
                    parity=sum(perm[u]>perm[v] for u in range(3) for v in range(u+1,3))%2
                    need((count>0)==(parity==0),'three-colour orientation coset exact')
        rows.append({'colours':colours,'states':len(triples),'one_link_support':sum(map(len,neighbours)),
                     'two_link_support':support2,'minimum_extensions':minimum,'maximum_extensions':maximum})
    # Exact Q(sqrt3) arithmetic: a+b sqrt3 represented by two Fractions.
    def add(x,y):return x[0]+y[0],x[1]+y[1]
    def mul(x,y):return x[0]*y[0]+3*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
    def sub(x,y):return x[0]-y[0],x[1]-y[1]
    z=[((Q(0),Q(0)),(Q(0),Q(0))),((Q(1),Q(0)),(Q(0),Q(0))),
       ((Q(1,2),Q(0)),(Q(0),Q(1,2)))]
    vertices=[(add(x,(Q(3*j,5),Q(0))),add(y,(Q(4*j,5),Q(0)))) for j in range(13) for x,y in z]
    need(len(set(vertices))==len(vertices),'triangle strip vertices distinct')
    edges=0
    for a,b in combinations(range(len(vertices)),2):
        dx=sub(vertices[a][0],vertices[b][0]);dy=sub(vertices[a][1],vertices[b][1])
        distance=add(mul(dx,dx),mul(dy,dy))
        expected=(a//3==b//3) or (b//3-a//3==1 and a%3==b%3)
        need((distance==(1,0))==expected,'induced unit-distance strip exact quadratic coordinates')
        edges+=distance==(1,0)
    return {'colour_transfers':rows,'geometry_vertices':len(vertices),'geometry_edges':edges}


def collatz():
    words=0;physical=0;multi=defaultdict(list)
    for length in range(1,7):
        for word in product(range(1,5),repeat=length):
            K=sum(word);B=0
            H=0
            for a in word:
                B=3*B+2**H;H+=a
            residue=((2**K-B)*pow(3**length,-1,2**(K+1)))%2**(K+1)
            need(residue%2==1,'Collatz word cylinder odd')
            for n in (residue,residue+2**(K+1),residue+3*2**(K+1)):
                x=n;clock=0
                for a in word:
                    x=3*x+1;valuation=0
                    while x%2==0:x//=2;valuation+=1
                    need(valuation==a,'cylinder realizes every exact valuation')
                    clock+=valuation
                need(clock==K and 2**K*x==3**length*n+B,'source endpoint and integer clock')
                need((x<n)==(2**K>3**length and n*(2**K-3**length)>B),
                     'descent iff exact source threshold')
                physical+=1
            multi[tuple(sorted(word))].append((B,word))
            words+=1
    for key,rows in multi.items():
        least=min(b for b,w in rows);most=max(b for b,w in rows)
        need(next(b for b,w in rows if w==key)==least,'increasing word minimizes carry')
        need(next(b for b,w in rows if w==key[::-1])==most,'decreasing word maximizes carry')
    # Necessity, separately: actual source walks must recover their cylinders.
    for n in range(1,2002,2):
        x=n;B=0;K=0
        for m in range(1,13):
            B=3*B+2**K
            x=3*x+1;a=0
            while x%2==0:x//=2;a+=1
            K+=a
            residue=((2**K-B)*pow(3**m,-1,2**(K+1)))%2**(K+1)
            need(n%2**(K+1)==residue,'actual sources recover exact cylinders')
    return {'formal_words':words,'constructed_physical_prefixes':physical,'multisets':len(multi),
            'actual_source_prefixes':1001*12}


def prime_costs():
    rows=[]
    for p in (5,7,11,13,17):
        cost=3*(p*p-1)//8-p
        for ell in range(1,302,2):
            charge=3*(ell*ell-1)//8-ell if ell>=3 else 0
            power=0;d=p
            while d<=ell:power+=ell//d;d*=p
            need(charge>=cost*power,'factorial valuation paid by orbit slack')
        adj=[sum(1<<j for j in range(p) if 1<=(j-i)%p<=(p-1)//2) for i in range(p)]
        # Independent internal-order witness; exponential DP through p=13.
        backwards=sum(bool(adj[j]>>i&1) for i in range(p) for j in range(i+1,p))
        need(backwards==(p*p-1)//8,'circular tournament attains feedback lower bound')
        if p<=13:need(old.feedback(adj)[0]==backwards,'circular optimum subset-DP control')
        rows.append({'p':p,'slack_cost':cost,'feedback':backwards})
    # Direct divisor definitions, including the prime properness boundary.
    for n in range(2,1001):
        divs=[d for d in range(2,n) if n%d==0]
        sf=[d for d in divs if all(d%(p*p) for p in old.primes_to(int(d**.5))) ]
        primes=[d for d in divs if old.largest_prime(d)==d]
        z=n;exponents=[]
        for p in old.primes_to(n):
            if z%p:continue
            e=0
            while z%p==0:z//=p;e+=1
            exponents.append(e)
        need((len(divs)==len(sf)+len(primes))==(sorted(exponents) in ([1],[3],[1,1,2])),
             'inherited proper divisor balance direct definitions')
    need((3**2*5*7,2*1+3+6,2*3+5+7)==(315,11,18),
         'mixed balanced group-order witness exact invariants')
    return rows


def main():
    result={'status':'FINITE-EXACT; same-agent checks, no independent-agent audit',
            'census':census(),'synchronized':synchronized_ports(),'grammar':grammar(),
            'prime_mass':smooth_and_prime(),'geometry':geometry(),
            'collatz':collatz(),'prime_costs':prime_costs()}
    result['checks']=dict(CHECKS);result['total_checks']=sum(CHECKS.values())
    result['dependencies']={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in (PREV/'verify.py',PREV/'independent_controls.py')}
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=='__main__':main()
