"""Actual ordinary tournament DAG: reverse binary potential, then K-bit stored-residue transfer.

No residue is consulted during structural graph construction.  This module does
not change clone/excess semantics: it deliberately implements ordinary e=0 only.
All validation gates survive python -O.  No network or third-party dependencies.
"""
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
from math import factorial, gcd, prod
from pathlib import Path
import argparse
import json


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def valid_n(n):
    if type(n) is not int or n < 0:
        raise ValueError('n must be a nonnegative integer')


def valid_k(K):
    if type(K) is not int or K < 1:
        raise ValueError('K must be a positive integer')


def lower_bound(D):
    return max(((x + 2*i - 1)//(2*i) for i,x in enumerate(D) if i), default=0)


def vfact(n):
    return n - n.bit_count()


def odd_factorials(n, modulus=None):
    values = [1]
    for j in range(1,n+1):
        value = values[-1] * (j // (j & -j))
        values.append(value if modulus is None else value % modulus)
    return values


@dataclass(frozen=True)
class Edge:
    target: int
    ell: int
    k: int
    u: int
    shift: int


@dataclass
class OrdinaryDAG:
    n: int
    vertices: list
    edges: list
    layers: list
    potential: list
    prefix_cost: list

    @property
    def baseline(self):
        return self.potential[0]

    @property
    def terminal(self):
        return len(self.vertices)-1

    def reduced_cost(self, source, edge):
        return edge.shift + self.potential[edge.target] - self.potential[source]

    def statistics(self, K=None):
        out = {'vertices_including_terminal':len(self.vertices),
               'input_states':len(self.vertices)-1,
               'edges':sum(map(len,self.edges)), 'baseline':self.baseline,
               'maximum_suffix_potential':max(self.potential)}
        if K is not None:
            valid_k(K)
            active = {v for v,c in enumerate(self.prefix_cost) if c < K}
            out.update({'K':K,'active_vertices_including_terminal':len(active),
                        'active_input_states':len(active-{self.terminal}),
                        'active_edges':sum(self.prefix_cost[v]+self.reduced_cost(v,e)<K
                                           for v in active for e in self.edges[v])})
        return out


def build_dag(n):
    """Generate full reachable stage/(S,D) DAG before looking at any weights."""
    valid_n(n)
    if n == 0:
        return OrdinaryDAG(0,[(-1,0,())],[[]],[[0]],[0],[0])
    top = n if n%2 else n-1
    vertices = [(top,n,(0,)*((top+1)//2))]
    edges = [[]]
    layers = [[0]]
    current = {(n,vertices[0][2]):0}
    for ell in range(top,0,-2):
        upcoming = {}
        for (s,D),source in sorted(current.items()):
            need(vertices[source] == (ell,s,D), 'stage is part of vertex identity')
            selected = D[(ell-1)//2] if ell <= s else 0
            b = lower_bound(D)
            for k in range(s//ell+1):
                rest = s-ell*k
                if ell == 1 and rest:
                    continue
                Dp = tuple(D[(d-1)//2]+k*(d-gcd(ell,d))
                           for d in range(1,min(ell-2,rest)+1,2))
                bp = lower_bound(Dp)
                shift = (k*(ell-1)//2 + ell*k*(k-1)//2 - vfact(k)
                         + k*(ell*b-selected) + rest*(k+b-bp))
                need(shift >= 0 and bp <= b+k, 'integral gauge edge')
                key = (rest,Dp)
                if key not in upcoming:
                    upcoming[key] = len(vertices)
                    vertices.append((ell-2,rest,Dp)); edges.append([])
                destination = upcoming[key]
                edges[source].append(Edge(destination,ell,k,n-s,shift))
        current = upcoming
        layers.append(sorted(current.values()))
    need(len(current)==1 and (0,()) in current, 'one structural terminal')
    need(current[(0,())] == len(vertices)-1, 'terminal last in topological order')
    V = [None]*len(vertices); V[-1]=0
    for source in range(len(vertices)-2,-1,-1):
        need(bool(edges[source]),'every reachable state has an all-singleton suffix')
        V[source] = min(e.shift+V[e.target] for e in edges[source])
    P = [None]*len(vertices); P[0]=0
    for source,outgoing in enumerate(edges):
        need(P[source] is not None,'structural prefix reachability')
        for e in outgoing:
            reduced = e.shift+V[e.target]-V[source]
            need(reduced >= 0,'nonnegative reduced cost')
            candidate = P[source]+reduced
            P[e.target] = candidate if P[e.target] is None else min(P[e.target],candidate)
    baseline = n//2 if n>=3 else 0
    need(V[0]==baseline,('proved termwise baseline',n,V[0],baseline))
    return OrdinaryDAG(n,vertices,edges,layers,V,P)


def edge_factor(edge, odd, modulus=None):
    numerator = odd[edge.u+edge.ell*edge.k]
    if modulus is None:
        denominator = odd[edge.u]*edge.ell**edge.k*odd[edge.k]
        factor,remainder = divmod(numerator,denominator)
        need(remainder==0 and factor>0 and factor%2==1,'exact positive odd transfer')
        return factor
    denominator = odd[edge.u]*pow(edge.ell,edge.k,modulus)%modulus
    denominator = denominator*odd[edge.k]%modulus
    need(denominator%2==1,'modular denominator is an odd unit')
    return numerator*pow(denominator,-1,modulus)%modulus


def normalized_modular(dag,K,prune=True):
    """Store only K-bit residues; no exact counts or B+K-bit absolute quotient.

    Python multiplication/shift temporaries can exceed K bits before reduction.
    """
    valid_k(K)
    modulus = 1<<K; odd=odd_factorials(dag.n,modulus)
    values = [0]*len(dag.vertices); values[0]=1
    operations = 0; skipped = 0
    for source,outgoing in enumerate(dag.edges):
        if prune and dag.prefix_cost[source]>=K:
            skipped += len(outgoing)
            continue
        for e in outgoing:
            shift=dag.reduced_cost(source,e)
            if prune and dag.prefix_cost[source]+shift>=K:
                skipped += 1
                continue
            contribution = 0 if shift>=K else ((edge_factor(e,odd,modulus)*values[source])<<shift)%modulus
            values[e.target]=(values[e.target]+contribution)%modulus
            operations += 1
    return values[-1]*pow(odd[dag.n],-1,modulus)%modulus, {
        'edge_updates':operations, 'valuation_filtered_edges':skipped,
        'residue_precision_bits':K}, values


def exact_audit(dag):
    """Independent exact Y/Z passes certify every reachable coordinate."""
    odd=odd_factorials(dag.n)
    Y=[0]*len(dag.vertices); Z=[0]*len(dag.vertices); Y[0]=Z[0]=1
    bits_y=bits_z=negative_coordinate_scales=0
    for source,outgoing in enumerate(dag.edges):
        exponent=dag.potential[source]-dag.baseline
        expected=Fraction(Y[source])*(Fraction(2)**exponent)
        need(expected.denominator==1 and expected==Z[source],('exact prefix gauge',dag.n,source))
        negative_coordinate_scales += exponent<0
        bits_y+=Y[source].bit_length(); bits_z+=Z[source].bit_length()
        for e in outgoing:
            factor=edge_factor(e,odd)
            Y[e.target]+=(factor*Y[source])<<e.shift
            Z[e.target]+=(factor*Z[source])<<dag.reduced_cost(source,e)
    need(Y[-1]==Z[-1]<<dag.baseline,'exact terminal scaling')
    ordinary,remainder=divmod(Y[-1],odd[-1])
    normalized,rem2=divmod(Z[-1],odd[-1])
    need(remainder==rem2==0,'exact final odd division')
    need(ordinary==normalized<<dag.baseline,'count normalization')
    return ordinary,normalized,Y,Z,{'prefix_coordinates':len(Y),
        'negative_coordinate_scale_exponents':negative_coordinate_scales,
        'summed_exact_Y_value_bits':bits_y,'summed_exact_Z_value_bits':bits_z,
        'max_exact_Y_bits':max(x.bit_length() for x in Y),
        'max_exact_Z_bits':max(x.bit_length() for x in Z)}


def odd_partitions(n,low=1):
    if not n:
        yield ()
    for d in range(low,n+1,2):
        for tail in odd_partitions(n-d,d):
            yield (d,)+tail


def partition_oracle(n):
    total=Fraction(0); valuations=[]
    for part in odd_partitions(n):
        lengths=set(part); counts={d:part.count(d) for d in lengths}
        exponent=sum((d-1)//2 for d in part)+sum(gcd(a,b) for a,b in combinations(part,2))
        denominator=prod(d**k*factorial(k) for d,k in counts.items())
        total+=Fraction(1<<exponent,denominator)
        valuations.append(exponent-sum(vfact(k) for k in counts.values()))
    need(total.denominator==1,'unmerged Burnside integrality')
    return int(total),min(valuations)


def class_enumeration(n):
    """Every labelled tournament, no structural filter, explicit S_n orbits."""
    pairs=list(combinations(range(n),2)); index={p:i for i,p in enumerate(pairs)}
    transformations=[]
    for permutation in permutations(range(n)):
        transform=[]
        for i,j in pairs:
            a,b=permutation[i],permutation[j]
            transform.append((index[tuple(sorted((a,b)))],a>b))
        transformations.append(transform)
    seen=set(); count=0
    for mask in range(1<<len(pairs)):
        if mask in seen:
            continue
        count+=1
        for transform in transformations:
            image=0
            for source,(target,reverse) in enumerate(transform):
                image|=(((mask>>source)&1)^reverse)<<target
            seen.add(image)
    need(len(seen)==1<<len(pairs),'full labelled universe exhausted')
    return count,1<<len(pairs)


def load_reference(repo):
    path=repo/'05-knowledge/results/port_residual_integral_20261025/counting/check.py'
    spec=spec_from_file_location('inherited_integral_reference',path)
    module=module_from_spec(spec); spec.loader.exec_module(module)
    return module


def default_repo():
    here=Path(__file__).resolve()
    for candidate in [here.parent.parent/'session64-repo',*here.parents]:
        if (candidate/'05-knowledge/results/port_residual_integral_20261025/counting/check.py').is_file():
            return candidate
    raise FileNotFoundError('supply --repo pointing to the research repository')


def controls(repo):
    inherited=load_reference(repo)
    known=[1,1,1,2,4,12,56,456,6880,191536,9733056,903753248,154108311168]
    prefixes=modular_comparisons=partition_comparisons=0; examples=[]
    samples=[]; seen_same_key_different_stages=False; retained_zero_residue=False
    for n in list(range(51))+[60,80,100]:
        dag=build_dag(n)
        count,normalized,Y,Z,detail=exact_audit(dag)
        prior,stats=inherited.integral_target(n,False)
        need(count==prior,('inherited exact engine',n))
        if n:
            need(stats['keys']==len(dag.vertices)-1,('unchanged ordinary input graph',n))
        prefixes+=detail['prefix_coordinates']
        if n<len(known): need(count==known[n],('known unrestricted prefix',n))
        if n<=27:
            direct,baseline=partition_oracle(n)
            need((count,dag.baseline)==(direct,baseline),('independent unmerged partitions',n))
            partition_comparisons+=1
        if n<=6:
            classes,labelled=class_enumeration(n)
            need(classes==count,('unrestricted full class census',n))
            examples.append({'n':n,'classes':classes,'labelled_tournaments':labelled})
        key_stages=defaultdict(set)
        for ell,s,D in dag.vertices: key_stages[(s,D)].add(ell)
        seen_same_key_different_stages |= any(len(v)>1 for v in key_stages.values())
        for K in range(1,9):
            residue,work,values=normalized_modular(dag,K)
            full,fullwork,fullvalues=normalized_modular(dag,K,False)
            need(residue==full==normalized%(1<<K),('normalized K-bit value',n,K))
            need(values==fullvalues==[x%(1<<K) for x in Z],('every modular intermediate value',n,K))
            need(work['edge_updates']==dag.statistics(K)['active_edges'],'exact filtered-operation count')
            retained_zero_residue |= any(y and y%(1<<K)==0 and z%(1<<K)!=0 for y,z in zip(Y,Z))
            modular_comparisons+=1
            if n in (30,50,60,80,100) and K in (2,8):
                samples.append(dict(n=n,**dag.statistics(K),**work))
        if n in (30,50,100):
            detail.update({'n':n}); examples.append(detail)
    need(seen_same_key_different_stages,'stage identity hostile exists in actual graph')
    need(retained_zero_residue,'absolute-zero normalized-nonzero hostile exercised')
    hostile_dag=build_dag(50)
    value,_,_=normalized_modular(hostile_dag,2)
    need(value==2 and hostile_dag.baseline==25,'order50 known precision hostile')
    naive=any(e.shift<hostile_dag.baseline for e in hostile_dag.edges[0])
    need(naive,'initial-only baseline subtraction would make a negative shift')
    for args in [(-1,), (True,), (1.5,)]:
        try: build_dag(*args)
        except ValueError: pass
        else: raise RuntimeError('invalid target accepted')
    for K in [0,-1,True,1.5]:
        try: normalized_modular(hostile_dag,K)
        except ValueError: pass
        else: raise RuntimeError('invalid precision accepted')
    return {'status':'PASS','ordinary_orders':{'consecutive':[0,50],'additional':[60,80,100]},
        'normalized_precision_bits':[1,8], 'prefix_coordinate_checks':prefixes,
        'modular_comparisons':modular_comparisons,'unmerged_partition_comparisons':partition_comparisons,
        'full_class_enumeration_and_exact_bit_samples':examples,'operation_samples':samples,
        'hostiles':{'same_key_distinct_stages':True,'absolute_zero_normalized_nonzero':True,
                    'initial_only_baseline_subtraction_invalid':naive,
                    'order50_baseline':25,'order50_normalized_mod4':value},
        'scope':'Unrestricted ordinary e=0; same structural graph, reverse pass plus optional exact valuation filter; no clone/core or runtime improvement claim.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--repo',type=Path)
    args=parser.parse_args()
    print(json.dumps(controls(args.repo or default_repo()),indent=2,sort_keys=True))
