#!/usr/bin/env python3
"""Independent, bounded Session 32 audit. Standard library; assertions not used.
Reads producer JSON but imports/runs no producer implementation. Writes nothing:
  python3 audit_graphs.py > certificate.json
  python3 -O audit_graphs.py > certificate.optimized.json
"""
from pathlib import Path
from itertools import combinations, product
from math import isqrt, comb, prod
from collections import Counter
import hashlib, json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
checks=0

def need(value, message):
    global checks
    checks+=1
    if not value: raise RuntimeError(message)


def edges(n):
    # Construct by square-reflection pairs, not pairwise square testing.
    return sorted((u,q*q-u) for q in range(2,isqrt(2*n-1)+1)
                  for u in range(max(1,q*q-n),(q*q-1)//2+1)
                  if q*q-u<=n)


def partition(V,E):
    p={v:v for v in V}
    def root(v):
        while p[v]!=v:
            p[v]=p[p[v]];v=p[v]
        return v
    cycles=0
    for a,b in E:
        a,b=root(a),root(b)
        if a==b:cycles+=1
        else:p[a]=b
    blocks={}
    for v in V:blocks.setdefault(root(v),[]).append(v)
    return sorted(tuple(sorted(x)) for x in blocks.values()),cycles


def degree(V,E):
    c=Counter(x for e in E for x in e)
    return {v:c[v] for v in V}


def normal_path(p):return min(tuple(p),tuple(reversed(p)))


def solve_degrees(V,E,endpoints):
    """Independent exact Hamilton edge-set solver, not vertex-order DFS.
    Propagate degree obligations, reject closed cycles, branch an edge.
    """
    V=tuple(V);E=tuple(E)
    incidence={v:tuple(i for i,e in enumerate(E) if v in e) for v in V}
    targets={v:1 if v in endpoints else 2 for v in V}
    answers=[];nodes=0
    def recur(state):
        nonlocal nodes
        nodes+=1
        while True:
            changes={}
            for v in V:
                unknown=[i for i in incidence[v] if state[i]<0]
                taken=sum(state[i]==1 for i in incidence[v])
                remaining=targets[v]-taken
                if remaining<0 or remaining>len(unknown):return
                if remaining in (0,len(unknown)):
                    value=int(remaining>0)
                    for i in unknown:
                        if i in changes and changes[i]!=value:return
                        changes[i]=value
            if not changes:break
            state=list(state)
            for i,value in changes.items():state[i]=value
        chosen=[E[i] for i,v in enumerate(state) if v==1]
        parts,cycles=partition(V,chosen)
        if cycles:return
        if all(v>=0 for v in state):
            if len(parts)==1:answers.append(tuple(chosen))
            return
        # Fail if remaining edges cannot reconnect the currently chosen pieces.
        if len(partition(V,[E[i] for i,v in enumerate(state) if v!=0])[0])!=1:return
        candidates=[(sum(state[j]<0 for j in incidence[v]),v) for v in V
                    if any(state[j]<0 for j in incidence[v])]
        _,v=min(candidates)
        i=next(j for j in incidence[v] if state[j]<0)
        for val in (0,1):
            new=list(state);new[i]=val;recur(new)
    recur([-1]*len(E))
    return answers,nodes


def path_from_edges(V,E):
    d=degree(V,E);ends=sorted(v for v in V if d[v]==1)
    need(len(ends)==2,'two reconstructed endpoints')
    a={v:[] for v in V}
    for u,v in E:a[u].append(v);a[v].append(u)
    result=[ends[0]]
    while len(result)<len(V):
        nxt=[v for v in a[result[-1]] if len(result)<2 or v!=result[-2]]
        need(len(nxt)==1,'unique path extension')
        result.append(nxt[0])
    need(len(set(result))==len(V),'path is spanning')
    return normal_path(result)

square=json.loads((ROOT/'square-sum/certificates.json').read_text())
ports=json.loads((ROOT/'ports/certificate.json').read_text())
source_files=['square-sum/README.md','square-sum/verify_square_sum.py',
              'square-sum/certificates.json','ports/verify_ports.py',
              'ports/certificate.json','ports/PORT_TRANSFER.md',
              'CAT_MAP_PATH_FORESTS.md',
              'check_catmap_forests.py']
if (ROOT/'ports/README.md').exists():source_files.append('ports/README.md')
pins={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in source_files}

# All-order degree-two proof reduces to seven bands. Beyond r=8 the polynomial
# r^2-6r-10 is positive and increasing, so 2*r^2-1 >= (r+3)^2.
bands=[]
for r in range(1,8):
    lo=max(r*r,((r+2)**2+2)//2)
    hi=min((r+1)**2-1,(r+3)**2//2)
    if lo<=hi:bands.append([lo,hi])
need(bands==[[13,15],[19,24],[25,32],[36,40],[49,50]],'all degree-two bands')
need(8*8-6*8-10>0 and 2*8+1-6>0,'infinite lower bound and monotonicity')
need([n for a,b in bands for n in range(a,b+1)]==square['degree_two_births_all'],'birth list')
arrivals=[]
for n in range(1,65):
    E=edges(n);V=range(1,n+1);P,_=partition(V,E)
    oldE=edges(n-1) if n>1 else [];oldP,_=partition(range(1,n),oldE)
    new=sorted(u for u,v in E if v==n)
    touch=sum(any(v in C for v in new) for C in oldP)
    d=isqrt(2*n-1)-isqrt(n)
    need(len(new)==d,'arrival degree '+str(n))
    need(len(P)==len(oldP)+1-touch,'component response')
    beta=len(E)-n+len(P);oldbeta=len(oldE)-(n-1)+len(oldP)
    need(beta==oldbeta+d-touch,'cycle response')
    need((len(P)==1)==(n==1 or n>=14),'connectivity control')
    deg=degree(V,E)
    for x in V:
        need(deg[x]==isqrt(x+n)-isqrt(x)-int(isqrt(2*x)**2==2*x),'full degree')
    rec=dict(n=n,neighbors=new,degree=d,touched_components=touch,
             components=len(P),cycle_rank=beta)
    arrivals.append(rec)
    if n<=36:need(rec==square['arrival'][n-1],'producer arrival ledger')

# Verify all listed paths, and independently establish exact counts through25
# using degree propagation, plus elementary cut certificates.
census=[]
for entry in square['census']:
    n=entry['n'];E=set(edges(n));V=set(range(1,n+1))
    for p in entry['paths']:
        need(len(p)==n and set(p)==V,'certificate spans vertices')
        need(all(tuple(sorted(e)) in E for e in zip(p,p[1:])),'certificate uses edges')
    need(len({normal_path(p) for p in entry['paths']})==entry['count'],'distinct certificate paths')
for c in square['cut_certificates']:
    n=c['n'];cut=set(c['cut'])
    parts,_=partition(set(range(1,n+1))-cut,[e for e in edges(n) if not set(e)&cut])
    need(parts==sorted(tuple(C) for C in c['remaining_components']),'cut components')
    need(len(parts)>len(cut)+1,'cut obstruction')
solver_nodes=0
for n in [14,15,16,17,23,24,25]:
    V=list(range(1,n+1));E=edges(n);deg=degree(V,E)
    leaves={v for v in V if deg[v]==1}
    found=set();pair_counts={}
    for pair in combinations(V,2):
        if not leaves<=set(pair):continue
        sols,nodes=solve_degrees(V,E,pair);solver_nodes+=nodes
        if sols:pair_counts[str(pair)]=len(sols)
        for selected in sols:found.add(path_from_edges(V,selected))
    expected={normal_path(p) for p in square['census'][n-1]['paths']}
    need(found==expected,'independent degree-solver exact set at '+str(n))
    census.append(dict(n=n,count=len(found),endpoint_counts=pair_counts))

# A G23 degree-two endpoint demonstrates why endpoint-blind suppression is wrong.
hostile=next(p for p in square['census'][22]['paths'] if 22 in (p[0],p[-1]))
need(degree(range(1,24),edges(23))[22]==2,'degree-two endpoint witness')
for n,a,b in [(15,1,10),(23,2,13),(25,11,24)]:
    need(square['census'][n-2]['count']==0,'old graph nontraceable')
    for p in square['census'][n-1]['paths']:
        i=p.index(n)
        need(0<i<n-1 and {p[i-1],p[i+1]}=={a,b},'forced ear')

# For n24 record each proposed endpoint and its exact failure under degree
# constraints. Some failures arise only once proper cycles are forbidden.
cycle24=tuple(square['proper_cycle24'])
need(set(cycle24).isdisjoint({18,2,9,11,22}),'proper-cycle vertices all internal')
need(all(tuple(sorted(e)) in set(edges(24)) for e in
         zip(cycle24,cycle24[1:]+cycle24[:1])),'proper-cycle square edges')

# Local signature reconstruction by a single edge-subset bit mask, independent
# of the producer's nested internal-forest/port enumeration and implementation.
allE=set(edges(25));block_states={};signature_report={}
def signature(V,internal,selected):
    internal=tuple(sorted(internal));selected=tuple(sorted(selected))
    total=degree(V,internal+selected)
    if any(total[v] not in (1,2) for v in V):return None
    endpoints=tuple(sorted(v for v in V if total[v]==1))
    if len(endpoints)>2:return None
    parts,cycles=partition(V,internal)
    if cycles:return None
    pairs=[]
    for C in parts:
        boundary=tuple(sorted(e for e in selected if e[0] in C))
        ends=tuple(v for v in endpoints if v in C)
        if not boundary:return None
        need(len(boundary)+len(ends)==2,'independent terminal identity')
        pairs.append((boundary,ends))
    return selected,endpoints,tuple(sorted(pairs))

for name,b in ports['blocks'].items():
    V=tuple(b['vertices']);S=set(V)
    internal=sorted(e for e in allE if set(e)<=S)
    boundary=sorted(e if e[0] in S else e[::-1] for e in allE if len(set(e)&S)==1)
    candidates=internal+boundary;new={};fullstates=[]
    for mask in range(1<<len(candidates)):
        selected=[e for i,e in enumerate(candidates) if mask>>i&1]
        ie=tuple(e for e in selected if set(e)<=S)
        pe=tuple(e for e in selected if len(set(e)&S)==1)
        sig=signature(V,ie,pe)
        if sig is None:continue
        need(sig not in new,'tree signature has unit multiplicity')
        new[sig]=tuple(sorted(tuple(sorted(e)) for e in selected))
        fullstates.append((sig,new[sig]))
    old={}
    for st in b['states']:
        pe=[tuple(b['ports'][i]) for i in st['selected_ports']]
        ie=[tuple(e) for e in st['internal_edges']]
        sig=signature(V,ie,pe)
        need(sig is not None,'producer local state admissible')
        advertised=tuple(sorted((tuple(sorted(tuple(b['ports'][i]) for i in p['ports'])),
                                 tuple(sorted(p['endpoints']))) for p in st['pairs']))
        need(advertised==sig[2] and tuple(st['endpoints'])==sig[1],'producer terminal pairing')
        old[sig]=tuple(sorted(tuple(sorted(e)) for e in ie+pe))
    need(new==old,'complete local signatures '+name)
    need(len(new)==b['profile_count'],'profile count '+name)
    block_states[name]=fullstates
    anonymous=Counter(tuple((pair[0],len(pair[1])) for pair in sig[2]) for sig in new)
    hist={str(k):v for k,v in sorted(Counter(anonymous.values()).items())}
    need(len(anonymous)==b['anonymous_endpoint_profile_count'],'anonymous signature count')
    need(hist==b['anonymous_multiplicity_histogram'],'anonymous realization weights')
    signature_report[name]={'count':len(new),'internal_edges':internal,'boundary_edges':boundary,
                            'anonymous_count':len(anonymous),'anonymous_multiplicity_histogram':hist}

# Glue using the FULL expanded edge sets and DSU, not the compressed virtual graph.
core=ports['core'];coreE=[tuple(e) for e in ports['core_edges']]
need(set(coreE)=={e for e in allE if set(e)<=set(core)},'core edge inventory')
# Independently contract the named tree blocks, preserving multiedges.
owner={v:26+i for i,b in enumerate(ports['blocks'].values()) for v in b['vertices']}
quotient=sorted(tuple(sorted((owner.get(u,u),owner.get(v,v)))) for u,v in allE
                if owner.get(u,u)!=owner.get(v,v))
Qvertices=core+[26,27,28]
qparts,qbeta=partition(Qvertices,quotient)
need(len(Qvertices)==13 and len(quotient)==20 and qbeta==8,'multigraph contraction')
need(len(qparts)==1,'quotient remains connected')
simple=set(quotient)
need(len(simple)==18 and partition(Qvertices,simple)[1]==6,'simplification loses two cycles')
need({e:m for e,m in Counter(quotient).items() if m>1}=={(15,28):2,(19,28):2},'precise parallel pairs')
coreparts,_=partition(core,coreE)
need(len(coreparts)==6>3+1,'quotient Hamiltonian cut obstruction')
need(quotient==sorted(tuple(e) for e in ports['bare_tree_quotient']['edges']),'producer quotient')
wholeV=list(range(1,26));budget=feasible=0;solutions=set();bad=[]
for triple in product(*block_states.values()):
    ends=[v for sig,_ in triple for v in sig[1]]
    if len(ends)>2:continue
    budget+=1
    partial=[e for _,selected in triple for e in selected]
    for mask in range(1<<len(coreE)):
        E=partial+[e for i,e in enumerate(coreE) if mask>>i&1]
        d=degree(wholeV,E)
        if any(d[v] not in (1,2) for v in wholeV):continue
        if sum(d[v]==1 for v in wholeV)!=2:continue
        feasible+=1
        parts,cyc=partition(wholeV,E)
        if len(parts)>1:
            bad.append({'components':parts,'cycle_rank':cyc});continue
        need(cyc==0,'connected degree-feasible graph is path')
        solutions.add(path_from_edges(wholeV,E))
need(budget==ports['endpoint_budget_profile_triples']==1543,'endpoint budget count')
need(feasible==ports['degree_feasible_gluings']==19,'degree feasible count')
need(len(bad)==ports['disconnected_degree_feasible_gluings']==9,'disconnected count')
need(solutions=={normal_path(s['path']) for s in ports['solutions']},'global transfer set')
need(len(solutions)==10,'ten paths')

# Independent matrix-forest count: each acyclic edge subset contributes the
# product of component sizes, with coefficient indexed by selected edge count.
def rooted_forest_polynomial(V,E):
    coeff=[0]*len(V)
    for mask in range(1<<len(E)):
        selected=[e for i,e in enumerate(E) if mask>>i&1]
        parts,cyc=partition(V,selected)
        if not cyc:coeff[len(selected)]+=prod(len(C) for C in parts)
    return coeff

def path_formula(m):return [comb(2*m-k-1,k) for k in range(m)]
def det_bareiss(matrix):
    matrix=[row[:] for row in matrix];previous=1
    for k in range(len(matrix)-1):
        pivot=matrix[k][k]
        need(pivot!=0,'positive principal minor in Bareiss')
        for i in range(k+1,len(matrix)):
            for j in range(k+1,len(matrix)):
                value=matrix[i][j]*pivot-matrix[i][k]*matrix[k][j]
                need(value%previous==0,'exact Bareiss division')
                matrix[i][j]=value//previous
        for i in range(k+1,len(matrix)):matrix[i][k]=0
        previous=pivot
    return matrix[-1][-1]
def forest_determinant(m,t):
    M=[[int(i==j) for j in range(m)] for i in range(m)]
    for i in range(m-1):
        M[i][i]+=t;M[i+1][i+1]+=t;M[i][i+1]-=t;M[i+1][i]-=t
    return det_bareiss(M)
def mul(A,B):return [[sum(A[i][k]*B[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
def matpow(A,n):
    R=[[1,0],[0,1]]
    for _ in range(n):R=mul(R,A)
    return R
A=[[2,1],[1,1]];B=[[3,-1],[1,0]];U=[[1,-1],[0,1]]
need(mul(A,U)==mul(U,B),'integral conjugacy')
need(U[0][0]*U[1][1]-U[0][1]*U[1][0]==1,'unimodular conjugator')
fib=[0,1]
for k in range(2,65):fib.append(fib[-1]+fib[-2])
forest_rows=[]
for m in range(1,17):
    f=path_formula(m)
    if m<=12:
        need(f==rooted_forest_polynomial(list(range(m)),[(i,i+1) for i in range(m-1)]),'path composition formula')
        # Both polynomials have degree at most m-1, so m distinct t values
        # certify equality of the full determinant and forest polynomial.
        for t in range(m):
            need(forest_determinant(m,t)==sum(c*t**k for k,c in enumerate(f)),
                 'exact determinant interpolation control')
    if m>=3:
        prev=path_formula(m-1);pp=path_formula(m-2)
        recur=[0]*m
        for k,c in enumerate(prev):recur[k]+=c;recur[k+1]+=2*c
        for k,c in enumerate(pp):recur[k+2]-=c
        need(f==recur,'coefficientwise forest recurrence')
    need(sum(f)==fib[2*m],'even Fibonacci forest count')
    need(matpow(B,m-1)[0][0]==sum(f),'boundary matrix observable')
    forest_rows.append({'m':m,'coefficients':f,'rooted_forest_count':sum(f),
                        'cat_map_trace':sum(matpow(A,m)[i][i] for i in range(2))})
need(sum(matpow(A,2)[i][i] for i in range(2))==7 and sum(path_formula(4))==21,'distinct observables')
blocks_forest={name:rooted_forest_polynomial(b['vertices'],[tuple(e) for e in b['internal_edges']])
               for name,b in ports['blocks'].items()}
need(blocks_forest=={'A':[1,8,20,18,5],'B':[1,6,10,4],'C':[1,10,36,56,35,6]},'block forest polynomials')
need(prod(sum(path_formula(m)) for m in (3,4,5))==9240,'G12 forest count')

out={'status':'PASS independent bounded structural audit', 'explicit_checks':checks,
     'source_sha256':pins, 'scope':{'arrival_control':'n=1..64; all-order degree-two bound proved separately',
      'independent_exact_path_sets':[14,15,16,17,23,24,25], 'all_listed_paths_validated':'n=1..25',
      'forest_subset_check':'path orders1..12 and the three G25 blocks',
      'not_audited':'Gerbicz large finite base/glue tables; FGST representation claims'},
     'degree_two_bands':bands,'arrival_12_to_15':arrivals[11:15],
     'degree_solver_nodes':solver_nodes,'independent_census':census,
     'degree_two_endpoint22_witness':hostile,'local_signatures':signature_report,
     'gluing':{'endpoint_budget_triples':budget,'degree_feasible':feasible,
               'disconnected':len(bad),'hamiltonian_paths_mod_reversal':len(solutions),
               'canonical_paths':sorted(solutions),'disconnected_expanded_examples':bad[:3]},
     'contraction':{'vertices':13,'multigraph_edges':20,'multigraph_cycle_rank':8,
       'simple_edges':18,'simple_cycle_rank':6,'core_components':coreparts},
     'cat_map':{'A':A,'B':B,'U':U,'AU':mul(A,U),'UB':mul(U,B)},
     'forest_rows':forest_rows,'block_forest_polynomials':blocks_forest,
     'G12_rooted_forest_count':9240}
print(json.dumps(out,indent=2))
