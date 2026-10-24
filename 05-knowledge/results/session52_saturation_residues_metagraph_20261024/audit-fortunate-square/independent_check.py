#!/usr/bin/env python3
"""Independent Session52 audit. No producer imports; standard library only.
State space consists of acyclic degree-at-most-two edge sets, not vertex paths.
"""
from collections import Counter
from copy import deepcopy
from functools import lru_cache
from itertools import combinations, groupby
import hashlib
import json
from math import gcd, isqrt, prod
from pathlib import Path

ROOT=Path(__file__).resolve().parent
SNAP=ROOT/'frozen_candidate'
def need(ok,msg):
    if not ok: raise ValueError(msg)
def canon(o):return json.dumps(o,sort_keys=True,indent=2)+'\n'
def parts(V,E):
    labels={v:v for v in V}
    for a,b in E:
        x,y=labels[a],labels[b]
        labels={v:(x if t==y else t) for v,t in labels.items()}
    return sorted([sorted(v for v in V if labels[v]==k) for k in set(labels.values())])
def forest_ok(V,E):
    degree=Counter(x for e in E for x in e)
    if any(degree[v]>2 for v in V):return False
    return len(E)==len(V)-len(parts(V,E))
def to_paths(V,E):
    need(forest_ok(V,E),'not a linear forest')
    A={v:set() for v in V}
    for a,b in E:A[a].add(b);A[b].add(a)
    out=[]
    for C in parts(V,E):
        ends=[v for v in C if len(A[v])<2];need(bool(ends),'cycle')
        path=[min(ends)]
        while True:
            nxt=A[path[-1]]-set(path)
            if not nxt:break
            path.append(min(nxt))
        need(set(path)==set(C),'bad forest traversal');out.append(path)
    return out

def solve(V,E,enumerate_opt=False):
    """Exact edge-subset optimization by vertex elimination and connectivity DP.
    Each edge is decided once at its first endpoint. The state records degrees
    and connectivity via previously processed vertices. Reject a choice joining
    vertices already in one component. Only acyclic max-degree-two sets remain.
    """
    V=tuple(sorted(V)); E=tuple(sorted(tuple(sorted(e)) for e in E))
    A={v:set() for v in V}
    for a,b in E:A[a].add(b);A[b].add(a)
    order=[];remaining=set(V);boundary=set()
    while remaining:
        def key(v):
            R=remaining-{v};F=(boundary|A[v])&R
            return len(F),len(A[v]&R),v
        v=min(remaining,key=key);order.append(v);remaining.remove(v);boundary=(boundary|A[v])&remaining
    ix={v:i for i,v in enumerate(order)}
    later={i:tuple(sorted(ix[w] for w in A[v] if ix[w]>i)) for i,v in enumerate(order)}
    def canonical(labels):
        d={};return tuple(d.setdefault(x,len(d)) for x in labels)
    def choices(v,degs,labels):
        eligible=[j-v for j in later[v] if degs[j-v]<2]
        for size in range(min(2-degs[0],len(eligible))+1):
            for chosen in combinations(eligible,size):
                joined=[labels[0]]+[labels[i] for i in chosen]
                if len(set(joined))!=len(joined):continue
                D=list(degs);L=list(labels)
                for i in chosen:D[i]+=1
                for i in range(len(L)):
                    if L[i] in joined:L[i]=labels[0]
                yield chosen,(v+1,tuple(D[1:]),canonical(L[1:]))
    @lru_cache(None)
    def opt(v,degs,labels):
        if not degs:return 0,1
        best=-1;count=0
        for chosen,state in choices(v,degs,labels):
            score,cnt=opt(*state);score+=len(chosen)
            if score>best:best,count=score,cnt
            elif score==best:count+=cnt
        return best,count
    initial=(0,(0,)*len(V),tuple(range(len(V))));maximum,count=opt(*initial)
    def walk(state,edges):
        if not state[1]:yield tuple(sorted(tuple(sorted(e)) for e in edges));return
        for chosen,after in choices(*state):
            if opt(*after)[0]+len(chosen)==opt(*state)[0]:
                extra=[(order[state[0]],order[state[0]+i]) for i in chosen]
                yield from walk(after,edges+extra)
    witness=next(walk(initial,[]))
    out=dict(minimum_path_cover=len(V)-maximum,optimal_cover_count=count,
             path_cover=to_paths(V,witness),selected_edges=list(witness),dp_states=opt.cache_info().currsize)
    if enumerate_opt:
        covers=list(walk(initial,[]));need(len(covers)==len(set(covers))==count,'cover enumeration count')
        out['all_optimal_path_covers']=[to_paths(V,es) for es in covers]
    return out

def brute(V,E):
    best=-1;covers=[]
    for mask in range(1<<len(E)):
        es=[e for i,e in enumerate(E) if mask>>i&1]
        if not forest_ok(V,es):continue
        if len(es)>best:best=len(es);covers=[es]
        elif len(es)==best:covers.append(es)
    return len(V)-best,covers

def certify_prime_pack(doc):
    certs=doc['primality_certificates'];known=set()
    for key in sorted(certs,key=int):
        N=int(key);c=certs[key]
        need(N>=2,'invalid integer')
        if N==2:need(c=={'base':1,'factors':[]},'base case');known.add(2);continue
        a=c['base'];fac=c['factors'];qs=[q for q,e in fac]
        need(type(a) is int and 1<a<N,'Lucas base domain')
        need(bool(fac) and len(set(qs))==len(qs),'factor list')
        need(all(type(q) is int and type(e) is int and q in known and e>0 for q,e in fac),'recursive factors')
        need(prod(q**e for q,e in fac)==N-1,'full factorization')
        need(pow(a,N-1,N)==1,'Fermat condition')
        need(all(gcd(pow(a,(N-1)//q,N)-1,N)==1 for q in qs),'order conditions')
        known.add(N)
    primes=[n for n in range(2,74) if all(n%d for d in range(2,isqrt(n)+1))]
    P=1;rows=[];tests=Counter()
    need(len(doc['rows'])==19,'19 primorial rows')
    for i,r in enumerate(doc['rows']):
        P*=primes[i];need((r['P'],r['n'],r['p'])==(P,i+1,primes[i]),'primorial row')
        for sign in ('plus','minus'):
            if i==0 and sign=='minus':need(r[sign] is None,'undefined L1');continue
            d0=r[sign]['offset'];ws=r[sign]['witnesses']
            need(type(d0) is int and d0>=2 and len(ws)==d0-1,'minimum interval length')
            for d,w in enumerate(ws,2):
                N=P+d if sign=='plus' else P-d;tests[sign]+=1
                if d==d0:need(w=={'prime':True} and N in known,'terminal exact primality');continue
                need(not w.get('prime',False),'premature prime')
                if 'divisor' in w:
                    q=w['divisor'];need(type(q) is int and 1<q<N and N%q==0,'proper divisor')
                elif 'mr_base' in w:
                    a=w['mr_base'];need(type(a) is int and 1<a<N and N>3 and N%2,'MR base domain')
                    k=N-1;s=0
                    while not k%2:k//=2;s+=1
                    chain=[pow(a,k*(2**j),N) for j in range(s)]
                    need(chain[0]!=1 and N-1 not in chain,'MR witness must fail')
                else:raise ValueError('missing exact compositeness witness')
        F=r['plus']['offset'];q=primes[i+1]
        need(F>=q,'wheel lower bound')
        rows.append(dict(n=i+1,F=F,q=q,P=P,bit=F==q,raw_gap=F-q,literal_gap=F-primes[i]-1))
    return rows,dict(all_lucas_certificates=len(known),positive_endpoint_tests=tests['plus'],negative_endpoint_tests=tests['minus'])

def square_graph(V):
    # Generate one reflection matching at a time, rather than pairwise isqrt.
    V=set(V);E=set()
    if not V:return V,E
    for k in range(1,isqrt(2*max(V))+1):
        for a in V:
            b=k*k-a
            if b in V and a<b:E.add((a,b))
    return V,E

def rle(xs):return [dict(value=v,length=len(list(g))) for v,g in groupby(xs)]

def check_candidate(data,arithmetic,graphs):
    need(len(data['prime_rows'])==19,'arithmetic output length')
    for a,r in zip(arithmetic,data['prime_rows']):
        need((r['n'],r['fortunate'],r['p_next'],r['primorial'],r['candidate_prime'],r['next_prime_gap'],r['literal_prime_plus_one_gap'])==(a['n'],a['F'],a['q'],a['P'],a['bit'],a['raw_gap'],a['literal_gap']),'arithmetic output disagreement')
        N=a['P']+a['q'];need(r['candidate']==N,'candidate integer')
        if not a['bit']:
            q=r['candidate_proper_factor'];need(type(q) is int and 1<q<N and N%q==0 and r['candidate_cofactor']==N//q,'first candidate factor')
    need(len(data['graph_rows'])==28,'graph output length')
    for r,a in zip(data['graph_rows'],graphs):
        n=a['N'];V,E=square_graph(range(1,n+1));C=parts(V,E)
        for key in ('N','edges','components','minimum_path_cover','path_cover_deficit','global_traceable','componentwise_traceable','hamilton_path_count'):
            need(r[key]==a[key],'graph output disagreement '+str(n)+' '+key)
        need(r['component_vertex_sets']==C,'components')
        degree=Counter(x for e in E for x in e);born=sorted(v for v in V if tuple(sorted((v,n))) in E)
        oldV,oldE=square_graph(range(1,n));oldC=parts(oldV,oldE)
        touched=sum(bool(set(part)&set(born)) for part in oldC)
        need(r['birth_neighbors']==born and r['birth_degree']==len(born)==isqrt(2*n-1)-isqrt(n),'birth neighbors')
        need(r['birth_touched_components']==touched,'touched components')
        need(r['cycle_rank']==len(E)-n+len(C),'cycle rank')
        need(r['leaves']==sorted(v for v in V if degree[v]==1),'leaves')
        need(r['branch_vertices']==sorted(v for v in V if degree[v]>=3),'branches')
        cc=[]
        for part in C:
            answer=solve(part,[e for e in E if set(e)<=set(part)])
            cc.append(answer['optimal_cover_count'] if answer['minimum_path_cover']==1 else 0)
        need(r['component_hamilton_path_counts']==cc,'component path counts')
        hp=data['hamilton_paths'][str(n)];need(len(hp)==a['hamilton_path_count'],'all-path count')
        canonical_paths=set()
        for path in hp:
            need(len(path)==n and set(path)==V,'path vertex set')
            es={tuple(sorted(e)) for e in zip(path,path[1:])};need(es<=E and len(es)==n-1,'path adjacency')
            canonical_paths.add(tuple(min(path,path[::-1])))
        need(len(canonical_paths)==len(hp),'duplicate reverse paths')
        paths=r['path_cover'];flat=[v for p in paths for v in p]
        need(len(flat)==n and set(flat)==V and len(paths)==a['minimum_path_cover'],'minimum cover partition')
        need(all(p and all(tuple(sorted(e)) in E for e in zip(p,p[1:])) for p in paths),'minimum cover adjacency')
    bits=[x['bit'] for x in arithmetic[:17]];graphbits=[x['componentwise_traceable'] for x in graphs[9:26]]
    need(bits==graphbits==data['exact_word'],'finite matching word')
    need([r['length'] for r in rle(bits)]==data['run_lengths']==[3,2,3,5,1,1,2],'finite run lengths')
    need(arithmetic[17]['bit']!=graphs[26]['componentwise_traceable'],'first extension failure')
    need(data['first_extension_failure']['F']==arithmetic[17]['F']==109,'first failure F')
    need(data['first_extension_failure']['graph_hamilton_path_count']==graphs[26]['hamilton_path_count']==35,'first failure graph')


def main():
    manifest=json.loads((SNAP/'MANIFEST.json').read_text())
    for f in manifest['files']:
        b=(SNAP/f['path']).read_bytes();need(len(b)==f['bytes'] and hashlib.sha256(b).hexdigest()==f['sha256'],'snapshot manifest')
    doc=json.loads((SNAP/'inputs/prime_offsets_certificates.json').read_text())
    arithmetic,prime_counts=certify_prime_pack(doc)
    graphs=[]
    for n in range(1,29):
        V,E=square_graph(range(1,n+1));result=solve(V,E,n==24);c=len(parts(V,E));pc=result['minimum_path_cover']
        row=dict(N=n,edges=len(E),components=c,path_cover_deficit=pc-c,
                 global_traceable=pc==1,componentwise_traceable=pc==c,
                 hamilton_path_count=result['optimal_cover_count'] if pc==1 else 0,**result)
        graphs.append(row)
    data=json.loads((SNAP/'results.json').read_text());check_candidate(data,arithmetic,graphs)
    # Cross-check the independent engine against literal edge-mask enumeration
    # on every simple labelled graph through five vertices.
    controls=0;births=0
    for n in range(0,6):
        V=set(range(n));full=list(combinations(range(n),2))
        for mask in range(1<<len(full)):
            E=[e for i,e in enumerate(full) if mask>>i&1];pc,covers=brute(V,E);dp=solve(V,E)
            need((pc,len(covers))==(dp['minimum_path_cover'],dp['optimal_cover_count']),'brute cross-check');controls+=1
            if n>4:continue
            c=len(parts(V,E))
            for neighbors in range(1<<n):
                W={v for v in V if neighbors>>v&1};EE=E+[(v,n) for v in W]
                pc2=solve(V|{n},EE)['minimum_path_cover'];cc=len(parts(V|{n},EE));r=sum(bool(set(C)&W) for C in parts(V,E))
                available=False
                for cover in covers:
                    D=Counter(x for e in cover for x in e);C=parts(V,cover)
                    hits=sum(any(v in W and D[v]<=1 for v in part) for part in C)
                    if hits>=2:available=True;break
                need(pc-1<=pc2<=pc+1 and (pc2==pc-1)==available,'birth endpoint iff')
                need(cc==c+1-r and pc2-cc-(pc-c)==pc2-pc-1+r,'deficit birth identity');births+=1
    mutants=[]
    def rejected(label,fun):
        try:fun()
        except (ValueError,KeyError,TypeError,IndexError):mutants.append(label);return
        raise RuntimeError('accepted semantic mutation '+label)
    for label,change in [
        ('ordinary/componentwise substitution',lambda d:d['graph_rows'][9].__setitem__('componentwise_traceable',False)),
        ('invented Hamilton path at24',lambda d:d['graph_rows'][23].__setitem__('hamilton_path_count',1)),
        ('wrong minimum cover at24',lambda d:d['graph_rows'][23].__setitem__('minimum_path_cover',1)),
        ('suppressed first extension failure',lambda d:d['first_extension_failure'].__setitem__('F',67)),
        ('duplicated path vertex',lambda d:d['hamilton_paths']['27'][0].__setitem__(1,18)),
        ('raw gaps substituted for bits',lambda d:d.__setitem__('exact_word',[r['next_prime_gap'] for r in d['prime_rows'][:17]])),
        ('false composite candidate factor',lambda d:d['prime_rows'][17].__setitem__('candidate_proper_factor',67)),
        ('forged wheel equality',lambda d:d['prime_rows'][17].__setitem__('candidate_prime',True))]:
        bad=deepcopy(data);change(bad);rejected(label,lambda bad=bad:check_candidate(bad,arithmetic,graphs))
    for label,change in [
        ('invalid Lucas base',lambda d:d['primality_certificates']['5'].__setitem__('base',1)),
        ('incomplete Lucas factorization',lambda d:d['primality_certificates']['5'].__setitem__('factors',[[2,1]])),
        ('missing minimality interval',lambda d:d['rows'][17]['plus']['witnesses'].pop()),
        ('shifted primorial',lambda d:d['rows'][0].__setitem__('P',3)),
        ('MR base equals endpoint',lambda d:d['rows'][3]['plus']['witnesses'].__setitem__(9,{'mr_base':221}))]:
        bad=deepcopy(doc);change(bad);rejected(label,lambda bad=bad:certify_prime_pack(bad))
    need(not forest_ok({0,1,2},[(0,1),(1,2),(0,2)]),'triangle cycle accepted as cover')
    need(solve({0,1,2},[(0,1),(1,2),(0,2)])['optimal_cover_count']==3,'triangle endpoint hostile')
    out=dict(status='PASS',method='Independent edge-subset connectivity DP; no producer imports',
             candidate_manifest_sha256=hashlib.sha256((SNAP/'MANIFEST.json').read_bytes()).hexdigest(),
             candidate_results_sha256=hashlib.sha256((SNAP/'results.json').read_bytes()).hexdigest(),
             arithmetic=arithmetic,prime_certificate_checks=prime_counts,graphs=graphs,
             complete_brute_graph_controls=controls,complete_vertex_birth_controls=births,
             semantic_mutations_rejected=mutants,
             prime_word=rle([r['bit'] for r in arithmetic[:17]]),
             graph_word=rle([r['componentwise_traceable'] for r in graphs[9:26]]),
             full_componentwise_word=rle([r['componentwise_traceable'] for r in graphs]))
    body=canon(out);(ROOT/'independent_results.json').write_text(body)
    print(canon(dict(status='PASS',candidate_manifest_sha256=out['candidate_manifest_sha256'],
                    independent_results_sha256=hashlib.sha256(body.encode()).hexdigest(),
                    prime_certificate_checks=prime_counts,brute_graph_controls=controls,
                    vertex_birth_controls=births,mutations_rejected=len(mutants),
                    minimum_cover_counts={r['N']:r['optimal_cover_count'] for r in graphs},
                    all_24_minimum_covers=len(graphs[23]['all_optimal_path_covers']))))
if __name__=='__main__':main()
