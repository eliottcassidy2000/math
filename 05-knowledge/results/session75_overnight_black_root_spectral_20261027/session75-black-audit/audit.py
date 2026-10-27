#!/usr/bin/env python3
"""Independent Session75 audit; stdlib, explicit failures, no producer imports.
Enumeration is a Gray walk through every labelled tournament. Moves use sets
of ordered arcs; anti-automorphisms are checked on every ordered arc.
"""
import hashlib, itertools, json, pathlib, sys
ROOT=pathlib.Path(__file__).resolve().parent

def need(ok, message):
    if not ok: raise RuntimeError(message)

def encode(E,n): return [sum(1<<v for u,v in E if u==i) for i in range(n)]
def decode(rows):
    n=len(rows); E=frozenset((u,v) for u in range(n) for v in range(n) if rows[u]>>v&1)
    need(all(((u,v) in E)+((v,u) in E)==1 for u in range(n) for v in range(u+1,n)), 'not tournament')
    need(not any(u==v for u,v in E),'loop')
    return E

def degrees(E,n): return [sum((u,v) in E for v in range(n)) for u in range(n)]
def converse(E): return frozenset((v,u) for u,v in E)
def reach(E,S,start):
    seen={start}
    while True:
        more=seen | {v for u,v in E if u in seen and v in S}
        if more==seen: return seen
        seen=more

def strong(E,S):
    if not S:return False
    root=min(S)
    return reach(E,S,root)==S and reach(converse(E),S,root)==S

def target(E,n):
    S=set(range(n));d=degrees(E,n)
    return [v for v in S if d[v] in (0,n-1) and strong(E,S-{v})]

def paths(E,n,end=None,closed=False):
    neigh={u:tuple(v for v in range(n) if (u,v) in E) for u in range(n)}
    def grow(p,remaining):
        if not remaining:
            if (end is None or p[-1]==end) and (not closed or (p[-1],p[0]) in E):yield tuple(p)
            return
        for v in neigh[p[-1]]:
            if v in remaining and (v!=end or len(remaining)==1):
                yield from grow(p+[v],remaining-{v})
    for v in range(n):
        if end is None or v!=end:yield from grow([v],set(range(n))-{v})

def anti_failure(E,p):
    r=dict(zip(p,p[::-1]))
    return next(((u,v,r[v],r[u]) for u,v in sorted(E) if (r[v],r[u]) not in E), None)

def reversal(E,p):
    B=set(E)
    for u,v in zip(p,p[1:]):
        need((u,v) in B,'R path not directed')
        B.remove((u,v)); B.add((v,u))
    return frozenset(B)

def complement_nonpath(E,p):
    pairs={frozenset((u,v)) for u,v in zip(p,p[1:])}
    return frozenset((u,v) if frozenset((u,v)) in pairs else (v,u) for u,v in E)

def cuts(E,n):
    return [sum(bool(s>>u&1) and not bool(s>>v&1) for u,v in E) for s in range(1,(1<<n)-1)]

def inspect_move(E,p,n,black=True):
    need(sorted(p)==list(range(n)),'labels')
    need(all((u,v) in E for u,v in zip(p,p[1:])),'directedness')
    failure=anti_failure(E,p)
    if black:need(failure is not None,'blue path supplied as black')
    F=reversal(E,p)
    need((anti_failure(F,p[::-1]) is None)==(failure is None),'R colour change')
    d,e=degrees(E,n),degrees(F,n)
    need(all(e[v]-d[v]==(v==p[-1])-(v==p[0]) for v in range(n)),'score current')
    a,b=cuts(E,n),cuts(F,n)
    need(all(b[s-1]-a[s-1]==((s>>p[-1]&1)-(s>>p[0]&1)) for s in range(1,(1<<n)-1)),'cut current')
    K=complement_nonpath(E,p)
    need(K==converse(F),'K/R identity')
    need(all((u,v) in K for u,v in zip(p,p[1:])),'K word lost')
    need((anti_failure(K,p) is None)==(failure is None),'K colour change')
    return F

def walk_gray(n,balanced=False):
    pairs=list(itertools.combinations(range(n),2)); E=set((v,u) for u,v in pairs)
    d=list(range(n)); old=0
    lo=(n-1)//2;hi=n//2
    for i in range(1<<len(pairs)):
        if i:
            g=i^(i>>1); change=g^old; k=change.bit_length()-1;u,v=pairs[k]
            if g>>k&1:E.remove((v,u));E.add((u,v));d[u]+=1;d[v]-=1
            else:E.remove((u,v));E.add((v,u));d[u]-=1;d[v]+=1
            old=g
        if not balanced or all(lo<=x<=hi for x in d):yield frozenset(E)

def construct(E,n):
    # A separate path search, not the producer's certificate or implementation.
    d=degrees(E,n)
    p=next((p for p in paths(E,n,closed=True) if d[p[-1]]>=d[p[0]] and anti_failure(E,p)),None)
    need(p is not None, 'no descending black closed path')
    dual=n%2==0 and d[p[-1]]!=n//2
    if dual:E=converse(E);p=p[::-1]
    v=p[-1];start=p[0];route=[p];F=reversal(E,p)
    if n>=7:need(strong(F,set(range(n))-{v}),'deleted core after first move not strong')
    for step in range(1,(n-1)//2):
        q=next((q for q in paths(F,n,end=v) if step!=1 or q[0]!=start),None)
        need(q is not None,'named second start / fixed endpoint missing')
        need(anti_failure(F,q) is not None,'chosen extension is blue')
        route.append(q);F=reversal(F,q)
    return [p[::-1] if dual else p for p in route],dual

def inspect_route(E,route,n,dual=False):
    original=E;Kstate=E; named=[]; records=[]
    for j,p in enumerate(route):
        need(strong(E,set(range(n))),'route intermediate not strong')
        q=p[::-1] if j%2 else p
        expected_source=converse(E) if j%2 else E
        need(Kstate==expected_source,'literal K gauge start')
        F=inspect_move(E,p,n)
        need(anti_failure(Kstate,q) is not None,'literal K blue')
        need(all((u,v) in Kstate for u,v in zip(q,q[1:])),'literal K directedness')
        Kstate=complement_nonpath(Kstate,q)
        need(Kstate==(converse(F) if (j+1)%2 else F),'literal K target gauge')
        records.append({'R_word':list(p),'K_word':list(q),'failed_anti_R':anti_failure(E,p),'R_target_rows':encode(F,n),'K_target_rows':encode(Kstate,n)})
        named.append(p[-1] if dual else p[0]);E=F
    d=degrees(original,n);lower=min(min(s,n-1-s) for s in d)
    need(len(route)==lower,'semidegree bound not tight')
    need(target(E,n) and target(Kstate,n),'final core not strong')
    if len(route)>=2:need(named[0]!=named[1],'two named starts lost')
    return {'initial_rows':encode(original,n),'dual':dual,'distance':len(route),'terminal_extreme':target(E,n),'moves':records}

def binding_check():
    b=json.loads((ROOT/'SOURCE_BINDINGS.json').read_text())
    for x in b['sources']:
        need(hashlib.sha256((ROOT/'sources'/x['path']).read_bytes()).hexdigest()==x['sha256'],'source binding '+x['path'])
    need((ROOT/'sources/producer/normal.json').read_bytes()==(ROOT/'sources/producer/optimized.json').read_bytes(),'producer normal/-O differ')
    need((ROOT/'sources/producer/audit-normal.json').read_bytes()==(ROOT/'sources/producer/audit-optimized.json').read_bytes(),'producer independent normal/-O differ')
    return {'pin':b['repository_pin'],'bound_files':len(b['sources']),'bindings_sha256':hashlib.sha256((ROOT/'SOURCE_BINDINGS.json').read_bytes()).hexdigest()}

def allcolour_route(E,n):
    cut=cuts(E,n);lam=min(cut);s=cut.index(lam)+1;current=E;ps=[]
    for j in range(lam):
        need(strong(current,set(range(n))),'allcolour route lost strong too early')
        # Search a literal cycle cut from A to Ac. No colour filter.
        p=next(p for p in paths(current,n,closed=True) if s>>p[0]&1 and not s>>p[-1]&1)
        current=reversal(current,p);ps.append(p)
    need(not strong(current,set(range(n))) and cuts(current,n)[s-1]==0,'allcolour route failed')
    return lam

def join(blocks):
    E=set();off=0;n=sum(len(x) for x in blocks)
    for rows in blocks:
        E|={(u+off,v+off) for u,v in decode(rows)}
        E|={(u,v) for u in range(off,off+len(rows)) for v in range(off+len(rows),n)}
        off+=len(rows)
    return frozenset(E)

def main():
    report={'status':'PASS','source_bindings':binding_check(),'independence':'No producer imports. Gray enumeration; ordered-arc operations; explicit anti-map on all arcs; graph reachability for target.'}
    balanced={}; total_moves=0
    for n in (4,6,7):
        count=0;digest=hashlib.sha256();dualcount=0;first=None
        for E in walk_gray(n,balanced=True):
            need(strong(E,set(range(n))),'balanced input not strong')
            route,dual=construct(E,n); rec=inspect_route(E,route,n,dual)
            if first is None:first=rec
            digest.update(json.dumps(rec,sort_keys=True).encode()+b'\n');count+=1;dualcount+=dual;total_moves+=len(route)
        need(count=={4:24,6:2640,7:2640}[n],'balanced universe count')
        balanced[n]={'labelled_count':count,'witness_digest':digest.hexdigest(),'converse_selected_count':dualcount,'first_literal_witness':first}
    report['complete_balanced_universes']=balanced
    report['balanced_R_and_K_moves_checked']=total_moves
    # Challenge inherited descending theorem on EVERY labelled strong n<=6.
    descent={};lambda_distributions={}
    for n in (3,4,5,6):
        cnt=blueonly=0;hist={}
        for E in walk_gray(n):
            if not strong(E,set(range(n))):continue
            cnt+=1;d=degrees(E,n)
            witness=next((p for p in paths(E,n,closed=True) if d[p[-1]]>=d[p[0]] and anti_failure(E,p)),None)
            if witness is None:
                need(all(anti_failure(E,p) is None for p in paths(E,n)),'black-active descending lemma failure')
                blueonly+=1
            lam=allcolour_route(E,n);hist[lam]=hist.get(lam,0)+1
        descent[n]={'strong_labelled':cnt,'pure_blue_exceptions':blueonly,'descending_black_witnesses':cnt-blueonly}
        lambda_distributions[n]=hist
    report['complete_strong_descending_check']=descent
    report['complete_allcolour_lambda_check']=lambda_distributions
    # Replay all saved producer routes, but independently recompute R and K.
    prod=json.loads((ROOT/'sources/producer/normal.json').read_text());pc=0
    for rec in prod['seeds']+prod['examples']+[prod['nonregular_target_separation']]:
        E=decode(rec['initial_rows']);n=len(rec['initial_rows']);route=[tuple(m['R_path']) for m in rec['moves']]
        got=inspect_route(E,route,n,rec['dual_for_construction'])
        for j,(m,g) in enumerate(zip(rec['moves'],got['moves'])):
            need(m['K_path']==g['K_word'],'producer K word')
            need(m['K_target_rows']==g['K_target_rows'],'producer K rows')
            need(m['K_source_rows']==(rec['initial_rows'] if j==0 else rec['moves'][j-1]['K_target_rows']),'producer K chain')
        need(encode(decode(rec['terminal_R_rows']),n)==got['moves'][-1]['R_target_rows'],'producer R terminal')
        pc+=len(route)
    report['producer_literal_moves_independently_checked']=pc
    # Repeated named-start hostile: actual anti maps and exact SCC shape.
    h=prod['same_start_hostile'];E=decode(h['initial_rows']);colors=[]
    for p in h['R_paths']:
        colors.append('blue' if anti_failure(E,p) is None else 'black'); E=inspect_move(E,p,7,black=False)
    d=degrees(E,7);src=[i for i in range(7) if d[i]==6];sink=[i for i in range(7) if d[i]==0];middle=set(range(7))-set(src)-set(sink)
    need(colors==['black','blue','blue'] and src==[0] and sink==[6] and len(middle)==5,'hostile colours/extremes')
    need(strong(E,middle) and all(sum((u,v) in E for v in middle)==2 for u in middle),'middle is not R5')
    need(not target(E,7),'same-start hostile accidentally core target')
    report['same_start_hostile']={'colors':colors,'source':src,'sink':sink,'middle':sorted(middle),'middle_regular_strong':True,'terminal_rows':encode(E,7)}
    # n10: exhaust ALL Hamilton words at depth 1; enumerate cuts and exact shape.
    r=prod['nonregular_target_separation'];E=decode(r['initial_rows']);n=10;d=degrees(E,n)
    need(strong(E,set(range(n))) and min(cuts(E,n))==1 and min(min(x,n-1-x) for x in d)==2,'n10 source invariants')
    first_exits=[];hpc=bc=0;core1=0
    for p in paths(E,n):
        hpc+=1;black=anti_failure(E,p) is not None;bc+=black
        F=reversal(E,p)
        if target(F,n):core1+=1
        if black and not strong(F,set(range(n))):
            if len(first_exits)<1:first_exits.append({'R_word':list(p),'failed_anti':anti_failure(E,p),'terminal_rows':encode(F,n)})
    need(core1==0 and first_exits,'n10 target separation')
    report['n10_separation']={'scores':d,'lambda':1,'delta':2,'hamilton_words':hpc,'black_words':bc,'one_step_core_targets_all_colours':core1,'black_nonstrong_exit':first_exits[0],'black_core_distance':2}
    controls=[]
    for h in prod['pure_blue_controls']:
        E=decode(h['rows']);n=len(h['rows']);ps=[p for p in itertools.permutations(range(n)) if all((u,v) in E for u,v in zip(p,p[1:]))]
        need(len(ps)==h['hamilton_paths'] and all(anti_failure(E,p) is None for p in ps),'pure blue permutations')
        controls.append({'name':h['name'],'hamilton_words':len(ps),'black_words':0})
    report['pure_blue_all_permutation_controls']=controls
    report['checker_sha256']=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()
    print(json.dumps(report,sort_keys=True,indent=2))

if __name__=='__main__':main()
