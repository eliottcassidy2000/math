#!/usr/bin/env python3
"""Full-arc converse checks on the exhaustive inherited positive-error slice.
No producer imports and no unrestricted order-ten tournament enumeration.
"""
from collections import Counter, defaultdict
from itertools import combinations, product, permutations
from pathlib import Path
import argparse, hashlib, json, math

HERE=Path(__file__).resolve().parent
CHECKS=Counter()
def need(ok,msg):
    CHECKS[msg]+=1
    if not ok: raise RuntimeError(msg)
def validate(T):
    n=len(T)
    need(all(0<=r<1<<n and not(r>>i&1) for i,r in enumerate(T)),'no loops or external bits')
    need(all(((T[i]>>j)^(T[j]>>i))&1 for i,j in combinations(range(n),2)),'exactly one arc per pair')
def reverse(T): return tuple(((1<<len(T))-1)^r^(1<<i) for i,r in enumerate(T))
def colors(T):
    c=(0,)*len(T)
    while True:
        sig=tuple((c[i],)+tuple(sum((r>>j)&1 for j in range(len(T)) if c[j]==k) for k in range(max(c)+1)) for i,r in enumerate(T))
        index={s:k for k,s in enumerate(sorted(set(sig)))}
        nxt=tuple(index[s] for s in sig)
        if len(set(nxt))==len(set(c)): return nxt
        c=nxt
def cells(T):
    c=colors(T)
    return tuple(tuple(i for i,x in enumerate(c) if x==k) for k in range(max(c)+1))
def mappings(A,B,first=False):
    """Every exact directed isomorphism; colors merely prune the full search."""
    if len(A)!=len(B): return []
    ca,cb=colors(A),colors(B)
    if Counter(ca)!=Counter(cb): return []
    candidates={v:tuple(w for w in range(len(B)) if cb[w]==ca[v]) for v in range(len(A))}
    order=sorted(range(len(A)),key=lambda v:(len(candidates[v]),-A[v].bit_count(),v))
    assigned={}; found=[]
    def dfs(k,used):
        if k==len(A):
            p=tuple(assigned[v] for v in range(len(A)))
            need(all(((A[i]>>j)&1)==((B[p[i]]>>p[j])&1) for i in range(len(A)) for j in range(len(A))),'completed full arc map')
            found.append(p);return first
        v=order[k]
        for w in candidates[v]:
            if w in used:continue
            if any(((A[v]>>u)&1)!=((B[w]>>z)&1) for u,z in assigned.items()):continue
            assigned[v]=w
            if dfs(k+1,used|{w}): return True
            del assigned[v]
        return False
    dfs(0,set());return found
def error(T):
    group=mappings(T,T);sizes={v:len(C) for C in cells(T) for v in C};todo=set(range(len(T)));orbits=[]
    while todo:
        v=min(todo);orb={g[v] for g in group};need(orb<=todo,'automorphism orbit partition');todo-=orb;orbits.append(sorted(orb))
    return sum(sizes[o[0]]-len(o) for o in orbits),group,orbits
def strong(T):
    for start in range(len(T)):
        seen={start};queue=[start]
        while queue:
            v=queue.pop()
            for w in range(len(T)):
                if T[v]>>w&1 and w not in seen: seen.add(w);queue.append(w)
        if len(seen)!=len(T):return False
    return True
def ordinal(*parts):
    out=[];offset=0;n=sum(map(len,parts))
    for P in parts:
        tail=((1<<n)-1)^((1<<(offset+len(P)))-1)
        out.extend((r<<offset)|tail for r in P);offset+=len(P)
    return tuple(out)
def cyc(n):return tuple(sum(1<<((v+d)%n) for d in range(1,(n+1)//2)) for v in range(n))
def substitute(Q,parts):
    starts=[];offset=0;out=[]
    for P in parts:starts.append(offset);out.extend(r<<offset for r in P);offset+=len(P)
    for i,j in combinations(range(len(parts)),2):
        a,b=(i,j) if Q[i]>>j&1 else (j,i)
        mask=((1<<len(parts[b]))-1)<<starts[b]
        for v in range(starts[a],starts[a]+len(parts[a])):out[v]|=mask
    return tuple(out)
def triangles9(signs):
    out=[r<<(3*i) for i in range(3) for r in cyc(3)]
    for (i,j),sign in zip(combinations(range(3),2),signs):
        for x in range(3):
            for y in range(3):
                a,b=(3*i+x,3*j+y) if y==(sign*x)%3 else (3*j+y,3*i+x)
                out[a]|=1<<b
    return tuple(out)
def permutation_power(p,k):
    q=tuple(range(len(p)))
    for _ in range(k):q=tuple(p[x] for x in q)
    return q
def permutation_order(p):
    seen=set();ans=1
    for v in range(len(p)):
        if v in seen:continue
        k=0;u=v
        while u not in seen:seen.add(u);k+=1;u=p[u]
        ans=math.lcm(ans,k)
    return ans

def anti_certificate(T,p):
    squared=tuple(p[p[i]] for i in range(len(p)));m=permutation_order(squared)
    need(m%2==1,'square of anti-map has odd order')
    q=permutation_power(p,m)
    need(all(q[q[v]]==v for v in range(len(T))),'odd power gives involution')
    need(all(((T[i]>>j)&1)==((T[q[j]]>>q[i])&1) for i,j in combinations(range(len(T)),2)),'involutory power is anti-map')
    C=cells(T);lookup={frozenset(c):i for i,c in enumerate(C)}
    cp=tuple(lookup[frozenset(q[v] for v in c)] for c in C)
    fixed_vertices=[v for v in range(len(T)) if q[v]==v]
    fixed_cells=[i for i in range(len(C)) if cp[i]==i]
    need(len(fixed_vertices)==len(T)%2,'anti-involution vertex fixed count')
    need(len(fixed_cells)==len(T)%2,'anti-involution stable-cell fixed count')
    return {'anti_map':p,'anti_involution':q,'anti_map_order':permutation_order(p),'stable_cell_permutation':cp,'fixed_vertices':fixed_vertices,'fixed_cells':fixed_cells}

def load(name):return json.loads((HERE/'sources'/name).read_text())
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',default=str(HERE/'normal.json'));args=ap.parse_args()
    pins=load('pins.json')
    for rec in pins:
        need(hashlib.sha256((HERE/'sources'/rec['file']).read_bytes()).hexdigest()==rec['sha256'],'inherited input SHA256')
    inherited=load('order-ten-audit-receipt.json')
    records=[{'source':'strong_'+r['sector'],'degree':r.get('r'),'rows':tuple(r['rows']),'expected_error':r['error']} for r in inherited['representatives'] if r['error']]
    need(len(records)==120,'exhaustive strong ten positive slice size')
    BAD7=(74,84,97,22,37,11,56)
    regular=load('regular-nine-receipt.json');regular9=[]
    for component in regular['components']:
        P=regular['states'][component['states'][0]]['parent']
        T=tuple(r|((1<<8) if r.bit_count()==3 else 0) for r in P)+(sum(1<<v for v,r in enumerate(P) if r.bit_count()==4),)
        validate(T);e,g,o=error(T)
        need(e==component['error'],'regular nine full automorphism error agrees with parent-reroot receipt')
        need(all(r.bit_count()==4 for r in T),'regular nine unique extension')
        regular9.append((T,e))
    need(len(regular9)==15,'regular nine input has fifteen classes')
    need(all(not mappings(A,B,True) for i,(A,_) in enumerate(regular9) for B,_ in regular9[:i]),'regular nine representatives independently distinct')
    strong9=[(T,e,'regular9') for T,e in regular9 if e]
    strong9.append((substitute(cyc(3),(BAD7,(0,),(0,))),14,'seven9'))
    for signs in product((-1,1),repeat=3):
        if math.prod(signs)==-1:strong9.append((triangles9(signs),18,'triangles9'))
    need(len(strong9)==17,'seventeen positive strong nine classes')
    for T,e,source in strong9:
        need(strong(T),'all reconstructed nine cores strong');need(error(T)[0]==e,'all reconstructed nine errors exact')
        records.extend([{'source':'left_singleton_'+source,'rows':ordinal((0,),T),'expected_error':e},{'source':'right_singleton_'+source,'rows':ordinal(T,(0,)),'expected_error':e}])
    small={0:[()],1:[(0,)],2:[(2,0)],3:[(6,4,0),cyc(3)]}
    for left in range(4):
        for L in small[left]:
            for R in small[3-left]:records.append({'source':'seven_with_three_padding','rows':ordinal(L,BAD7,R),'expected_error':14})
    need(len(records)==160,'all 160 positive-error ordinary representatives')
    hist=Counter();sc=Counter();strong_hist=Counter();sc_by_sector=Counter();certificates=[]
    for i,r in enumerate(records):
        T=r['rows'];validate(T);e,g,orbits=error(T)
        need(e==r['expected_error'],'every full ten error agrees with inherited classification')
        need({frozenset(C) for C in cells(T)}=={frozenset(C) for C in cells(reverse(T))},'actual stable partitions invariant under converse')
        need(all(len(C)%2 for C in cells(T)),'actual stable cells odd')
        need(error(reverse(T))[0]==e,'error invariant under converse')
        anti=mappings(T,reverse(T));hist[e]+=1
        if strong(T):strong_hist[e]+=1
        cert={'index':i,'source':r['source'],'degree':r.get('degree'),'rows':T,'error':e,'stable_cells':cells(T),'automorphism_order':len(g),'vertex_orbits':orbits,'self_converse':bool(anti)}
        if anti:
            sc[e]+=1;sc_by_sector[(r['source'],r.get('degree'))]+=1
            need(len(anti)==len(g),'anti-isomorphism coset same size as automorphism group')
            cert['anti_certificates']=[anti_certificate(T,p) for p in anti]
        certificates.append(cert)
    # Compare actual representatives against every possible converse partner.
    partners=[];isomorphism_pairs=0
    for i,r in enumerate(records):
        T=r['rows'];candidate=[]
        for j,s in enumerate(records):
            if s['expected_error']!=r['expected_error']:continue
            if mappings(T,s['rows'],True):
                isomorphism_pairs+=1;need(i==j,'160 ordinary representatives pairwise distinct')
            p=mappings(reverse(T),s['rows'],True)
            if p:candidate.append((j,p[0]))
        need(len(candidate)==1,'each converse has exactly one representative partner')
        j,p=candidate[0];partners.append(j);certificates[i]['converse_partner_index']=j;certificates[i]['converse_to_partner_map']=p
    need(all(partners[partners[i]]==i for i in range(160)),'converse partner involution')
    need(all((partners[i]==i)==certificates[i]['self_converse'] for i in range(160)),'partner fixed points equal direct self-converse tests')
    merged=Counter(records[i]['expected_error'] for i,j in enumerate(partners) if i<=j)
    need(all(2*merged[e]==hist[e]+sc[e] for e in hist),'coefficientwise converse Burnside identity')
    # The regular-block matrix formula includes zero-error classes; check the
    # complete five-cell slice separately so the subtraction is visible.
    five_counts=defaultdict(Counter)
    for rec in inherited['representatives']:
        if rec['sector']!='five':continue
        T=tuple(rec['rows']);e,group,_=error(T);is_sc=bool(mappings(T,reverse(T),True))
        five_counts[rec['r']][(e,is_sc)]+=1
    # Sharpness and insufficiency controls: odd central bad core, doubled bad
    # core at n=14, and an actual size-paired non-self-converse ten-tournament.
    controls=[]
    for label,T,expected_e,expected_sc in [('central_bad7_order9',ordinal((0,),BAD7,(0,)),14,True),('paired_bad7_order14',ordinal(BAD7,reverse(BAD7)),28,True)]:
        e,group,_=error(T);anti=mappings(T,reverse(T));need(e==expected_e and bool(anti)==expected_sc,'SCC error refinement boundary control')
        controls.append({'label':label,'rows':T,'error':e,'stable_sizes':list(map(len,cells(T))),'self_converse':bool(anti),'certificate':anti_certificate(T,anti[0])})
    counter=next(r for r in certificates if r['source']=='strong_five' and not r['self_converse'])
    need(sorted(map(len,counter['stable_cells']))==[5,5],'paired stable size condition not sufficient')
    controls.append({'label':'paired_stable_sizes_not_sufficient','representative_index':counter['index'],'rows':counter['rows'],'error':counter['error'],'stable_sizes':[5,5],'self_converse':False})
    # All small labelled tournaments independently test the all-order lemma.
    small_tests=Counter()
    for n in range(1,6):
        pairs=list(combinations(range(n),2))
        for bits in range(1<<len(pairs)):
            T=[0]*n
            for k,(i,j) in enumerate(pairs):
                a,b=(i,j) if bits>>k&1 else (j,i);T[a]|=1<<b
            T=tuple(T)
            for anti in mappings(T,reverse(T)):
                anti_certificate(T,anti);small_tests[n]+=1
    need(hist==Counter({14:12,18:46,36:4,40:84,72:14}),'ordinary ten spectrum reproduced')
    result={'status':'PASS','scope':'complete positive-error order-ten slice, not all order-ten tournaments or metagraph edges','ordinary_positive_histogram':sorted(hist.items()),'strong_positive_histogram':sorted(strong_hist.items()),'self_converse_positive_histogram':sorted(sc.items()),'merged_positive_histogram':sorted(merged.items()),'ordinary_positive_classes':sum(hist.values()),'self_converse_positive_classes':sum(sc.values()),'merged_positive_nodes':sum(merged.values()),'ordinary_error_sum':sum(e*c for e,c in hist.items()),'self_converse_error_sum':sum(e*c for e,c in sc.items()),'merged_error_sum':sum(e*c for e,c in merged.items()),'merged_raw_second_moment':sum(e*e*c for e,c in merged.items()),'self_converse_by_sector':[[list(k),v] for k,v in sorted(sc_by_sector.items())],'five_cell_complete_checks':{str(d):[[list(k),v] for k,v in sorted(histogram.items())] for d,histogram in five_counts.items()},'small_labelled_anti_maps_checked':dict(small_tests),'controls':controls,'representatives':certificates,'checks':dict(CHECKS)}
    Path(args.output).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('representatives','checks','controls')},indent=2,sort_keys=True))
if __name__=='__main__':main()
