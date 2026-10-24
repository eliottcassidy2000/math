#!/usr/bin/env python3
"""Independent Session51 converse audit; no producer implementation imports."""
from collections import Counter
from itertools import combinations, permutations, product
from pathlib import Path
import hashlib, json, math
HERE = Path(__file__).resolve().parent
CHECKS = Counter()
def need(ok, what):
    CHECKS[what] += 1
    if not ok: raise RuntimeError(what)
def arc(a,u,v): return (a[u] >> v) & 1
def validate(a):
    n=len(a)
    need(all(0 <= r < (1<<n) and not arc(a,i,i) for i,r in enumerate(a)), 'legal tournament rows')
    need(all(arc(a,i,j)+arc(a,j,i)==1 for i,j in combinations(range(n),2)), 'exact tournament orientation')
def reverse(a):
    n=len(a); full=(1<<n)-1
    return tuple(full ^ r ^ (1<<i) for i,r in enumerate(a))
def induced(a,vertices):
    return tuple(sum(arc(a,u,v)<<j for j,v in enumerate(vertices)) for u in vertices)
def stable(a):
    cells=[tuple(range(len(a)))] if a else []
    while True:
        new=[]
        for cell in cells:
            groups={}
            for v in cell:
                sig=tuple(sum(arc(a,v,w) for w in c) for c in cells)
                groups.setdefault(sig,[]).append(v)
            new.extend(tuple(groups[s]) for s in sorted(groups))
        if new==cells:return tuple(cells)
        cells=new
def signatures(a):
    # A recursively refined signature is invariant under ordinary isomorphism.
    labels=[0]*len(a)
    while True:
        keys=[(labels[i],tuple(sorted(labels[j] for j in range(len(a)) if arc(a,i,j)))) for i in range(len(a))]
        tab={k:i for i,k in enumerate(sorted(set(keys)))}
        nxt=[tab[k] for k in keys]
        if nxt==labels:return tuple(labels)
        labels=nxt
def maps(a,b,first=False):
    if len(a)!=len(b):return []
    n=len(a); da=tuple(x.bit_count() for x in a); db=tuple(x.bit_count() for x in b)
    if sorted(da)!=sorted(db):return []
    la,lb=signatures(a),signatures(b)
    if sorted(la)!=sorted(lb):return []
    cand=[{j for j in range(n) if (da[i],la[i])==(db[j],lb[j])} for i in range(n)]
    results=[]; chosen={}; used=set()
    def visit():
        if len(chosen)==n:
            results.append(tuple(chosen[i] for i in range(n)));return first
        remaining=[]
        for i in range(n):
            if i in chosen:continue
            feasible=[j for j in cand[i]-used if all(arc(a,i,k)==arc(b,j,v) for k,v in chosen.items())]
            if not feasible:return False
            remaining.append((len(feasible),i,feasible))
        _,i,feasible=min(remaining)
        for j in sorted(feasible):
            chosen[i]=j;used.add(j)
            if visit():return True
            used.remove(j);del chosen[i]
        return False
    visit();return results

def defect(a,aut=None):
    aut=maps(a,a) if aut is None else aut
    remaining=set(range(len(a)));orbits=[]
    while remaining:
        u=min(remaining);o={p[u] for p in aut};remaining-=o;orbits.append(o)
    cells=stable(a)
    return sum(next(len(c) for c in cells if min(o) in c)-len(o) for o in orbits)
def strong(a):
    if not a:return False
    for rows in (a,reverse(a)):
        reached=1;old=0
        while old!=reached:
            old=reached
            for i in range(len(a)):
                if (reached>>i)&1:reached|=rows[i]
        if reached!=(1<<len(a))-1:return False
    return True

def unique(rows):
    reps=[]
    for a in rows:
        if not any(maps(a,b,True) for b in reps):reps.append(a)
    return reps

def join(parts):
    rows=[];total=sum(map(len,parts));start=0
    for a in parts:
        following=((1<<total)-1)^((1<<(start+len(a)))-1)
        rows.extend((r<<start)|following for r in a);start+=len(a)
    return tuple(rows)
def substitute_cycle(core):
    n=len(core);rows=[r|(1<<n) for r in core]+[1<<(n+1),(1<<n)-1]
    return tuple(rows)
def power(p,k):
    q=list(range(len(p)))
    for _ in range(k):q=[p[x] for x in q]
    return tuple(q)
def anti_control(a,antis):
    cells=stable(a);n=len(a)
    need(all(len(c)%2 for c in cells),'actual stable cells odd')
    for p in antis:
        order=1;q=p
        while q!=tuple(range(n)):
            q=tuple(p[x] for x in q);order+=1
        need((n<=1 and order==1) or order%4==2,'anti automorphism order twice odd except trivial orders')
        tau=p if n<=1 else power(p,order//2)
        need(power(tau,2)==tuple(range(n)),'odd power anti involution')
        need(all(arc(a,i,j)==arc(a,tau[j],tau[i]) for i,j in combinations(range(n),2)), 'extracted map remains anti')
        fixed=sum(tau[i]==i for i in range(n))
        need(fixed==n%2,'fixed vertex parity and at most one')
        images=[frozenset(tau[v] for v in c) for c in cells]
        need(all(image in [frozenset(c) for c in cells] for image in images),'anti involution maps actual stable cells')
        fixedcells=sum(image==frozenset(c) for image,c in zip(images,cells))
        need(fixedcells==n%2,'fixed actual cell parity and at most one')

def small_classes(maxn=6):
    levels=[[()]]
    for n in range(maxn):
        proposals=[]
        for a in levels[-1]:
            for mask in range(1<<n):
                b=tuple(r|(((mask>>i)&1)<<n) for i,r in enumerate(a))+(((1<<n)-1)^mask,)
                proposals.append(b)
        levels.append(unique(proposals))
    results=[]
    expected_a=[1,1,1,2,4,12,56]
    expected_f=[1,1,1,2,2,8,12]
    for n,reps in enumerate(levels):
        f=0;s=0;fs=0
        for a in reps:
            validate(a);anti=maps(a,reverse(a));f+=bool(anti);s+=strong(a);fs+=strong(a) and bool(anti)
            if a:anti_control(a,anti)
        need(len(reps)==expected_a[n],'small ordinary class count')
        need(f==expected_f[n],'small self converse count')
        results.append({'n':n,'ordinary':len(reps),'self_converse':f,'strong':s,'strong_self_converse':fs})
    # Verify SCC series identities coefficientwise, independently of generation.
    for n in range(maxn+1):
        expected= (results[n]['strong_self_converse'] if n else 1)
        expected+=sum(results[k]['strong']*results[n-2*k]['self_converse'] for k in range(1,n//2+1))
        need(expected==results[n]['self_converse'],'small SCC self converse series')
    return levels,results

def matrices(d):
    options=[sum(1<<j for j in c) for c in combinations(range(5),d)]
    cols=[0]*5;rows=[]
    def rec():
        k=len(rows)
        if k==5:
            yield tuple(rows);return
        for r in options:
            nxt=[cols[j]+((r>>j)&1) for j in range(5)]
            if any(c>d or c+4-k<d for c in nxt):continue
            old=cols[:];cols[:]=nxt;rows.append(r)
            yield from rec()
            rows.pop();cols[:]=old
    return list(rec())
def move(m,a,b):
    return tuple(sum(((m[(i-a)%5]>>((j-b)%5))&1)<<j for j in range(5)) for i in range(5))
def anti_matrix(m):
    return tuple(sum(((m[-j%5]>>(-i%5))&1)<<j for j in range(5)) for i in range(5))
def build5(m):
    rows=[]
    for i in range(5):rows.append((1<<((i+1)%5))|(1<<((i+2)%5))|(m[i]<<5))
    for j in range(5):rows.append((1<<(5+(j+1)%5))|(1<<(5+(j+2)%5))|sum((1-((m[i]>>j)&1))<<i for i in range(5)))
    return tuple(rows)
def matrix_checks():
    result=[];rigid_reps=[]
    for d in (1,2):
        universe=matrices(d);pending=set(universe);summary=Counter();fixed=Counter();sym=sum(all(((m[i]>>j)&1)==((m[j]>>i)&1) for i,j in combinations(range(5),2)) for m in universe)
        while pending:
            m=min(pending);orb={move(m,a,b) for a,b in product(range(5),repeat=2)}
            need(orb<=pending,'matrix translation orbits partition');pending-=orb
            a=build5(m);aut=maps(a,a);antis=maps(a,reverse(a));e=defect(a,aut)
            need(len(aut)*len(orb)==25,'matrix full vertex automorphism orbit stabilizer')
            need(bool(antis)==(anti_matrix(m) in orb),'matrix converse agrees full isomorphism')
            summary[(e,bool(antis))]+=1
            anti_control(a,antis)
            if e:rigid_reps.append(a)
        for a,b in product(range(5),repeat=2):
            count=sum(move(anti_matrix(m),a,b)==m for m in universe);fixed[count]+=1
            need(count==(sym if a==b else math.comb(2,d//2)),'affine anti fixed set formula')
        need(sym==(26 if d==1 else 112),'symmetric regular matrix coefficient')
        need(summary[(40,True)]==(4 if d==1 else 20),'rigid self converse split')
        need(sum(summary[(e,True)] for e in (0,40))==(6 if d==1 else 24),'all five sector self converse split')
        result.append({'degree':d,'labelled':len(universe),'symmetric':sym,'classes':[{'defect':e,'self_converse':sc,'count':c} for (e,sc),c in sorted(summary.items())], 'anti_fixed_histogram':sorted(fixed.items())})
    return result,rigid_reps

def ten_checks(levels,rigid_reps):
    inherited=json.loads((HERE/'sources/order-ten-audit-receipt.json').read_text())
    strong10=[tuple(x['rows']) for x in inherited['representatives'] if x['error']]
    need(len(strong10)==120,'inherited strong ten positive representative support')
    need(all(any(maps(a,b,True) for b in strong10) for a in rigid_reps),'independently generated five rigid classes match inherited atlas')
    seven_example=next(tuple(x['rows']) for x in inherited['representatives'] if x['sector']=='seven')
    r7=induced(seven_example,max(stable(seven_example),key=len));need(defect(r7)==14,'extracted unique defective regular seven')
    core9=[]
    for x in inherited['representatives']:
        if x['sector']!='triangles' or not x['error']:continue
        a=tuple(x['rows']);singleton=next(c[0] for c in stable(a) if len(c)==1)
        b=induced(a,[i for i in range(10) if i!=singleton])
        if len(stable(b))==3:core9.append(b)
    rigid9=unique(core9);need(len(rigid9)==4,'rigid three triangle strong nine classes')
    r9=json.loads((HERE/'sources/regular-nine-receipt.json').read_text());regular9=[]
    for comp in r9['components']:
        if not comp['error']:continue
        p=tuple(r9['states'][comp['states'][0]]['parent'])
        a=tuple(r|((r.bit_count()==3)<<8) for r in p)+(sum((r.bit_count()==4)<<i for i,r in enumerate(p)),)
        need(defect(a)==comp['error'],'full automorphism regular nine defect agrees inherited receipt')
        regular9.append(a)
    strong9=rigid9+regular9+[substitute_cycle(r7)]
    need(len(strong9)==17 and all(strong(a) for a in strong9),'positive strong nine representatives')
    full=list(strong10)
    for a in strong9:full.extend([join([(0,),a]),join([a,(0,)])])
    for i in range(4):
        for a in levels[i]:
            for b in levels[3-i]:full.append(join([a,r7,b]))
    need(len(full)==160,'full ten positive representative construction')
    hist=Counter();selfhist=Counter();proof=[]
    for i,a in enumerate(full):
        validate(a);aut=maps(a,a);e=defect(a,aut);hist[e]+=1;antis=maps(a,reverse(a));selfhist[e]+=bool(antis)
        need(e>0,'all constructed ten representatives positive error')
        need(all(not maps(a,b,True) for b in full[:i]),'full ten positive representatives ordinary pairwise nonisomorphic')
        mates=[j for j,b in enumerate(full) if maps(reverse(a),b,True)]
        need(len(mates)==1,'full ten positive universe uniquely converse closed')
        need(defect(full[mates[0]])==e,'converse invariant defect')
        anti_control(a,antis)
        proof.append({'rows':a,'defect':e,'automorphism_order':len(aut),'anti_automorphism_order':len(antis),'converse_partner':mates[0]})
    need(hist==Counter({14:12,18:46,36:4,40:84,72:14}),'ordinary full ten spectrum')
    need(+selfhist==Counter({40:24}),'self converse full ten spectrum')
    merged={e:(c+selfhist[e])//2 for e,c in hist.items()}
    need(merged=={14:6,18:23,36:2,40:54,72:7},'merged full ten spectrum')
    need(sum(merged.values())==92,'merged positive support')
    need(sum(e*c for e,c in merged.items())==3234,'merged positive first moment')
    need(sum(e*c for e,c in hist.items())==5508,'ordinary positive first moment')
    return {'ordinary_spectrum':sorted(hist.items()),'self_converse_spectrum':sorted((+selfhist).items()),'merged_spectrum':sorted(merged.items()),'representatives':proof},r7

def weighted_word_controls(r7,levels):
    # Converse-closed strong alphabet includes both self-converse and paired types.
    nonself=next(a for a in levels[6] if strong(a) and not maps(a,reverse(a),True))
    pool=[(0,),(2,4,1),r7,nonself,reverse(nonself)]
    dual=[0,1,2,4,3]
    weights=[(len(a),defect(a)) for a in pool]
    need(all(maps(reverse(a),pool[dual[i]],True) for i,a in enumerate(pool)),'selected strong SCC alphabet converse closed')
    # Full words through 16, so two defective components and the squared t substitution are tested.
    cap=16;words=[((),0,0)];cursor=0
    while cursor<len(words):
        word,n,e=words[cursor];cursor+=1
        for k,(size,defect_) in enumerate(weights):
            if n+size<=cap:words.append((word+(k,),n+size,e+defect_))
    ordinary=Counter();fixed=Counter()
    for word,n,e in words:
        ordinary[(n,e)]+=1
        if word==tuple(dual[i] for i in word[::-1]):fixed[(n,e)]+=1
    for n in range(cap+1):
        for e in (0,14,28):
            expected=(n,e)==(0,0)
            expected+=sum((n,e)==w for i,w in enumerate(weights) if dual[i]==i)
            expected+=sum(fixed[(n-2*s,e-2*d)] for s,d in weights if n>=2*s and e>=2*d)
            need(fixed[(n,e)]==expected,'weighted SCC palindrome recurrence with squared t')
    need(fixed[(14,28)]==1,'paired defective SCC positive control')
    need(fixed[(10,14)]==0,'even ten single defective SCC hostile')
    need(maps(join([r7,r7]),reverse(join([r7,r7])),True),'actual fourteen vertex paired defect SC')
    wrong_without_t_square=sum(fixed[(14-2*s,28-d)] for s,d in weights if 14>=2*s and 28>=d)
    need(wrong_without_t_square!=fixed[(14,28)],'missing t square mutation rejected')
    return {'selected_strong_weights':weights,'maximum_order':cap,'direct_words':len(words),'self_converse_words':sum(fixed.values()),'two_defective_component_coefficient':fixed[(14,28)]}

def source_pins():
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
    for path,sha in pins.items():
        need(hashlib.sha256((HERE/path).read_bytes()).hexdigest()==sha,'inherited source pin '+path)

def run():
    source_pins()
    levels,small=small_classes();matrices_,rigid=matrix_checks();ten,r7=ten_checks(levels,rigid)
    weighted=weighted_word_controls(r7,levels)
    return {'status':'PASS','scope':'Independent full isomorphism tests of all 160 classified positive-error order-ten classes; no unrestricted order-ten enumeration','small_orders':small,'five_matrix_sectors':matrices_,'order_ten':ten,'weighted_SCC_controls':weighted,'checks':dict(sorted(CHECKS.items()))}
if __name__=='__main__':print(json.dumps(run(),sort_keys=True,indent=2))
