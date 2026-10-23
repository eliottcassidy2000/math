"""Literal controls for full-orbit incidence obstructions; no center census.

--discover writes a small explicit control bank.  Default/--write/--verify
only consume its words and do not run Hamilton-path discovery.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import itertools as it
import json
from math import gcd
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
GATES=0
PARAMETERS=((3,5),(3,7),(4,5),(4,7),(6,3),(8,3),(2,9))
PINS=(
 '05-knowledge/results/odd_shore_deficit_centers_20261021/session22-centers/PROOF.md',
 '05-knowledge/results/center_transitivity_20261021c/README.md',
 '05-knowledge/results/center_transitivity_20261021c/VT15.md',
 '05-knowledge/results/center_transitivity_20261021c/receipt15.json',
 '05-knowledge/results/deficit_six_centers_20261021b/COMBINED_CUTOFF.md',
 '05-knowledge/results/deficit_twelve_20261021c/receipt.json',
 '05-knowledge/results/deficit_fourteen_20261021c/receipt.json',
 '05-knowledge/results/center_transitivity_20261021c/receipt35.json')


def need(value,why):
    global GATES
    GATES+=1
    if not value: raise ValueError(why)


def integer(x,lo=0,hi=None):
    if type(x) is not int or x<lo or hi is not None and x>hi:
        raise ValueError('exact integer and range')
    return x


def lf(p): return p.read_bytes().replace(b'\r\n',b'\n')
def sha(b): return hashlib.sha256(b).hexdigest()
def serial(x): return (json.dumps(x,sort_keys=True,indent=2)+'\n').encode()
def arc(t,a,b): return (t[a]>>b)&1


def source(t):
    if type(t) is not tuple or not t: raise ValueError('immutable source')
    n=len(t)
    for i,row in enumerate(t):
        integer(row,0,(1<<n)-1)
        if row>>i&1: raise ValueError('source loop')
        for j in range(i):
            if arc(t,i,j)+arc(t,j,i)!=1: raise ValueError('not tournament')
    return t


def permutation(p,vertices):
    if type(p) is not tuple or any(type(x) is not int for x in p) or sorted(p)!=sorted(vertices):
        raise ValueError('literal vertex permutation')
    return p


def path(t,p,vertices=None):
    permutation(p,range(len(t)) if vertices is None else vertices)
    if any(not arc(t,a,b) for a,b in zip(p,p[1:])): raise ValueError('non-native path')
    return p


def is_anti(t,p):
    permutation(p,range(len(t)))
    return all(arc(t,p[a],p[b])==arc(t,b,a) for a,b in it.combinations(range(len(t)),2))


def key(t,p):
    path(t,p)
    return sum(arc(t,p[a],p[b])<<k for k,(a,b) in enumerate(it.combinations(range(len(t)),2)))


def blue(t,p):
    path(t,p)
    rho=[0]*len(t)
    for a,b in zip(p,p[::-1]): rho[a]=b
    return is_anti(t,tuple(rho))


def target(t,p):
    path(t,p)
    kept={frozenset((a,b)) for a,b in zip(p,p[1:])}
    return sum((arc(t,a,b) if frozenset((a,b)) in kept else arc(t,b,a))<<k
               for k,(a,b) in enumerate(it.combinations(range(len(t)),2)))


def frame(q,m):
    integer(q,2);integer(m,3)
    if m%2!=1: raise ValueError('odd shore')
    n=q*m;t=[0]*n
    for a,b in it.combinations(range(n),2):
        s,i=a%q,a//q;r,j=b%q,b//q
        if s==r: forward=0<(j-i)%m<=(m-1)//2
        elif s<r: forward=(r==s+1 and i==j)
        else: forward=not(s==r+1 and i==j)
        x,y=(a,b) if forward else (b,a);t[x]|=1<<y
    return tuple(t)


def blocks(q,m):
    ans=[]
    for start in range(q):
        word=[];j=0
        for d in range(q):
            s=(start+d)%q
            word.extend(s+q*((j+k)%m) for k in range(m))
            j=(j+m-1+(s==q-1))%m
        ans.append(tuple(word))
    return tuple(ans)


def runs(w):
    ans=[]
    for a in w:
        if ans and ans[-1][0]==a: ans[-1]=(a,ans[-1][1]+1)
        else: ans.append((a,1))
    return tuple(ans)


def discover_path(t,vertices,start):
    """Producer only: capped discovery, never called by verify()."""
    left=set(vertices);left.remove(start);nodes=0
    def dfs(word):
        nonlocal nodes
        nodes+=1
        if nodes>1000000: raise ValueError('discovery cap, no claimed proof')
        if not left: return tuple(word)
        for v in sorted(left):
            if arc(t,word[-1],v):
                left.remove(v);ans=dfs(word+[v]);left.add(v)
                if ans is not None: return ans
        return None
    ans=dfs([start])
    if ans is None: raise ValueError('no discovery witness')
    return ans


def produce():
    records=[]
    for q,m in PARAMETERS:
        t=frame(q,m);n=q*m
        cut=[tuple((s+j)%n for j in range(n)) for s in range(q)]
        block=blocks(q,m);pref=None
        if m>=5:
            r=8 if q==2 else 4
            indices=tuple(range(r)) if m!=7 else (0,1,3,5)
            cyc=tuple(q*j+q-1 for j in indices)
            v=0;exits=tuple(x for x in cyc if arc(t,x,v))[:r//2+1]
            suffix=discover_path(t,set(range(n))-set(cyc),v)
            words=[]
            for e in exits:
                k=cyc.index(e);words.append(cyc[k+1:]+cyc[:k+1]+suffix)
            pref={'cycle':cyc,'terminal':v,'exits':exits,'suffix':suffix,'words':words}
        allwords=cut+list(block)+(pref['words'] if pref else [])
        records.append({'q':q,'m':m,'source':t,'anti':tuple(range(n-1,-1,-1)),
                        'cuts':cut,'blocks':block,'prefix':pref,
                        'targets':[target(t,p) for p in allwords]})
    (HERE/'inputs.json').write_bytes(serial({'schema':'incidence-input-v1','records':records,
        'pins':{p:sha(lf(ROOT/p)) for p in PINS}}))


def record_tuple(raw):
    if type(raw) is not dict or set(raw)!={'q','m','source','anti','cuts','blocks','prefix','targets'}:
        raise ValueError('record schema')
    r=dict(raw)
    for name in ('source','anti','targets'): r[name]=tuple(r[name])
    for name in ('cuts','blocks'): r[name]=tuple(tuple(p) for p in r[name])
    if r['prefix'] is not None:
        p=r['prefix']
        if type(p) is not dict or set(p)!={'cycle','terminal','exits','suffix','words'}: raise ValueError('prefix schema')
        r['prefix']={k:tuple(tuple(x) for x in v) if k=='words' else tuple(v) if k!='terminal' else v for k,v in p.items()}
    return r


def audit_record(raw):
    r=record_tuple(raw);q=integer(r['q'],2);m=integer(r['m'],3);t=source(r['source']);n=len(t)
    need((q,m) in PARAMETERS and n==q*m and t==frame(q,m),'literal source and allowed control')
    need(is_anti(t,r['anti']),'actual source anti-involution')
    need(all(arc(t,i,(i+1)%n) for i in range(n)),'actual marked Hamilton cycle')
    need(all(arc(t,a,b)==arc(t,(a+q)%n,(b+q)%n) for a,b in it.combinations(range(n),2)),'rotation automorphism')
    degrees=tuple(x.bit_count() for x in t)
    need(len(set(degrees))==q and all(degrees[i]==degrees[i%q] for i in range(n)),
         'distinct degree classes plus transitive rotation certify exact full orbits')
    need(r['cuts']==tuple(tuple((s+j)%n for j in range(n)) for s in range(q)),'all declared cut representatives')
    need(r['blocks']==blocks(q,m),'source-directed whole-shore paths')
    groups={};allwords=[]
    for name in ('cuts','blocks'):
        words=r[name]
        for p in words:path(t,p)
        need(len({key(t,p) for p in words})==q,'distinct actual full-Aut orbits')
        need(sum(blue(t,p) for p in words)<=gcd(2,q),'blue bound')
        groups[name]={key(t,p) for p in words if not blue(t,p)}
        allwords+=list(words)
    need(not groups['cuts']&groups['blocks'],'cut/block orbit disjointness')
    p=r['prefix'];prefix_count=0
    if p is not None:
        cyc=p['cycle'];rr=8 if q==2 else 4
        need(len(cyc)==rr and len(set(cyc))==rr and all(type(x) is int and 0<=x<n and x%q==q-1 for x in cyc),'literal proper prefix support')
        need(all(arc(t,a,b) for a,b in zip(cyc,cyc[1:]+cyc[:1])),'directed prefix cycle')
        v=integer(p['terminal'],0,n-1);path(t,p['suffix'],set(range(n))-set(cyc))
        need(p['suffix'][0]==v and len(p['exits'])==rr//2+1 and len(set(p['exits']))==len(p['exits']),'common suffix and exit count')
        expected=[]
        for e in p['exits']:
            need(type(e) is int and e in cyc and arc(t,e,v),'literal exit')
            k=cyc.index(e);expected.append(cyc[k+1:]+cyc[:k+1]+p['suffix'])
        need(p['words']==tuple(expected),'exact prefix/suffix splices')
        mates=tuple(tuple(r['anti'][x] for x in word[::-1]) for word in p['words'])
        for word in p['words']+mates:path(t,word)
        pg={key(t,w) for w in p['words']+mates if not blue(t,w)}
        need(len(pg)>=rr,'half-supported black-orbit and mate lower bound')
        if q>=3:
            need(not pg&(groups['cuts']|groups['blocks']),'all three actual black families disjoint')
            for word in p['words']+mates:
                run=runs(tuple(x%q for x in word))
                need(max(k for _,k in run)>=4 and len(run)>q,'split-orbit run sidecar')
        groups['prefix_and_mates']=pg;prefix_count=len(pg);allwords+=list(p['words'])
    need(len(r['targets'])==len(allwords),'target count')
    for word,claim in zip(allwords,r['targets']):
        integer(claim,0,(1<<(n*(n-1)//2))-1)
        need(claim==target(t,word),'source-bound complement-except-path target')
        # Independent target matrix: every actual path arc remains, all others flip.
        edges={frozenset((a,b)) for a,b in zip(word,word[1:])}
        for k,(a,b) in enumerate(it.combinations(range(n),2)):
            need((claim>>k&1)==(arc(t,a,b)^(frozenset((a,b)) not in edges)),'literal target edge')
    count=len(set().union(*groups.values()))
    need(count>=8,'eight actual black full-Aut path orbits')
    return {'q':q,'m':m,'n':n,'degrees':degrees[:q],'black_family_counts':{k:len(v) for k,v in groups.items()},'union_black_orbits':count,'source_sha256':sha(serial(t))}


def partitions(n):
    def go(w,maximum):
        if len(w)==n: yield tuple(w);return
        for x in range(maximum+2):yield from go(w+[x],max(maximum,x))
    yield from go([0],0)


def reflected_action(w):
    mapping={}
    for a,b in zip(w[::-1],w):
        if a in mapping and mapping[a]!=b:return None
        mapping[a]=b
    if len(set(mapping.values()))!=len(mapping):return None
    return tuple(sorted(mapping.items()))


def necklace_controls():
    totals=Counter();admissible={}
    for q in range(2,9):
        sizes=set()
        for w in partitions(q):
            c=Counter(w)
            if any(s%2==0 for s in c.values()):continue
            hist=Counter(c.values())
            # Odd full vertex orbits under an anti-involution: exactly n%2 fixed orbits.
            if sum(v%2 for v in hist.values())==q%2:sizes.add(tuple(sorted(c.values())))
            if 1 not in c.values():continue
            cut=[tuple(w[(s+i)%q] for i in range(3*q)) for s in range(q)]
            block=[tuple(a for j in range(q) for a in [w[(s+j)%q]]*3) for s in range(q)]
            need(len(set(cut))==q and len(set(block))==q,'singleton makes both necklaces faithful')
            need(not set(cut)&set(block),'singleton run separates both families')
            for family in (cut,block):
                cnt=Counter(a for word in family if (a:=reflected_action(word)) is not None)
                need(not cnt or max(cnt.values())<=gcd(2,q),'every possible anti label action has bounded blue starts')
            totals[q]+=1
        admissible[q]=sorted(sizes)
    need(admissible[4]==[(1,1,1,1)],'q4 anti-pairing orbit partitions')
    need(admissible[6]==[(1,1,1,1,1,1),(3,3)],'q6 anti-pairing orbit partitions')
    need(admissible[8]==[(1,1,1,1,1,1,1,1),(1,1,3,3)],'q8 anti-pairing orbit partitions')
    need(all(max(s)==7 or 1 in s for s in admissible[7]),'non-VT q7 has singleton')
    return {'singleton_necklaces':dict(totals),'paired_odd_orbit_partitions':admissible}


def refused(fn,*args):
    try:fn(*args)
    except (ValueError,TypeError):need(True,'hostile rejected')
    else:need(False,'hostile accepted')


def verify():
    global GATES
    GATES=0
    data=json.loads((HERE/'inputs.json').read_text())
    need(type(data) is dict and set(data)=={'schema','records','pins'} and data['schema']=='incidence-input-v1','input schema')
    need(data['pins']=={p:sha(lf(ROOT/p)) for p in PINS},'frozen exact mathematical inputs')
    need([(r['q'],r['m']) for r in data['records']]==list(PARAMETERS),'complete declared seven-control universe')
    controls=[audit_record(r) for r in data['records']]
    necklaces=necklace_controls()
    one=data['records'][0]
    def altered(k,v):return {**one,k:v}
    refused(audit_record,altered('q',3.0))
    refused(audit_record,altered('source',[float(x) for x in one['source']]))
    refused(audit_record,altered('anti',list(range(len(one['source'])))))
    wrong=[p[:] for p in one['blocks']];wrong[0][-1]=wrong[0][0]
    refused(audit_record,altered('blocks',wrong))
    wrong=one['targets'][:];wrong[0]^=1;refused(audit_record,altered('targets',wrong))
    wrong={**one['prefix'],'cycle':one['prefix']['cycle']*2}
    refused(audit_record,altered('prefix',wrong))
    wrong={**one['prefix'],'terminal':False};refused(audit_record,altered('prefix',wrong))
    # Exact small positive boundary: six-center has only six black orbits.
    t=(36,17,10,19,37,14);source(t);ks=set();bs=set()
    for p in it.permutations(range(6)):
        if all(arc(t,a,b) for a,b in zip(p,p[1:])):
            k=key(t,p);ks.add(k)
            if blue(t,p):bs.add(k)
    need(len(ks)-len(bs)==6,'known six-center positive boundary')
    need(2*(3-gcd(2,3))==4 and 2*(4-gcd(2,4))==4,'two families alone do not exclude q3 or q4')
    need(2*(6-gcd(2,6))==8 and 2*(8-gcd(2,8))==12,'higher singleton exclusion bounds')
    return {'schema':'center-incidence-20261022','checks':GATES,'status':'PROVED analytic exclusions with FINITE-EXACT literal controls, not a census',
      'excluded_orders':[18,20,24,28],'order21_necessary':'vertex-transitive',
      'nonVT_q3_m_ge5_excluded':True,'controls':controls,'necklaces':necklaces,
      'remaining_before_peer_order15_composition':[15,21],
      'scope':'No complete classification at order21; orbit-label equality is not an orbit-equality certificate.',
      'pins':data['pins'],'input_sha256':sha(lf(HERE/'inputs.json'))}


def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group()
    g.add_argument('--discover',action='store_true');g.add_argument('--write',action='store_true');g.add_argument('--verify',action='store_true')
    a=p.parse_args()
    if a.discover:produce();return
    payload=serial(verify())
    if a.write:(HERE/'receipt.json').write_bytes(payload)
    if a.verify and payload!=lf(HERE/'receipt.json'):raise ValueError('saved receipt mismatch')
    print(f'{GATES:,} exact checks; LF SHA256 {sha(payload)}')


if __name__=='__main__':main()
