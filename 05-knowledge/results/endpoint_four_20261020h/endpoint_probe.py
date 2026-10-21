"""Owned continuation of the frozen endpoint-three bank; literal sources.

No isomorphism quotient. Search outputs are certificates, not theorem oracles.
"""
from pathlib import Path
import sys
import json
from itertools import combinations
from collections import Counter
from hashlib import sha256
from dataclasses import dataclass
import random

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / 'endpoint_three_repairs_20261020h'
sys.path.insert(0, str(OLD))
import probe as P
import root_check as RC
import probe_complete as PC


def sources(n, k):
    for middle in P.OLD.blue_presentations(n-2):
        for incoming in combinations(range(2,n-1),k-1):
            outgoing = tuple(sorted(n-1-v for v in incoming))
            rows = [sum(1<<v for v in range(1,n-1))]
            rows += [(row<<1)|(1<<(n-1)) for row in middle]
            rows += [1]
            for u in incoming:
                rows[0] ^= 1<<u
                rows[u] ^= 1
            for v in outgoing:
                rows[v] ^= 1<<(n-1)
                rows[-1] ^= 1<<v
            yield tuple(rows), incoming


def arcs(t):
    return {(u,v) for u in range(len(t)) for v in range(len(t)) if t[u]>>v&1}


def success(t, words):
    start=t
    for q in words:
        if not P.path(t,q) or P.blue(t,q):
            return False
        t=P.complement(t,q)
    if t != P.complement(start,tuple(range(len(t)))):
        raise RuntimeError('wrong literal target')
    return True


def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate_source(t):
    require(type(t) is tuple and len(t) >= 2, 'exact source tuple required')
    n=len(t)
    require(all(type(row) is int and 0 <= row < (1<<n) for row in t),
            'exact integer bit rows required')
    require(all(not (t[u]>>u&1) for u in range(n)), 'no source loops')
    require(all(((t[u]>>v&1)+(t[v]>>u&1)) == 1
                for u in range(n) for v in range(u)), 'complete orientation required')
    return n


def guarded_layer(t,k):
    n=validate_source(t)
    require(type(k) is int and k >= 1, 'exact endpoint layer required')
    p=tuple(range(n))
    require(P.path(t,p) and P.blue(t,p) and P.arc(t,n-1,0),
            'literal closed standard-blue source required')
    require(t[0].bit_count()==n-1-k and t[-1].bit_count()==k,
            'wrong marked endpoint degrees')
    return tuple(v for v in range(1,n-1) if P.arc(t,v,0))


@dataclass(frozen=True)
class Repair:
    source: tuple
    branch: str
    words: tuple
    target: tuple


def audit(receipt, expected_source):
    validate_source(expected_source)
    require(type(receipt) is Repair, 'typed repair required')
    validate_source(receipt.source)
    require(receipt.source==expected_source, 'literal source mismatch')
    require(type(receipt.branch) is str and bool(receipt.branch), 'branch mark required')
    n=len(expected_source)
    require(type(receipt.words) is tuple and len(receipt.words)==3, 'three words required')
    for q in receipt.words:
        require(type(q) is tuple and len(q)==n and all(type(v) is int for v in q),
                'exact integer path tuple required')
        require(set(q)==set(range(n)), 'labelled permutation required')
    validate_source(receipt.target)
    require(receipt.target==P.complement(expected_source,tuple(range(n))),
            'wrong literal target')
    require(success(expected_source,receipt.words), 'not an actual three-black repair')
    # Independent directed-set implementation; also returns literal anti-failures.
    PC.independent_replay(expected_source,receipt.words)
    return receipt


def package(t,name,words):
    return audit(Repair(t,name,tuple(tuple(q) for q in words),
                        P.complement(t,tuple(range(len(t))))),t)


TEN_SPECIAL = {
 (438,724,808,339,676,842,149,812,594,73): (
  (0,1,2,3,4,5,8,6,7,9),(0,1,2,3,5,8,6,7,9,4),(0,1,2,3,5,6,7,8,9,4)),
 (414,628,616,466,420,841,657,806,582,25): (
  (0,1,2,3,4,5,8,6,7,9),(0,1,2,3,9,8,7,4,6,5),(0,1,2,3,9,7,4,6,8,5)),
 (126,44,120,752,994,704,898,775,559,7): (
  (0,1,2,3,4,6,7,8,5,9),(0,1,2,3,8,5,9,4,6,7),(0,1,2,3,8,9,4,5,6,7)),
}


def repair_ten(t):
    incoming=guarded_layer(t,3)
    require(len(t)==10, 'order-ten consumer only')
    name,words=RC.ledger(10,incoming,arcs(t))
    if success(t,words):
        return package(t,'inherited/'+name,words)
    require(t in TEN_SPECIAL, 'source absent from the three exceptional rows')
    return package(t,'finite/order10-exception',TEN_SPECIAL[t])


def fourth_words(n,incoming):
    require(type(n) is int and n>=7, 'order at least seven required')
    require(type(incoming) is tuple and all(type(v) is int for v in incoming),
            'exact incoming port tuple required')
    require(tuple(sorted(set(incoming)))==incoming and len(incoming)==3,
            'three distinct sorted incoming ports required')
    a,b,x,y,c,z=0,1,2,n-3,n-2,n-1
    require(y in incoming and c in incoming and all(2<=v<=c for v in incoming),
            'open end-pair plus one incoming port required')
    h=next(v for v in incoming if v not in (y,c))
    k=n-4 if h==n-4 else y
    left=tuple(range(3,k)); right=tuple(range(k,y+1))
    return ((a,b,c,z,x)+left+right,
            (b,c,z,x)+left+(a,)+right,
            (b,x)+left+(a,)+right+(c,z))


def repair_fourth_open(t):
    incoming=guarded_layer(t,4)
    n=len(t)
    require(n>=7 and P.arc(t,1,n-2), 'open middle required')
    return package(t,'proved/open-endpair-adaptive-cut',fourth_words(n,incoming))


def repair_fourth_xy(t):
    incoming=guarded_layer(t,4)
    n=len(t); a,b,x,y,c,z=0,1,2,n-3,n-2,n-1
    require(n>=14 and y in incoming and c in incoming, 'end-pair scope n>=14 required')
    h=next(v for v in incoming if v not in (y,c))
    require(h not in (x,n-4) and P.arc(t,c,x) and P.arc(t,x,y),
            'extra-port and c->x->y guards required')
    w=tuple(range(3,n-3))
    words=((a,x)+w+(y,c,z,b),(z,)+w+(y,x,b,a,c),(b,z)+w+(y,x,a,c))
    return package(t,'proved/endpair-xy-extension',words)


def ordered_chart(t,q):
    n=validate_source(t)
    require(type(q) is tuple and len(q)==n and all(type(v) is int for v in q),
            'exact marked path required')
    require(set(q)==set(range(n)) and P.path(t,q), 'actual Hamiltonian chart required')
    return tuple((i,j-1) for i in range(n) for j in range(i+2,n)
                 if P.arc(t,q[j],q[i]))


def chart_trace(receipt):
    audit(receipt,receipt.source)
    t=receipt.source; trace=[]
    for q in receipt.words:
        before=ordered_chart(t,q)
        after=P.complement(t,q)
        chart_after=ordered_chart(after,q)
        universe={(i,j) for i in range(len(t)-1) for j in range(i+1,len(t)-1)}
        require(set(before).isdisjoint(chart_after) and set(before)|set(chart_after)==universe,
                'marked chart complement failed')
        trace.append({'path':q,'before_edges':before,'after_edges':chart_after})
        t=after
    return trace


def verify():
    checks=0
    def check(condition,message):
        nonlocal checks
        checks+=1
        require(condition,message)
    def rejects(f):
        nonlocal checks
        checks+=1
        try: f()
        except (ValueError,TypeError): return
        raise ValueError('malformed or off-scope receipt accepted')
    ten=Counter(); branches=Counter(); h10=sha256(); specials=[]
    for t,incoming in sources(10,3):
        r=repair_ten(t)
        check(guarded_layer(t,3)==incoming,'source enumeration ports changed')
        ten['sources']+=1; branches[r.branch]+=1
        if r.branch.startswith('finite/'):
            ten['finite_exceptions']+=1
            specials.append({'rows':t,'incoming':incoming,'words':r.words,
                             'target':r.target,'anti_failures':PC.independent_replay(t,r.words),
                             'ordered_chart_trace':chart_trace(r)})
        else: ten['inherited_selected']+=1
        h10.update(json.dumps([t,incoming,r.branch,r.words],separators=(',',':')).encode()+b'\n')
    check(dict(ten)=={'sources':86016,'inherited_selected':86013,'finite_exceptions':3},
          'complete order-ten universe changed')
    check({tuple(row['rows']) for row in specials}==set(TEN_SPECIAL),'exception list incomplete')
    fourth={}; h4=sha256(); samples=[]
    for n in range(7,11):
        count=0
        for t,incoming in sources(n,4):
            if not (n-3 in incoming and n-2 in incoming and P.arc(t,1,n-2)): continue
            r=repair_fourth_open(t)
            check(r.source==t,'fourth-layer source changed')
            count+=1
            h4.update(json.dumps([t,incoming,r.words],separators=(',',':')).encode()+b'\n')
            if count==1: samples.append({'rows':t,'incoming':incoming,'words':r.words})
        fourth[str(n)]=count
    check(fourth=={'7':16,'8':96,'9':1024,'10':10240},'fourth-layer census changed')
    larger=Counter(); hlarger=sha256(); rng=random.Random(2026102004)
    for n in tuple(range(14,31))+(40,60):
        for h in range(3,n-4):
            incoming=(h,n-3,n-2)
            for rep in range(8):
                aa,_=RC.make_source(n,incoming,rng)
                if (n-2,2) not in aa or (2,n-3) not in aa: continue
                t=tuple(sum(1<<v for u,v in aa if u==i) for i in range(n))
                r=repair_fourth_xy(t)
                check(r.source==t,'xy-extension source changed')
                larger[str(n)]+=1
                hlarger.update(json.dumps([t,incoming,r.words],separators=(',',':')).encode()+b'\n')
    # Exact permutation-orbit controls support the all-order proof, not its scope.
    for n in range(7,129):
        q=fourth_words(n,(n-5,n-3,n-2))[0]
        rq={q[i]:q[n-1-i] for i in range(n)}
        g=tuple(rq[n-1-i] for i in range(n))
        orbit=[0]
        while g[orbit[-1]]!=0: orbit.append(g[orbit[-1]])
        if n%2: expected=set(range(n))
        elif n%4==2: expected=set(range(0,n,2))
        else: expected={0,2}|set(range(6,n-5,4))
        check(set(orbit)==expected,'first-word orbit formula failed')
        check(P.edges(fourth_words(n,(n-5,n-3,n-2))[0]) ^
              P.edges(fourth_words(n,(n-5,n-3,n-2))[1]) ^
              P.edges(fourth_words(n,(n-5,n-3,n-2))[2]) == P.edges(tuple(range(n))),
              'adaptive XOR ledger failed')
    t=(6,36,8,19,103,77,15)
    adapted=repair_fourth_open(t)
    old=((0,1,5,6,2,3,4),(1,5,6,2,3,0,4),(1,2,3,0,4,5,6))
    check(not P.path(P.complement(t,old[0]),old[1]),'old cut unexpectedly legal')
    check(not P.arc(P.complement(t,old[0]),3,0),'lost join witness changed')
    chart_trace(adapted)
    first=next(iter(TEN_SPECIAL)); rt=repair_ten(first)
    rejects(lambda: repair_ten(tuple(float(v) for v in first)))
    bad=list(first); bad[0]=True
    rejects(lambda: repair_ten(tuple(bad)))
    rejects(lambda: repair_ten(list(first)))
    rejects(lambda: audit(rt,next(t for t in TEN_SPECIAL if t!=first)))
    rejects(lambda: audit(Repair(first,rt.branch,
                                 (tuple(float(v) for v in rt.words[0]),)+rt.words[1:],rt.target),first))
    rejects(lambda: audit(Repair(first,rt.branch,rt.words,first),first))
    rejects(lambda: fourth_words(10,(4,7,7)))
    rejects(lambda: fourth_words(10,(True,7,8)))
    rejects(lambda: fourth_words(6,(2,3,4)))
    rejects(lambda: repair_fourth_open(first))
    rejects(lambda: repair_fourth_xy(first))
    rejects(lambda: ordered_chart(t,tuple(reversed(range(7)))))
    return {'status':'FINITE-EXACT census and PROVED scoped fourth-layer interface',
            'checks':checks,
            'order10_endpoint3':{'universe':'all literal closed standard-blue sources; no filters',
                'counts':dict(ten),'branches':dict(sorted(branches.items())),
                'certificate_stream_sha256':h10.hexdigest(),'exceptions':specials},
            'endpoint4_open_endpair':{'all_order_scope':'n>=7, I={h,n-3,n-2}, middle 1->n-2',
                'finite_census':fourth,'certificate_stream_sha256':h4.hexdigest(),'samples':samples},
            'endpoint4_xy_endpair':{'all_order_scope':'n>=14, I={h,y,c}, h notin {2,n-4}, c->2->y; no b/c guard',
                'deterministic_guarded_samples':dict(larger),
                'certificate_stream_sha256':hlarger.hexdigest()},
            'old_cut_hostile':{'rows':t,'incoming':(3,4,5),'failed_words':old,
                'first_failed_join':(3,0),'adaptive_words':adapted.words},
            'scope_exclusions':['not every fourth-layer branch','not black-core connectivity',
                                'no unmarked graph or tournament-isomorphism quotient']}


if __name__=='__main__':
    print(json.dumps(verify(),indent=2,sort_keys=True))
