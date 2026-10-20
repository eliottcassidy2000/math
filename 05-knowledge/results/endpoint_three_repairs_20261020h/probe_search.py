"""Exact searches only for minimal inherited-bank residuals and one hostile.

The labelled three-K search is complete when exhaustive=True. For the
ten-vertex control it stops on the first certificate and reports that scope.
"""
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import importlib.util
import json
import probe as P

OLD_SEARCH = P.OLD_PATH.with_name('probe_search.py')
spec = importlib.util.spec_from_file_location('endpoint_two_search',OLD_SEARCH)
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)


def find_one(t):
    n = len(t)
    pairs = list(combinations(range(n),2))
    bits = {pair:1<<i for i,pair in enumerate(pairs)}
    original = sum(bits[(i,i+1)] for i in range(n-1))
    partner = P.complement(t,tuple(range(n)))
    first,last = S.path_masks(t,bits),S.path_masks(partner,bits)
    examined = 0
    for q,qmask in first:
        middle = P.complement(t,q)
        for s,smask in last:
            examined += 1
            r = S.recover_path(original^qmask^smask,pairs,n)
            if r is None:
                continue
            if not P.path(middle,r):
                r = r[::-1]
            if P.path(middle,r) and not P.blue(middle,r):
                words=(q,r,s)
                evaluation=P.evaluate(t,words)
                P.need(evaluation['success'],'first certificate replay')
                return {'scope':'prefix search stopped at first certificate',
                        'first_black_paths':len(first),'last_black_paths':len(last),
                        'pairs_examined':examined,'words':words,'evaluation':evaluation}
    return {'scope':'complete labelled search; no three-K certificate',
            'first_black_paths':len(first),'last_black_paths':len(last),
            'pairs_examined':examined,'words':None}


def cyclic_cut_hostile():
    n,m=10,5
    t=[0]*n
    for i,j in combinations(range(n),2):
        if i%2==j%2:
            forward=((j//2-i//2)%m) in {1,2}
        else:
            even,odd = (i,j) if i%2==0 else (j,i)
            even_to_odd = odd//2 != (even//2-1)%m
            forward = even_to_odd if i==even else not even_to_odd
        t[i if forward else j] |= 1<<(j if forward else i)
    return tuple(t)


def main():
    records=[]
    # All six minimal order-six residuals, not one per coarse profile.
    for t,u,v in P.endpoint_three(6):
        other=P.complement(t,tuple(range(6)))
        if not P.strong(other) or any(P.evaluate(t,w)['success'] for _,w in P.bank(6)):
            continue
        answer=S.search(t)
        best=answer['best'][0] if answer['best'] else None
        records.append({'n':6,'rows':t,'incoming':u,'outgoing':v,
                        'scope':'complete labelled three-K search',
                        'counts':{k:answer[k] for k in ['first_black_paths','last_black_paths','pairs_examined','three_K_certificates']},
                        'best':best})
        print('minimal',u,t,'certificates',answer['three_K_certificates'],flush=True)
    # A bounded sample: first residual in each profile at order seven.
    diagnosed=json.loads(Path(P.__file__).with_suffix('.json').read_text())
    for row in diagnosed['one_residual_per_profile']:
        if row['n']!=7:
            continue
        t=tuple(row['rows'])
        answer=find_one(t)
        records.append({**row,**answer})
        print('profile7',row['incoming'],row['profile'],'found',answer['words'] is not None,flush=True)
    t=cyclic_cut_hostile()
    p=tuple(range(10))
    incoming=tuple(v for v in range(1,9) if P.arc(t,v,0))
    outgoing=tuple(v for v in range(1,9) if P.arc(t,9,v))
    P.source_checks(t,incoming,outgoing)
    cut=p[1:]+p[:1]
    P.need(P.path(t,cut) and P.blue(t,cut),'genuine cyclic-cut blue hostile')
    answer=None
    for name,words in P.bank(10):
        result=P.evaluate(t,words)
        if result['success']:
            answer={'scope':'inherited bank certificate','name':name,'words':words,'evaluation':result}
            break
    if answer is None:
        answer=find_one(t)
    hostile={'n':10,'rows':t,'incoming':incoming,'outgoing':outgoing,
             'cyclic_cut':cut,'cyclic_cut_is_blue':True,**answer}
    print('hostile10',hostile['scope'],'found',hostile['words'] is not None,flush=True)
    receipt={'status':'FINITE-EXACT bounded residual search and explicit cyclic-cut hostile',
        'source_hashes':{'current_search':sha256(Path(__file__).read_bytes()).hexdigest(),
                         'current_probe':sha256(Path(P.__file__).read_bytes()).hexdigest(),
                         'inherited_endpoint_two_search':sha256(OLD_SEARCH.read_bytes()).hexdigest()},
        'records':records,'cyclic_cut_hostile':hostile}
    Path(__file__).with_suffix('.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    main()
