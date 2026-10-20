"""Finite controls only: finish the eleven refined-bank residuals at n6..9.

There is no order-ten census and no all-order inference from these controls.
Every selected certificate gets an independent directed-arc-set replay.
"""
from collections import Counter
from hashlib import sha256
from pathlib import Path
import json
import probe as P
import probe_refined as R
import probe_search as S


def independent_replay(rows,words):
    n=len(rows)
    initial={(u,v) for u in range(n) for v in range(n) if rows[u]>>v&1}
    now=initial
    witnesses=[]
    for q in words:
        P.need(len(q)==n and set(q)==set(range(n)),'word permutation')
        retained=set(zip(q,q[1:]))
        P.need(retained<=now,'independent directedness')
        reversal=dict(zip(q,q[::-1]))
        failures=[(u,v,reversal[v],reversal[u]) for u,v in sorted(now)
                  if (reversal[v],reversal[u]) not in now]
        P.need(failures,'independent actual blackness')
        witnesses.append(failures[0])
        now=({(v,u) for u,v in now}-{(v,u) for u,v in retained})|retained
    original={(i,i+1) for i in range(n-1)}
    target=({(v,u) for u,v in initial}-{(v,u) for u,v in original})|original
    P.need(now==target,'independent literal target')
    return witnesses


def main():
    summaries,certificates,unresolved,excluded=[],[],[],[]
    digest=sha256()
    for n in range(6,10):
        stats=Counter()
        candidates=list(P.bank(n))+list(R.new_bank(n))
        for t,u,v in P.endpoint_three(n):
            P.source_checks(t,u,v)
            stats['universe']+=1
            other=P.complement(t,tuple(range(n)))
            if not P.strong(other):
                stats['excluded_nonstrong_partner']+=1
                P.need(n==6 and not any(not P.blue(other,q) for q in P.paths(other)),
                       'only excluded source has pure-blue n6 partner')
                excluded.append({'n':n,'rows':t,'partner_rows':other,'incoming':u,'outgoing':v,
                                 'reason':'pure-blue nonstrong partner C3->C3'})
                continue
            stats['retained']+=1
            for name,qs in candidates:
                if P.evaluate(t,qs)['success']:
                    stats['bank_covered']+=1
                    break
            else:
                result=S.find_one(t)
                if result['words'] is None:
                    unresolved.append({'n':n,'rows':t,'incoming':u,'outgoing':v,**result})
                    continue
                name,qs='finite_residual',tuple(tuple(q) for q in result['words'])
                stats['finite_certified']+=1
                certificates.append({'n':n,'rows':t,'partner_rows':other,'incoming':u,'outgoing':v,
                    'words':qs,'independent_anti_failures':independent_replay(t,qs),
                    'endpoint_excesses':[s['endpoint_excess'] for s in result['evaluation']['steps']],
                    'search_scope':result['scope'],'pairs_examined':result['pairs_examined']})
                print('finite certificate',n,u,'pairs',result['pairs_examined'],flush=True)
            independent_replay(t,qs)
            digest.update(json.dumps([t,u,name,qs],separators=(',',':')).encode()+b'\n')
        P.need(stats['universe']=={6:12,7:96,8:640,9:7680}[n],'complete finite universe')
        P.need(stats['retained']=={6:11,7:96,8:640,9:7680}[n],'explicit filter count')
        summaries.append({'n':n,**dict(stats)})
        print(json.dumps(summaries[-1]),flush=True)
    P.need(not unresolved,'finite residual lacks a three-black certificate')
    hostile=S.cyclic_cut_hostile()
    cut=tuple(range(1,10))+(0,)
    P.need(P.path(hostile,cut) and P.blue(hostile,cut),'sharp cyclic-cut control')
    ten=S.find_one(hostile)
    P.need(ten['words'] is not None,'ten-vertex hostile has a repair')
    ten_control={'rows':hostile,'cyclic_cut':cut,'cyclic_cut_blue':True,
                 'words':ten['words'],'independent_anti_failures':independent_replay(hostile,ten['words']),
                 'search_scope':ten['scope'],'pairs_examined':ten['pairs_examined']}
    report={'status':'FINITE-EXACT: every retained marked source at n6..9 has a literal three-black K detour',
            'scope':'Strong mixed endpoint filter; one n6 pure nonstrong partner excluded. No order-ten census.',
            'source_hashes':{**{p.name:sha256(p.read_bytes()).hexdigest() for p in
                             [Path(__file__),Path(P.__file__),Path(R.__file__),Path(S.__file__)]},
                             'inherited_endpoint_two_kernel':sha256(P.OLD_PATH.read_bytes()).hexdigest(),
                             'inherited_endpoint_two_search':sha256(S.OLD_SEARCH.read_bytes()).hexdigest()},
            'certificate_stream_sha256':digest.hexdigest(),'summaries':summaries,
            'finite_residual_certificates':certificates,'excluded':excluded,'unresolved':unresolved,
            'one_order_ten_hostile_control':ten_control}
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('PASS',sum(s['retained'] for s in summaries),'retained sources',len(certificates),'finite exceptions',flush=True)


if __name__=='__main__':
    main()
