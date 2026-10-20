"""Compact complete n6..9 replay after adding the new arbitrary-block ledgers.

All source filters and all failures are explicit. This is a finite bank
coverage statement, not a claim that a remaining source is a bridge.
"""
from collections import Counter,defaultdict
from hashlib import sha256
from pathlib import Path
import json
import probe as P


def new_bank(n):
    a,b,x,y,c,z=0,1,2,n-3,n-2,n-1
    d=tuple(range(3,n-3))
    raw=[('last_c_x',((a,b,x)+d+(z,y,c),(b,x)+d+(z,y,a,c),(b,x)+d+(y,a,c,z)))]
    raw.append(('open_end_pair',((a,b,c,z,x)+d+(y,),
                                 (b,c,z,x)+d+(a,y),(b,x)+d+(a,y,c,z))))
    for cut in range(3,n-2):
        head,tail=tuple(range(1,cut)),tuple(range(cut,n-2))
        raw.append(('closed_last_c_cut_'+str(cut),
                    ((a,)+head+(z,)+tail+(c,),head+(z,)+tail+(a,c),
                     head+tail+(a,c,z))))
    if n>=7:
        e=n-4
        interior=tuple(range(3,n-4))
        raw.append(('C_two_ports',((a,b,x)+interior+(z,e,y,c),
            (a,x)+interior+(z,e,y,b,c),(a,x)+interior+(e,y,b,c,z))))
        raw.append(('closed_end_pair_final',((a,x)+d+(y,c,z,b),
            (z,)+d+(y,x,b,a,c),(b,z)+d+(y,x,a,c))))
    if n>=8:
        w=tuple(range(3,n-3))
        h=(x,)+w[::-1]+(c,y,z)
        raw.append(('mirrored_TT3_nonmodule',((a,x)+w+(y,c,z,b),
                                             ((a,b)+h)[::-1],(b,)+h[::-1]+(a,))))
    for name,words in raw:
        for suffix,qs in [('',words),('_reflected',P.OLD.reflect_words(n,words))]:
            compiled=P.compile_ledger(n,qs)
            P.need(not compiled['guard_conflicts'] and compiled['path_XOR_is_standard'],'new ledger compiler')
            yield name+suffix,qs


def main():
    summaries,residuals=[],[]
    digest=sha256()
    for n in range(6,10):
        stats,marks=Counter(),defaultdict(Counter)
        candidates=list(P.bank(n))+list(new_bank(n))
        first={}
        for t,u,v in P.endpoint_three(n):
            P.source_checks(t,u,v)
            stats['universe']+=1
            other=P.complement(t,tuple(range(n)))
            if not P.strong(other):
                stats['excluded_nonstrong_partner']+=1
                continue
            P.need(any(not P.blue(t,q) for q in P.paths(t)) and any(not P.blue(other,q) for q in P.paths(other)),
                   'retained endpoints actually mixed')
            stats['strong_mixed']+=1
            kind='open' if P.arc(t,1,n-2) else 'closed'
            key=kind+' U='+','.join(map(str,u))
            marks[key]['total']+=1
            for name,words in candidates:
                if P.evaluate(t,words)['success']:
                    stats['covered']+=1
                    marks[key]['covered']+=1
                    digest.update(json.dumps([t,u,name,words],separators=(',',':')).encode()+b'\n')
                    break
            else:
                stats['residual']+=1
                marks[key]['residual']+=1
                first.setdefault(key,{'n':n,'rows':t,'partner_rows':other,'incoming':u,'outgoing':v,'kind':kind})
        P.need(stats['universe']=={6:12,7:96,8:640,9:7680}[n],'expected full universe')
        summaries.append({'n':n,**dict(stats),'port_patterns':{k:dict(v) for k,v in sorted(marks.items())}})
        residuals+=list(first.values())
        print(json.dumps({'n':n,**dict(stats)}),flush=True)
        print('residual patterns',{k:v['residual'] for k,v in marks.items() if v.get('residual')},flush=True)
    receipt={'status':'FINITE-EXACT coverage by the explicitly tested refined local bank',
             'universe':'All closed standard-blue (n-4,3) sources n6..9; strong mixed endpoint filter explicit',
             'source_hashes':{'refined':sha256(Path(__file__).read_bytes()).hexdigest(),
                              'probe':sha256(Path(P.__file__).read_bytes()).hexdigest()},
             'certificate_stream_sha256':digest.hexdigest(),
             'summaries':summaries,'one_residual_per_port_pattern':residuals}
    Path(__file__).with_suffix('.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    main()
