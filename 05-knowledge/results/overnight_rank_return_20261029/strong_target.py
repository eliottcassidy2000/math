"""Actual strong-target complement probes, with sources and integer minors."""
from collections import Counter
import json
from time import perf_counter

from complement_kernel import HERE, complement, native, paths, rank_mod, unit_certificate, need
from complement_verify import strong, verify_minors


def build():
    source_bank=json.loads((HERE/'complement_units.json').read_text())
    result=[]
    start=perf_counter()
    for number,record in enumerate(source_bank,1):
        rows=record['source']
        entry=dict(n=len(rows),source=rows)
        if record['path'] is None:
            entry['zero_steps']=True
        elif strong(complement(rows,record['path'])):
            entry['certificate']=dict(path=record['path'],minors=record['minors'])
        else:
            M=native(rows)[2]
            tested=[]
            orbit_codes=set()
            for P in paths(rows):
                code=tuple(sum((rows[u]>>v&1)<<j for j,v in enumerate(P)) for u in P)
                if code in orbit_codes: continue
                orbit_codes.add(code)
                U=complement(rows,P)
                if not strong(U): continue
                N=native(U)[2]
                cert=unit_certificate(M+N)
                if cert is not None:
                    entry['certificate']=dict(path=P,minors=cert)
                    break
                obstruction=next((p for p in (2,3,5,7,11,13,17,19,23,29,31,37,1009)
                                  if rank_mod(M+N,p)<len(M)),None)
                tested.append(dict(path=P,prime=obstruction,
                                   rank=rank_mod(M+N,obstruction) if obstruction else None))
            if 'certificate' not in entry:
                entry['failed_strong_orbits']=tested
                print('NO STRONG SUCCESS',len(rows),rows,len(tested),flush=True)
        result.append(entry)
        if number%30==0: print(number,round(perf_counter()-start,1),flush=True)
    (HERE/'strong_target.json').write_text(json.dumps(result,separators=(',',':'))+'\n')
    print(Counter((x['n'],'zero' if x.get('zero_steps') else 'one' if 'certificate' in x else 'none') for x in result))


def verify():
    source_bank=json.loads((HERE/'complement_units.json').read_text())
    records=json.loads((HERE/'strong_target.json').read_text())
    need(len(records)==len(source_bank)==396, 'complete inherited and independently certified source bank')
    need(all(not strong(complement([2,4,1],P)) for P in paths([2,4,1])),
         'positive-length strong-target requirement already fails on C3; zero steps are allowed')
    counts=Counter()
    for number,(old,record) in enumerate(zip(source_bank,records),1):
        rows=record['source']
        need(rows==old['source'] and record['n']==old['n'] and strong(rows),'same actual strong source')
        M=native(rows)[2]
        if record.get('zero_steps'):
            need(old['path'] is None,'zero steps only with inherited unit certificate')
            verify_minors(M,old['minors'])
            counts[(len(rows),'zero')]+=1
        else:
            cert=record['certificate']
            U=complement(rows,cert['path'])
            need(strong(U) and complement(U,cert['path'])==rows,'actual strong target and involution')
            verify_minors(M+native(U)[2],cert['minors'])
            counts[(len(rows),'one')]+=1
        if number%75==0: print('strong verified',number,flush=True)
    need(counts==Counter({(3,'zero'):1,(4,'zero'):1,(5,'zero'):3,(5,'one'):3,
                         (6,'zero'):7,(6,'one'):28,(7,'zero'):16,(7,'one'):337}), 'strong-target census')
    print('PASS strong-target integral recovery on all396 sources')


if __name__=='__main__':
    import sys
    build() if '--write' in sys.argv else verify()
