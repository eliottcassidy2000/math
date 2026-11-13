#!/usr/bin/env python3
"""Complete entry test for the six finite determinant outputs, fixed refill seed."""
from hashlib import sha256
from pathlib import Path
import json
from refill import AP, count, select_target, require

def main():
    base=Path(__file__).resolve().parent
    source=base.parent/'pair_growth2.json'
    payload=source.read_bytes(); data=json.loads(payload)
    require([] in select_target(count(1,AP))['addresses'],'positive selector control')
    entries=[]
    for hit in data['hits']:
        h,p=hit['M']; j,q=hit['N']
        target=36*p*q+6*h*q+18*p*j-h*j
        addresses=[]; examined=0; j0=1
        while 240*3**j0*128<=target:
            result=select_target(target,j0)
            examined+=result['bounded_states_examined']
            addresses.extend({'j0':j0,'word':w} for w in result['addresses'])
            j0+=1
        entries.append({'M':hit['M'],'N':hit['N'],'target':target,
                        'clock_addresses':addresses,'initial_j_tested':[1,j0-1],
                        'word_states_examined':examined})
    receipt={'status':'FINITE-EXACT',
             'scope':'Complete fixed A-seed wrapped-refill entry test with all j0>=1; not all clock modules',
             'source':'../pair_growth2.json','source_sha256':sha256(payload).hexdigest(),
             'entries':entries,
             'reproduce':'python3 -B 05-knowledge/results/symmetry_port_return_20261113/clock/entry_probe.py'}
    (base/'determinant_entry.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'targets':len(entries),'word_states':sum(x['word_states_examined'] for x in entries),
                      'hits':sum(len(x['clock_addresses']) for x in entries)}))

if __name__=='__main__': main()
