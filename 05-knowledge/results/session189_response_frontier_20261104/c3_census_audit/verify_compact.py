#!/usr/bin/env python3
"""Verify compact replay with no tower exports or producer source dependency.
Usage: verify_compact.py CENSUS_TSV PRODUCER_CANDIDATES_JSON AUDIT_OUTPUT_DIR
"""
import csv, hashlib, json, pathlib, sys

def need(ok,message):
    if not ok:
        raise RuntimeError(message)
need(len(sys.argv)==4,'usage: verify_compact.py CENSUS_TSV CANDIDATES_JSON OUTPUT_DIR')
census,original,out=map(pathlib.Path,sys.argv[1:])
producer=json.loads(original.read_text())
def table(path):
    with path.open() as f:
        return list(csv.DictReader(f,delimiter='\t'))
census_rows=table(census)
full=table(out/'all_rows.tsv')
need(len(full)==len(census_rows)==97144,'Census row count mismatch')
for i,(a,b) in enumerate(zip(census_rows,full)):
    need(int(a['node'])==int(b['node'])==i,'ID mismatch')
    need(int(a['code'])==int(b['source_code']),'Source code mismatch')
    need(int(b['n'])==9,'Wrong order')
summary=table(out/'summary.tsv')
need(len(summary)==1,'Expected one order in compact replay')
s={k:int(v) for k,v in summary[0].items()}
need(s['n']==9,'Wrong summary order')
need(s['labeled_mass']==s['expected_mass']==2**36,'Incomplete orbit mass')
need(s['ordinary']==2*s['merged']-s['self_converse'],'Incorrect merged partition')
for field,key in [('merged','total_merged_nodes'),('prime_merged','prime_nodes'),('p_nonzero','p_nonzero'),('p_zero','p_zero')]:
    need(s[field]==producer[key],f'Producer count mismatch: {field}')
candidates=[]
for x in table(out/'candidates.tsv'):
    candidates.append(dict(node=int(x['node']),code=int(x['code']),
                           rows=list(map(int,x['rows'].split(','))),
                           starts=list(map(int,x['starts'].split(','))),
                           ends=list(map(int,x['ends'].split(',')))))
expected=[{**x,'code':int(x['code'])} for x in producer['candidates']]
need(candidates==expected,'Candidate ID/code/row/start/end mismatch')
need(len(candidates)==s['p_zero'],'Candidate count mismatch')
keys=set();mass=0;primes=nonzero=zero=0
for r in full:
    aut=int(r['aut']);need(aut>0 and 362880%aut==0,'Invalid automorphism divisor')
    k,d=int(r['key']),int(r['converse_key'])
    need((k==d)==bool(int(r['self_converse'])),'Self-converse mismatch')
    for key in ([k] if k==d else [k,d]):
        need(key not in keys,'Duplicate ordinary key');keys.add(key);mass+=362880//aut
    prime=int(r['prime']);need(prime in (0,1),'Invalid prime flag')
    if prime:
        primes+=1
        starts=list(map(int,r['starts'].split(',')));ends=list(map(int,r['ends'].split(',')))
        need(len(starts)==len(ends)==9,'Wrong endpoint lengths')
        need(sum(starts)==sum(ends),'Start/end sum mismatch')
        mask=sum(((a+b)%2)<<i for i,(a,b) in enumerate(zip(starts,ends)))
        need(mask==int(r['p_mask']),'Parity mask mismatch')
        if mask:nonzero+=1
        else:zero+=1
need(len(keys)==s['ordinary'] and mass==2**36,'Row certificate mass mismatch')
need((primes,nonzero,zero)==(s['prime_merged'],s['p_nonzero'],s['p_zero']),'Row filter count mismatch')
paths=[census,original,pathlib.Path(__file__),pathlib.Path(__file__).with_name('audit.cpp'),
       out/'all_rows.tsv',out/'summary.tsv',out/'candidates.tsv']
receipt={'status':'PASS','summary':s,'candidate_rows_checked':len(candidates),
         'independent_label_mass':mass,'ordinary_keys_checked':len(keys),
         'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
(out/'receipt.compact.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ('status','candidate_rows_checked','independent_label_mass','ordinary_keys_checked')},indent=2))
