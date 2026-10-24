#!/usr/bin/env python3
"""Reference comparison only; never imported by independent audit.py."""
from pathlib import Path
from hashlib import sha256
import argparse,json
P=Path(__file__).resolve().parent
N=0
def need(ok,why):
    global N
    N+=1
    if not ok:raise RuntimeError(why)
parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
for entry in json.loads((P/'REFERENCE_ARCHIVES.json').read_text()):
    raw=(P/entry['archive']).read_bytes();need(len(raw)==entry['bytes'] and sha256(raw).hexdigest()==entry['sha256'],'reference archive bytes')
ours=json.loads((P/'receipt.json').read_text());source=P/'archives/producer_certificate.json';theirs=json.loads(source.read_text())
need(ours['base_sources']==theirs['base_sources'],'identical receiver inventory')
need(ours['pair_sources']==theirs['old_pair_sources'],'identical linked pair')
need(ours['old_pair_holes']==theirs['old_pair_holes'],'identical linked18-slot holes')
allbits=(1<<2048)-1;tables=0
for row in theirs['table_census']:
    p=str(row['source']);table=ours['relation_census'][p]['tables'][str(row['allowed_table_mask'])]
    need(table['complete_grades']==row['complete_grades'],'all producer complete-grade lists')
    need([r['count'] for r in table['grades']]==row['coverage'],'all producer grade coverages')
    need([r['first_hole'] for r in table['grades']]==row['first_missing'],'all producer hostile first holes')
    bitmap=b''.join((allbits^int(r['hole_bitmap_hex'],16)).to_bytes(256,'little') for r in table['grades'])
    need(sha256(bitmap).hexdigest()==row['support_bitmap_sha256'],'every actual producer support bitmap identical')
    need(ours['costs'][p]==row['cost'],'producer actual source cost');tables+=1
need(tables==105,'all105 producer nonempty tables compared')
labels=ours['labels'];base=ours['base_sources']
def plus(q,r):return 16*((q//16+r//16)%128)+((q^r)&15)
def charge(mask):
    q=0
    for i,s in enumerate(base):
        if mask>>i&1:q=plus(q,labels[str(s)])
    return q
pairs=0
for row in theirs['triple3015_exact_boundaries']:
    grade=row['grade'];own=ours['exact_branch_minima'][str(grade)]
    need(own['histogram']==row['exact_clock_histogram'],'all exact branch boundary histograms')
    minima={r['charge']:r for r in own['minimizing_pairs']}
    need(len(minima)==row['both_branch_charges']==len(row['minimizing_witnesses']),'same entire both-branch charge universe')
    for record in row['minimizing_witnesses']:
        q=record['charge'];a=record['U_branch_base_mask'];b=record['V_branch_base_mask']
        need(record['clock_minimum']==minima[q]['clocks'] and record['slot_minimum']==minima[q]['slots'],'every producer exact minimum independently recovered')
        need(0<=a<1<<16 and 0<=b<1<<16 and a.bit_count()==grade and b.bit_count()==grade-3,'producer minimizing witnesses actual exact grades')
        need(charge(a)==q and plus(charge(b),1263)==q,'producer minimizing witnesses actual charge')
        need(8+2*(a^b).bit_count()==record['clock_minimum'] and 3+(a^b).bit_count()==record['slot_minimum'],'producer minimizing witnesses attain actual native cost');pairs+=1
need({tuple(x) for x in ours['sharp_receiver_V_sets']}=={tuple(x) for x in theirs['sharp_fixed_exterior_fibre']['states_V']},'same sharp two-state fibre')
result={'status':'PASS','compared_tables':tables,'compared_support_grades':20*tables,'compared_exact_boundary_minima':pairs,'checks':N,'producer_certificate_sha256':sha256(source.read_bytes()).hexdigest(),'independent_receipt_sha256':sha256((P/'receipt.json').read_bytes()).hexdigest(),'scope':'All original producer support bitmaps, full complete-grade lists, first holes, source costs, all per-charge boundary minima and supplied minimizing masks match the independent audit.'}
text=json.dumps(result,indent=2,sort_keys=True)+'\n'
if args.write:(P/'PRODUCER_COMPARISON.json').write_text(text)
else:need(json.loads((P/'PRODUCER_COMPARISON.json').read_text())==result,'frozen comparison receipt exact')
print(json.dumps(result,sort_keys=True))
