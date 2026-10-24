"""Complete boundary-conditioned charge sections for the actual XOR tree."""
from itertools import combinations
from pathlib import Path
from hashlib import sha256
import json,sys
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parents[1];CHECKS=0
def need(v,m):
    global CHECKS
    CHECKS+=1
    if not v:raise RuntimeError(m)
def main():
    supply=RESULTS/'compatible_charged_helpers_20261023/supply/selected_supply.json';native=HERE/'native_2125_receipt.json'
    d=json.loads(supply.read_text());n=json.loads(native.read_text());by={z['source']:z for z in d['helpers']}
    ids=next(b['source_ids']for b in d['banks']if len(b['source_ids'])==53)
    labels={s:z['label']for s,z in by.items()};labels[2125]=n['label']
    def charge(ss):
        low=0
        for s in ss:low^=sum(v<<j for j,v in enumerate(labels[s][1:]))
        return((sum(labels[s][0]for s in ss)%128)<<4)|low
    rows=[]
    for expanded in(False,True):
        supplyids=ids+([2125]if expanded else[]);free=[s for s in supplyids if s not in(368,776,3139)]
        for bit in(0,1):
            fixed=[368]if bit else[776,3139];k=4-len(fixed);section={};count=0
            for ss in combinations(free,k):
                q=charge(fixed+list(ss));section.setdefault(q,list(ss));count+=1
            need(len(section)=={(False,0):946,(False,1):2047,(True,0):972,(True,1):2048}[expanded,bit],
                 'complete native phase-conditioned coverage count')
            rows.append(dict(expanded54=expanded,center_V=bit,fixed_V=fixed,free_source_count=len(free),
                             free_V_count=k,subset_count=count,charge_count=len(section),missing=sorted(set(range(2048))-set(section))))
            if expanded and bit:
                complete=section
                for q,ss in complete.items():need(len(set(ss))==3 and set(ss)<=set(free)and charge([368]+ss)==q,'literal all-charge fixed-phase section')
    costs={s:z['selected']['cost']for s,z in by.items()};costs[2125]=8
    six=sorted((costs[s],s)for s in ids+[2125]if s not in(368,776,3139))[-6:]
    need(sum(c for c,s in six)==43,'actual six-slot phase-fixed clock diameter bound')
    out=dict(status='FINITE-EXACT actual fixed-shared-owner-phase section',rows=rows,
             fixed_V=[368],fixed_U=[776,3139],shared_owner_value=3,full_helper_occupancy=4,
             fixed_phase_section={q:complete[q]for q in sorted(complete)},reset_support_bound=6,reset_clock_bound=43,
             largest_six_cost_slots=six,checks=CHECKS,
             inputs={p.name:sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()for p in(supply,native)})
    out=json.loads(json.dumps(out));p=HERE/'phase_sections_receipt.json'
    if '--write'in sys.argv:p.write_text(json.dumps(out,indent=2)+'\n',newline='\n')
    else:need(out==json.loads(p.read_text()),'frozen complete boundary section')
    print(f'PASS {CHECKS}; fixed centerV covers2048, fixed centerU only972; reset6slots/43clocks')
if __name__=='__main__':main()
