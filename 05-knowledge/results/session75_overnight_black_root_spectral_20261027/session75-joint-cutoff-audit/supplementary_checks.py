#!/usr/bin/env python3
"""Clock-border, metadata hostiles and scoped producer-code replay."""
import argparse
import copy
import importlib.util
import json
from pathlib import Path
import independent_audit as own

def load_module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--packet',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();root=args.packet
    V={int(line.split()[0]) for line in (root/'frontier.txt').read_text().splitlines() if line.strip() and not line.startswith('#')}
    rows,_=own.reconstruct(4096,V)
    before=own.direct(rows,4096,36,64);after=own.direct(rows,4096,36,66)
    arrivals=[r for r in rows if r['d']<=36 and 64<r['t']<=66]
    saved=own.read(root/'results/N4096_clock_border_mechanism.json')
    own.check(len(arrivals)==1 and arrivals[0]==saved['added_source'] and arrivals[0]['s']==13759,'unique clock-border source')
    z=arrivals[0]; interactions=[]
    for x in rows:
        if x['u']==z['u'] and x['d']<=36 and x['t']<=64:
            gap=x['m']^z['m']
            interactions.append({'source':x,'q':gap&-gap,'beta':own.B(x['v'],z['v'])})
    own.check(interactions==saved['same_terminal_pair_contributions_including_zero'],'every clock-border pair')
    nonzero=[(p['source']['s'],p['q'],p['beta']) for p in interactions if p['beta']]
    own.check(nonzero==[(9,1,-1),(2831,8,1),(10065,1,-1)],'three nonzero interactions')
    change=[y-x for x,y in zip(before['level_deltas'],after['level_deltas'])]
    own.check(change==[-2,0,0,1]+[0]*9 and (before['D_plus'],after['D_plus'])==(495,494),'clock-border exact level and debt change')
    def rejected(altered):
        try:own.check(altered==rows,'independent source reconstruction rejects source mutation')
        except RuntimeError:return True
        return False
    mutations={}
    for field in ('s','u','p','t','d','m','A','v'):
        altered=copy.deepcopy(rows)
        if field in ('A','v'):altered[-1][field][0]+=1
        else:altered[-1][field]+=2 if field=='t' else 1
        mutations[field]=rejected(altered)
    zero=next(r for r in rows if own.Q(r['A'])==0)
    moving=next(r for r in rows if r['d']>36 and own.Q(r['A'])>0)
    mutations['deleted_zero']=rejected([r for r in rows if r['s']!=zero['s']])
    mutations['deleted_moving']=rejected([r for r in rows if r['s']!=moving['s']])
    own.check(all(mutations.values()),'all complete-record hostiles rejected')
    original_norm=sum(own.Q(r['A']) for r in rows)
    mutated_norm=sum(own.Q(r['A']) for r in rows if r['s']!=moving['s'])
    own.check(original_norm==9949>mutated_norm,'moving-source denominator loss')
    # Reuse is explicitly reported here: these are the producer's implementations,
    # executed directly without writing to, copying, or mutating their packet.
    scan=load_module(root/'scan.py','producer_scan')
    literal=load_module(root/'literal_verify.py','producer_literal')
    rr,_=scan.records(4096,scan.frontier()); sweep=scan.sweep(rr,4096)
    recorded=own.read(root/'results/receipt_N4096.json')
    for key,value in sweep.items():
        own.check(recorded[key]==value,'producer scan replay '+key)
    reconstructed,h,peak=literal.reconstruct(4096)
    own.check(reconstructed==rows,'producer literal replay complete fields')
    expected=own.read(root/'results/literal_receipt_N4096.json')
    own.check(h==expected['detailed_source_routes_and_Bchains_sha256'] and peak==expected['max_exterior_state'],'producer route-stream replay')
    pair=literal.pair_birth_verify(rows,4096,recorded)
    own.check(pair==expected['all_rectangle_pair_birth_check'],'producer pair recurrence replay')
    witness=own.read(Path(__file__).resolve().parent/'LARGE_INTEGER_WITNESS.json')
    x=witness['s']; route=[]; clock=0
    while x not in V:
        y,dt=own.step(x);route.append({'x':x,'next':y,'edge_clock':dt});clock+=dt;x=y
    own.check(route==witness['route'] and x==witness['terminal'] and clock==witness['clock'],'full high-state witness route')
    own.check(max(z['x'] for z in route)==28167930305>2**32-1 and len(route)==116 and route[-1]['x']==9827,'high-state exact integer safeguard')
    out={'status':'PASS','packet_manifest_sha256':own.digest(root/'MANIFEST.json'),
         'clock_border':{'added_source':z,'all_pairs_checked':len(interactions),'nonzero':nonzero,'delta_change':change,'before_debt':495,'after_debt':494},
         'metadata_mutations_rejected':mutations,'denominator_control':{'original':original_norm,'wrong_after_moving_deletion':mutated_norm},
         'large_integer_witness':{'source':witness['s'],'route_edges_checked':len(route),'maximum_exterior_state':28167930305},
         'producer_replay':{'N':4096,'scope':'Producer scan full surface, literal full source reconstruction and pair-birth full surface; separate from independent reconstruction. Not a replay of the all-19-scale run_all.sh.','rectangles':5917,'pair_report':pair}}
    args.out.write_bytes(own.enc(out));print(json.dumps(out,sort_keys=True))

if __name__=='__main__':main()
