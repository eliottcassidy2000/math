#!/usr/bin/env python3
import argparse,importlib.util,json,pathlib
HERE=pathlib.Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('independent',HERE/'independent.py');ind=importlib.util.module_from_spec(spec);spec.loader.exec_module(ind)
def main():
 p=argparse.ArgumentParser();p.add_argument('--producer',required=True);p.add_argument('--output',required=True);a=p.parse_args();prod=json.loads(pathlib.Path(a.producer).read_text());own=json.loads((HERE/'normal.json').read_text());selected=prod['q1_classes']+[x for x in prod['q3_candidates'] if x['error']]+prod['q7_classes']+prod['q9_classes'];matches=[]
 for record in selected:
  r=tuple(record['rows']);g=ind.isomorphisms(r,r);sizes=sorted(map(len,ind.orbits(g)));ind.need(len(g)==record['full_automorphism_order'],'producer full group independently verified');ind.need(sizes==record['vertex_orbit_sizes'],'producer full orbit sizes independently verified')
  hit=[i for i,x in enumerate(own['representatives']) if ind.isomorphisms(r,tuple(x['rows']),True)];ind.need(len(hit)==1,'producer class matches exactly one independently generated class');matches+=hit
 ind.need(sorted(matches)==list(range(len(own['representatives']))),'producer and independent full class bijection')
 ind.need(prod['nonrigid_unpointing_deficit']==own['pointing_deficit'],'deficit agrees');ind.need(prod['nonrigid_pointed_states']==own['nonrigid_pointed_types'],'pointed mass agrees')
 out=dict(status='PASS',class_bijection=matches,classes=len(matches),pointing_deficit=own['pointing_deficit'],independent_target_count_input=False);pathlib.Path(a.output).write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
