#!/usr/bin/env python3
"""Normalized native tile capacity versus literal owner-labelled resources."""
import argparse
import csv
from hashlib import sha256
import importlib.util
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('literal_resource_checker',HERE/'check.py')
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)

def word(port):
    a,p=port;n=len(a)
    pairs=[(i,j) for i in range(n) for j in range(i+2,n)]
    return sum(a[p[i]][p[j]]<<b for b,(i,j) in enumerate(pairs))

def tile(port):
    return sorted({word(p) for p in base.orbit(port)})

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
    pin=next(p for p in pins['files'] if p['path'].endswith('merged_nodes.tsv'))
    data=(args.repo/pin['path']).read_bytes()
    base.need(sha256(data).hexdigest()==pin['sha256'],'original canonical source atlas pinned')
    rows=list(csv.DictReader(data.decode().splitlines(),delimiter='\t'))
    base.need(int(rows[43]['canon_a'])==17301538,'source43 canonical code')
    a=base.matrix(17301538);ps=[(7,6,5,3,2,1,0,4),(7,6,5,2,1,3,0,4)]
    resources=[base.orbit((a,p)) for p in ps]
    base.need(all(base.valid(a,p) for p in ps),'both hostile paths authentic')
    base.need([word((a,p)) for p in ps]==[344062,344062],'two full lifts have the same normalized word')
    base.need(len(set().union(*resources))==8,'two distinct raw resources')
    hostile_tile=[344062,1046512,1050639,1753089]
    base.need(all(tile((a,p))==hostile_tile for p in ps),'one normalized unit resource')
    b=base.matrix(int(rows[2242]['canon_a']))
    hall_paths=[(1,0,5,3,2,4,6,7),(1,5,4,3,2,0,6,7),(4,1,5,3,2,0,6,7)]
    hall_keys=[min(tile((b,p))) for p in hall_paths]
    base.need(len(set(hall_keys))==3,'native three-port target fibre has three distinct normalized tiles')
    base.need(sum(p[0]==1 for p in hall_paths)==2,'native shared pins leave two normalized resources for three demands')
    packet_keys={}
    for filename in ('packet-four-source.json','packet-three-source.json'):
        packet=json.loads((HERE/filename).read_text());keys=[]
        for slot in packet['slots']:
            physical=tuple(map(tuple,slot.get('physical_matrix',slot.get('adjacency'))))
            names=('before','after') if 'physical_matrix' in slot else ('old','new')
            for name in names:
                r=slot[name];p=tuple(r.get('physical_path',r['path']))
                base.need(base.valid(physical,p),'packet normalized from actual physical path')
                base.need(word((physical,p))==r['mask'],'packet uses original forward word')
                keys.append(min(tile((physical,p))))
        base.need(len(keys)==len(set(keys))==8,'all eight endpoint resources distinct after native normalization')
        packet_keys[filename]=keys
    out=dict(status='FINITE-EXACT PASS',checks=base.CHECKS,quotient_hostile=dict(source=43,canon=17301538,paths=ps,word=344062,raw_orbits=2,raw_port_union=8,normalized_tile=hostile_tile,normalized_unit_capacity=1),native_hall_hostile_normalized_keys=hall_keys,packet_normalized_keys=packet_keys)
    args.out.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(dict(status=out['status'],checks=base.CHECKS,receipt_sha256=sha256(args.out.read_bytes()).hexdigest())))

if __name__=='__main__':main()
