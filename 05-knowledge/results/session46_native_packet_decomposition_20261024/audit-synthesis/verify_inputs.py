#!/usr/bin/env python3
"""Verify explicit reviewed inputs and cross-check receipt claims, without replaying producers."""
import argparse
from hashlib import sha256
import json
from pathlib import Path


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser()
    where=parser.add_mutually_exclusive_group(required=True)
    where.add_argument('--workspace');where.add_argument('--result-root')
    parser.add_argument('--repo',required=True);parser.add_argument('--out',required=True)
    args=parser.parse_args()
    here=Path(__file__).resolve().parent
    pins=json.loads((here/'SOURCE_PINS.json').read_text())
    checks=0
    def require(ok,description):
        nonlocal checks
        checks+=1
        if not ok: raise RuntimeError(description)
    def locate(scope,path):
        if scope=='repository': return Path(args.repo)/path
        if args.result_root: return Path(args.result_root)/path
        first,slash,rest=path.partition('/')
        if first in ('SYNTHESIS.md','GUARD_AND_REUSE_CONTRACT.md'):
            return Path(args.workspace)/'session46-root'/path
        return Path(args.workspace)/('session46-'+first)/rest
    for item in pins['reviewed_sources']:
        path=locate(item['scope'],item['path'])
        require(path.is_file(),'missing reviewed input '+item['path'])
        require(digest(path)==item['sha256'],'changed reviewed input '+item['path'])
    reports=[]
    pairs=set()
    for rel in pins['packet_manifests']:
        mp=locate('result',rel);manifest=json.loads(mp.read_text())
        entries=manifest.get('files',manifest.get('artifacts'))
        if isinstance(entries,dict):
            entries=[dict(path=k,sha256=v) for k,v in entries.items()]
        for entry in entries:
            p=mp.parent/entry['path']
            require(p.is_file(),'missing manifest artifact '+str(p))
            require(digest(p)==entry['sha256'],'manifest digest mismatch '+str(p))
            if 'bytes' in entry:
                require(p.stat().st_size==entry['bytes'],'manifest size mismatch '+str(p))
            if 'normal' in p.name and p.suffix=='.json':
                other=p.with_name(p.name.replace('normal','optimized'))
                if other.is_file():
                    pairs.add((p,other))
        reports.append({'manifest':rel,'artifacts_verified':len(entries)})
    for a,b in sorted(pairs):
        require(a.read_bytes()==b.read_bytes(),'normal/optimized receipt mismatch '+str(a))
    def read(path): return json.loads(locate('result',path).read_text())
    producer=read('larger-circuits/cost-optimal18/normal.json')['circuits'][0]
    native=read('audit-larger/cost-optimal18/normal.json')['circuits'][0]
    for field,value in [('balanced_choices',48620),('closed_ports_per_endpoint',72),
                        ('closed_ports_union',144),('cost',37),('frozen_X_debt_at_U',0),
                        ('full_H256_accepted_choices',2),('full_H8_accepted_choices',754),
                        ('joint_choices',262144)]:
        require(producer[field]==value,'producer circuit claim '+field)
    for pfield,afield in [('balanced_choices','balanced_states'),
                          ('closed_ports_per_endpoint','closed_ports_per_endpoint'),
                          ('closed_ports_union','closed_ports_union'),
                          ('cost','ordinary_owner_clock_cost'),
                          ('frozen_X_debt_at_U','integer_X_target_debt'),
                          ('full_H256_accepted_choices','H256_accepted'),
                          ('full_H8_accepted_choices','H8_accepted'),
                          ('joint_choices','all_joint_states'),('plus','plus'),('minus','minus')]:
        require(producer[pfield]==native[afield],'producer/independent native mismatch '+pfield)
    original=read('audit-larger/normal.json')['circuits'][0]
    require((original['ordinary_owner_clock_cost'],original['H8_accepted'],
             original['H256_accepted'],abs(original['integer_X_target_debt']))==(43,792,2,512),
            'original eighteen-source claim')
    boundary=read('audit-larger/cost-optimal18/cost-boundary.normal.json')['results']
    for r in boundary:
        cost=r['cost']
        require(r['equal_charge_relations']=={36:1129,37:69358}[cost],'boundary relation count')
        require(len(r['canonical_primitive_partitions'])=={36:0,37:1}[cost],'boundary primitive count')
    moves=read('move-lattices/normal.json')
    require((moves['four_slot_moves'],moves['six_slot_moves'],moves['four_move_rank'],
             moves['integer_move_lattice_index'])==(27,1917,23,2048),'move lattice summary')
    require(moves['low_cardinality_fibers'][2]['component_distribution']==
            {'1':1851,'2':175,'3':21,'4':1},'cardinality-four components')
    require(moves['binary_hostile']['selected_sources']==
            [[176,1425,1549,2087],[1394,1923,2837,3369]],'full28 hostile source sets')
    require(sum((1851,175,21,1))==2048 and 175+21+1==197,'component census arithmetic')
    routes=read('audit-compensation/optimization.normal.json')['routes']
    expected={'buffer_3015':(6,4,50,13,4,50),'buffer_1394':(8,3,40,16,3,40),
              'buffer_1394_2905':(6,4,48,12,4,48),
              'eighteen_buffer_1394':(6,4,55,17,4,55),
              'eighteen37_buffer_2071':(6,4,49,12,5,57)}
    for route in routes:
        if route['name'] not in expected: continue
        s=route['at_optimal_slot_budget'];c=route['at_optimal_clock_budget']
        actual=(s['threshold'],s['minimum_events'],s['minimum_total_clock_work'],
                c['threshold'],c['minimum_events'],c['minimum_total_clock_work'])
        require(actual==expected[route['name']],'route tradeoff '+route['name'])
    theory=read('decomposition-theory/normal.json')['davenport']
    require(theory['conditional_native_support_bound']==264 and
            theory['abstract_minimal_zero_sum_length']==132 and
            not theory['native_sharpness_claimed'] and not theory['all_order_constant_claimed'],
            'Davenport scope receipt')
    report={'status':'PASS','checks':checks,'reviewed_sources':len(pins['reviewed_sources']),
            'packet_manifests':reports,'manifest_artifacts':sum(r['artifacts_verified'] for r in reports),
            'normal_optimized_pairs':len(pairs),
            'source_pins_sha256':digest(here/'SOURCE_PINS.json'),
            'scope':'Input hashes, byte equality and synthesis receipt cross-checks only; not producer replays'}
    Path(args.out).write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps(report,sort_keys=True))

if __name__=='__main__': main()
