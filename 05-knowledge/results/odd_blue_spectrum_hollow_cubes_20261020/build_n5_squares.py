#!/usr/bin/env python3
"""Extract exact two-square data after verify_odd_pairing.py's inherited check."""
import json
from pathlib import Path
p=Path(__file__).resolve().parent
d=json.loads((p/'pairing_n5.json').read_text()); packets=[]
for root in sorted(set(r['cube_packet'] for r in d['rows'])):
    rows=[r for r in d['rows'] if r['cube_packet']==root]
    ids={r['mask']:r['id'] for r in rows}; frame=rows[0]['mask']&~3
    packet={'packet_id':root,'hollow_frame_mask':frame,
            'blue_degree':rows[0]['blue'],'hollow_group_order':rows[0]['hollow_aut'],
            'pair_orbits':[[0],[1]],'vertices':[],'edges':[]}
    for r in rows:
        assert 'inherited_node' in r, 'Rerun verifier with --inherited pointing to the full census directory.'
        packet['vertices'].append({
            'completion_bits':[(r['mask']>>i)&1 for i in range(2)],
            'frame_class_id':r['id'],'normalized_frame_canonical_mask':r['mask'],
            'inherited_census_node':r['inherited_node'],
            'inherited_explorer_representative':r['inherited_representative'],
            'blue_degree':r['blue'],'total_mass':r['mass'],'c3':r['c3'],
            'type':r['type'],'H':r['H_total'],'Aut':r['aut_total']})
        for i in range(2):
            target=ids[r['mask']^(1<<i)]
            if r['id']<target:
                packet['edges'].append({'from':r['id'],'to':target,'toggled_pair_orbit':[i]})
    packets.append(packet)
out={'n':5,
     'label_convention':'x0,x1 -> f -> y0,y1; a_i=1 means x_i->y_i. A has x1->x0. Cross bit0: y1->x0 and y0->x1; cross bit1 reverses both. Canonical frame mask encodes [a0,a1,A01,B01].',
     'source':'pairing_n5.json, exact verifier cross-matched to census_n5.json',
     'edge_warning':'These are auxiliary matching-orbit toggle edges, NOT literal antipodal blue edges.',
     'packets':packets}
(p/'n5_auxiliary_squares.json').write_text(json.dumps(out,indent=2)+'\n')
print('PASS: two exact n5 squares, eight classes, eight auxiliary edges.')
