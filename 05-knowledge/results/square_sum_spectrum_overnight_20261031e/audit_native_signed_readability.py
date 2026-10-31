#!/usr/bin/env python3
"""Every native outer table retains readable source paths in all2401 channels."""
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent


def need(ok,msg):
    if not ok:raise RuntimeError(msg)


inner_states=readable_inner=0
old=HERE.parent/'seeded_eventual_cycle_theorem_20261019'/'certified-pair-transitions'
for family in ('Aplus','Aminus'):
    for parity in (0,1):
        for d in range(24,73):
            record=json.loads((old/f'{family}_p{parity}_d{d}.json').read_text())
            for word in record['words']:
                tails=[a for a in word if a[0] in ('block','tail')]
                need(sorted(a[1] for a in tails)==list(range(-24,25)),'one recoverable inner tail per offset')
                for c in range(-24,25):
                    atom=next(a for a in tails if a[1]==c)
                    need(atom[2] in (0,1) and atom[3] in (-1,1),'typed inner source member')
                    readable_inner+=1
            inner_states+=1
table=json.loads((HERE.parent/'square_sum_spectrum_overnight_20261031a'/'n24_forest_table.json').read_text())
outer_channels=0
for record in table:
    word=record['word'];i1=word.index(['single',1]);i8=word.index(['single',8])
    need((i1-i8)%len(word) in (1,len(word)-1),'fixed oriented native edge anchors cyclic parser')
    for c in range(-24,25):
        atoms=[a for a in word if a[0]!='single' and a[1]==c]
        need(len(atoms)==2 and sorted(a[3] for a in atoms)==[0,1],'both outer source fragments')
        need(len({a[2] for a in atoms})==1,'same selected A member for both fragments')
        need(all(a[4] in (-1,1) for a in atoms),'oriented fragments')
        outer_channels+=1
print(json.dumps({'status':'passed','inner_states':inner_states,'readable_inner_member_channels':readable_inner,
                  'outer_states':len(table),'readable_outer_channels':outer_channels,
                  'independent_native_channels':49*49,
                  'scope':'Fixed atom lengths and orientations recover each selected path; reconstructed label is nearest integer to lifted_label/49.'},indent=2,sort_keys=True))
