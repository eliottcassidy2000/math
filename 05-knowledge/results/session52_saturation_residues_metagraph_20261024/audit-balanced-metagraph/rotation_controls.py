#!/usr/bin/env python3
from audit import cycle_rotation_fixed
from pathlib import Path
import argparse,json
pa=argparse.ArgumentParser();pa.add_argument('--output',default='extra-rotation-controls.json');a=pa.parse_args()
rows=[]
for n,steps in [(8,[1,3,6]),(14,[1,3,5,6,10])]:
    for step in steps:
        value,orbits=cycle_rotation_fixed(n,step)
        rows.append({'n':n,'step':step,'fixed':value,'edge_orbits':orbits})
Path(a.output).write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps(rows,indent=2))
