#!/usr/bin/env python3
import json,subprocess,tempfile,sys
from pathlib import Path
from freeze_and_verify import adjacency,two_lane_count,validate
P=Path(__file__).resolve().parent
rows,blocks=adjacency([['D',8],['S',0]])
controls=[]
def reject(name,fn):
 try:fn()
 except ValueError:controls.append(name)
 else:raise ValueError('accepted '+name)
bad=rows.copy();bad[0]|=1
reject('self_loop',lambda:validate(bad))
bad2=rows.copy();bad2[0]^=1<<1
reject('one_sided_edge_mutation',lambda:validate(bad2))
reject('incomplete_partition',lambda:two_lane_count(rows,blocks[:-1]))
bad3=rows.copy();bad3[0]^=1<<4;bad3[4]^=1
reject('nonordinal_cross_edge',lambda:two_lane_count(bad3,blocks))
exe=Path(sys.argv[1])
text=(P/'adjacency_80557.txt').read_text()
with tempfile.TemporaryDirectory() as d:
 target=Path(d)/'wrong.txt';target.write_text(text.replace('19 80557','19 80559',1))
 p=subprocess.run([str(exe),str(target)],text=True,capture_output=True)
 if p.returncode==0 or 'count mismatch' not in p.stderr:raise ValueError('wrong target accepted')
 controls.append('wrong_target_subset_DP')
print(json.dumps({'status':'passed','rejected':controls},sort_keys=True))
