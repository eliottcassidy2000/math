import json, sys
from pathlib import Path
D=Path(__file__).resolve().parent
src=json.loads((D/'source_interface.json').read_text())
names=src['target_names']+[55,56];old=src['target'];n=len(names);rows=old+[0,0]
flips=[(4,31),(6,22),(39,31),(39,4)]
for a,b in flips:
 i,j=names.index(a),names.index(b);rows[i]^=1<<j;rows[j]^=1<<i
# new prefix chain50->51->55->56; root39 beats prefix; prefix beats suffix.
for a in [55,56]:
 i=names.index(a)
 for j,b in enumerate(names):
  if b==a or b in [55,56] and b>a:continue
  if b in [39,50,51] or b==55:rows[j]|=1<<i
  else:rows[i]|=1<<j
x=dict(source_H=170807,target_H=512423,names=names,rows=rows,root=39,inport=4,prefix=[50,51,55,56],suffix=[a for a in names if a not in [39,50,51,55,56]],old_edge_flips=flips,birth_names=[55,56],h_B=5655,F_w=31673,scope='Four old-edge reversals and two births. Independent native compilation passed (eight complements, four temporary helpers). No reusable closure claimed.')
if '--check' in sys.argv:
 if json.loads((D/'candidate.json').read_text()) != json.loads(json.dumps(x)): raise RuntimeError('candidate reconstruction mismatch')
else:
 (D/'candidate.json').write_text(json.dumps(x,indent=2)+'\n')
print(json.dumps(x))
