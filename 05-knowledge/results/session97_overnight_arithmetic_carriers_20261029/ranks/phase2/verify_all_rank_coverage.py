"""Join independent F2/F3 verifiers and check coverage of every claimed point list."""
import json,subprocess,sys
from pathlib import Path
P=Path(__file__).resolve().parent
outputs={}
for name in ['verify_highrank_snapshot.py','verify_odd_primary_certificates.py']:
 cmd=[sys.executable]+(['-O']if not __debug__ else[])+[str(P/name)]
 outputs[name]=json.loads(subprocess.check_output(cmd,text=True,cwd=P))
cs={c['id']:c for c in json.loads((P/'highrank_inputs.json').read_text())['curves']}
binary={z['curve_id']:z for z in json.loads((P/'all_highrank_character_certificates.json').read_text())}
ternary={z['curve']:z for z in json.loads((P/'odd_primary_certificates.json').read_text())}
if len(cs)!=66 or set(binary)!=set(cs):raise ValueError('snapshot coverage differs')
full2=full3=total=0
for i,c in cs.items():
 n=len(c['points'])
 if n!=c['rank_lower_bound']:raise ValueError('point count differs from rank claim')
 if binary[i]['character_rank']==n:full2+=1
 elif i in ternary and ternary[i]['rank']==n:full3+=1
 else:raise ValueError('unproved claimed rank')
 total+=n
print(json.dumps({'curves':len(cs),'binary_full':full2,'ternary_repairs':full3,'independent_point_directions':total,'verifier_checks':sum(z['checks']for z in outputs.values())},sort_keys=True))
