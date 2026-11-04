import json,itertools
from pathlib import Path
p=Path(__file__).parent
source=json.load(open(p.parent/'multiport_audit'/'source_interface.json'))
names=source['target_names'];rows=source['target'];suffix_names=[1,2,4,6,7,8,14,21,22,31,32,33,40,54]
es=list(itertools.combinations(suffix_names,2));flips=[es[63],es[75]];W=[suffix_names[i] for i in range(14) if 132>>i&1]
old=rows[:]
for a,b in flips:
 i,j=names.index(a),names.index(b);rows[i]^=1<<j;rows[j]^=1<<i
r=names.index(39)
for name in suffix_names:
 i=names.index(name)
 if bool(rows[i]>>r&1)!=(name in W):rows[i]^=1<<r;rows[r]^=1<<i
changes=[(names[i],names[j]) for i,j in itertools.combinations(range(len(names)),2) if (rows[i]^old[i])>>j&1]
obj={'source_H':170807,'target_H':512425,'names':names,'rows':rows,'root':39,'prefix':[50,51],'suffix_names':suffix_names,'suffix_flips':flips,'in_ports':W,'all_changed_pairs':changes,'h':18583,'J':3690,'F':124383,'root_E2':516115,'root_defect':3690,'scope':'FINITE-EXACT strong multi-in-port +4 witness; not zero-ready; not reusable closure; native guards independently certified'}
json.dump(obj,open(p/'candidate.json','w'),indent=2)
print(json.dumps(obj))
