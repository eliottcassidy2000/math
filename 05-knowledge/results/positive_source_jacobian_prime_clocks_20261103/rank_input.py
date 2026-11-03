from pathlib import Path
from carrier import forest_groups
import json
here=Path(__file__).resolve().parent
edges,groups,_=forest_groups(6)
lines=[str(len(groups))];meta=[]
for (isolated,pairs),terms in groups.items():
    blocks=[(v,v,1) for v in isolated]+[(a,b,(-1)**(parity+1)) for a,b,parity in pairs]
    lines.append(' '.join(map(str,[len(blocks)]+[z for b in blocks for z in b])))
    meta.append(dict(isolated=isolated,pairs=pairs,terms=terms))
(here/'rank_input.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')
(here/'rank_types.json').write_text(json.dumps(meta,separators=(',',':'))+'\n',encoding='utf-8')
