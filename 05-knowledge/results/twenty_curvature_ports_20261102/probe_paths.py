"""Scout actual square-sum path compatibility before selecting a proof scope."""
from pathlib import Path
from collections import Counter
from itertools import combinations
import importlib.util,json
HERE=Path(__file__).resolve().parent
source=HERE.parent/'joint_incidence_clocked_frames_20261102/path_kernel.py'
spec=importlib.util.spec_from_file_location('inherited_paths',source)
K=importlib.util.module_from_spec(spec);spec.loader.exec_module(K)
rows=[]
for n in range(25,31):
    edges,paths=K.graph_paths(n,True);pats=K.patterns(edges,paths)
    pairs=[(i,j) for i,j in combinations(range(len(paths)),2)
           if not(pats[i][0]&pats[j][1] or pats[i][1]&pats[j][0])]
    deg=Counter(v for e in pairs for v in e)
    unseen=set(range(len(paths)));components=[]
    while unseen:
        start=min(unseen);part={start};todo=[start]
        while todo:
            v=todo.pop()
            for i,j in pairs:
                if v not in (i,j):continue
                u=j if v==i else i
                if u not in part:part.add(u);todo.append(u)
        unseen-=part
        pair_data=[(i,j,((pats[i][0]|pats[i][1]) ^ (pats[j][0]|pats[j][1])).bit_count()//2)
                   for i,j in pairs if i in part]
        components.append({'vertices':sorted(part),'pair_d':pair_data})
    spectrum,visits=K.moment_spectrum(pats,len(edges))
    row={'n':n,'edges':len(edges),'directed_paths':len(paths),'compatible_pairs':len(pairs),
         'degree_histogram':dict(Counter(deg[i] for i in range(len(paths)))),
         'guard_spectrum':spectrum,'clique_visits':visits,'components':components}
    rows.append(row);print(json.dumps(row,sort_keys=True),flush=True)
(HERE/'path_scout.json').write_text(json.dumps(rows,sort_keys=True,indent=2)+'\n',encoding='utf-8')
