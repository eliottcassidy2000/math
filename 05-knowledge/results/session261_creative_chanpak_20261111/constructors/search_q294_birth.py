from module_tools import *
import json
exec(Path(__file__).with_name('explore_symmetric_reselection.py').read_text().split('out=[]\nfor ports')[0].split('def graph')[1].join(['def graph','']))
T=(6,4,8,19,103,79,15);hits=[];near=[]
for b in range(128):
 Z=tuple(T[i]|((1<<7) if b>>i&1 else 0) for i in range(7))+(sum(1<<i for i in range(7) if not b>>i&1),)
 h,rs=pairdata(Z)
 for r in rs:
  if r['U']==294:
   d={'orientation_mask':b,'rows':Z,'H':h,**r};hits.append(d);print(d,flush=True)
Path(__file__).with_name('q294_births.json').write_text(json.dumps({'source_rows':T,'extensions_tested':128,'selected_subsets':'all pairs','target_q':294,'hits':hits},indent=2))
print('END',len(hits))
