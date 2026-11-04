from pathlib import Path
import json, sys
D=Path(__file__).resolve().parent;s=json.loads((D.parent/'candidate.json').read_text());names=[x for x in s['names'] if x not in [55,56]];idx=[s['names'].index(x) for x in names];rows=[sum(((s['rows'][i]>>j)&1)<<k for k,j in enumerate(idx)) for i in idx]
edges=[(i,j) for i in range(14) for j in range(i+1,14)];suffix=s['suffix'];flips=[tuple(suffix[i] for i in edges[k]) for k in [43,59]]+[(39,4),(39,8)]
for a,b in flips:
 i,j=names.index(a),names.index(b);rows[i]^=1<<j;rows[j]^=1<<i
x=dict(source_H=512423,target_H=1537271,names=names,rows=rows,root=39,inport=8,prefix=[50,51],suffix=suffix,old_edge_flips=flips,deleted_names=[55,56],preserved_original_names=names,h_B=100007,F_w=359316,scope='Deletion-allowed continuation. Deletes the prior two prefix births, retains original17 names. No closed return claimed.')
candidate=x
print(json.dumps(candidate))
pairs={tuple(map(int,l.split()[1:3])) for l in (D.parent/'prefixes.txt').read_text().splitlines()};ls=[json.loads(l) for l in (D/'hits.jsonl').read_text().splitlines()];ls=[x for x in ls if x.get('lattice')];large=[x for x in ls if x['P2A']>=256]
if len(ls)!=16 or max(x['P2A'] for x in ls)!=336: raise RuntimeError('candidate universe changed')
if len(large)!=1 or not large[0]['P2A']<9*large[0]['hA']: raise RuntimeError('order8 exclusion failed')
match=[x for x in ls if (x['hA'],x['P2A']) in pairs];print(match)
if len(match)!=1: raise RuntimeError('unique realization failed')
v=match[0]
if tuple(v[k] for k in ('flip_a','flip_b','mark','hA','P2A','target')) != (43,59,5,1,4,1537271): raise RuntimeError('wrong survivor')
classification=dict(lattice_candidates=len(ls),max_P2=max(x['P2A'] for x in ls),excluded_order8=large,unique_match=match)
for name,value in [('candidate.json',candidate),('classification.json',classification)]:
 if '--check' in sys.argv:
  if json.loads((D/name).read_text()) != json.loads(json.dumps(value)): raise RuntimeError(name+' reconstruction mismatch')
 else:
  (D/name).write_text(json.dumps(value,indent=2)+'\n')
