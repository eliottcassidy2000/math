import json,itertools,pathlib
D=pathlib.Path(__file__).parent
c=json.load(open(D/'source.json'))
names=c['names']; sn=c['suffix_names']; edges=list(itertools.combinations(sn,2))
for digit,sel in [(2,(14,83,128,3)),(4,(66,73,6144,3))]:
 hits=[json.loads(l)for l in open(D/f'hits{digit}.jsonl')];z=next(z for z in hits if tuple(z[k] for k in ['a','b','W','s'])==sel)
 rows=c['rows'][:];changes=[edges[z['a']],edges[z['b']]];ports=[v for i,v in enumerate(sn)if z['W']>>i&1];changes +=[(39,v)for v in set(c['in_ports'])^set(ports)]
 for a,b in changes:
  i,j=names.index(a),names.index(b);rows[i]^=1<<j;rows[j]^=1<<i
 # Fresh prefix55 follows prefix50,51 and beats suffix; root39 beats it.
 cut=sum(1<<names.index(v)for v in [39,50,51]);rows=[r|(((cut>>i)&1)<<len(rows))for i,r in enumerate(rows)]+[((1<<len(rows))-1)^cut]
 out=dict(source_H=512425,target_H=z['target'],names=names+[55],rows=rows,root=39,prefix=[50,51,55],suffix_names=sn,in_ports=ports,all_changed_pairs=changes,suffix_flips=changes[:2],h=z['h'],J=z['J'],F=z['F'],root_E2=z['h']+8*z['F'],root_defect=z['J'],scope='Finite exact candidate, pending independent native and phase audit')
 (D/f'candidate{digit}.json').write_text(json.dumps(out,indent=2)+'\n');print(digit,changes,ports,out['root_E2'])
