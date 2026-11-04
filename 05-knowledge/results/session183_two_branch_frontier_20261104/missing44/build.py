import json,pathlib,itertools
D=pathlib.Path(__file__).parent
src=json.load(open(D/'source.json'));sn=src['suffix_names'];names=src['names'][:];edges=list(itertools.combinations(sn,2));z=next(json.loads(l)for l in open(D/'hits.jsonl')if 'target' in json.loads(l));ports=[v for i,v in enumerate(sn)if z['W']>>i&1];changes=[edges[z[k]]for k in ['a','b','c']]+[(39,v)for v in sorted(set(src['in_ports'])^set(ports))];rows=src['rows'][:];prefix=src['prefix'][:]
for a,b in changes:
 i,j=names.index(a),names.index(b);rows[i]^=1<<j;rows[j]^=1<<i
while len(prefix)<z['s']:
 x=max(names)+1;cut=sum(1<<names.index(v)for v in [39]+prefix);rows=[r|(((cut>>i)&1)<<len(rows))for i,r in enumerate(rows)]+[((1<<len(rows))-1)^cut];names.append(x);prefix.append(x)
out=dict(source_H=src['target_H'],target_H=z['target'],names=names,rows=rows,root=39,prefix=prefix,suffix_names=sn,in_ports=ports,all_changed_pairs=changes,suffix_flips=changes[:3],h=z['h'],J=z['J'],F=z['F'],root_E2=z['h']+(1<<z['s'])*z['F'],root_defect=z['J'],scope='Finite exact physical candidate; native guards require independent audit')
(D/'candidate44.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
