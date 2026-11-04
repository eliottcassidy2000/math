import json,pathlib,itertools
D=pathlib.Path(__file__).parent
for word,sel in [('22',(55,85,7104,3)),('24',(12,89,8,3)),('42',(66,72,704,3))]:
 src=json.load(open(D.parent/'closure'/f'candidate{word[0]}.json'));sn=src['suffix_names'];names=src['names'];edges=list(itertools.combinations(sn,2));hits=map(json.loads,(D/f'grandchildren{word[0]}.jsonl').read_text().splitlines());z=next(z for z in hits if tuple(z[k]for k in ['a','b','W','s'])==sel);ports=[v for i,v in enumerate(sn)if z['W']>>i&1];changes=[edges[z['a']],edges[z['b']]]+[(39,v)for v in sorted(set(src['in_ports'])^set(ports))];rows=src['rows'][:]
 for a,b in changes:
  i,j=names.index(a),names.index(b);rows[i]^=1<<j;rows[j]^=1<<i
 out=dict(source_H=src['target_H'],target_H=z['target'],names=names,rows=rows,root=39,prefix=src['prefix'],suffix_names=sn,in_ports=ports,all_changed_pairs=changes,suffix_flips=changes[:2],h=z['h'],J=z['J'],F=z['F'],root_E2=z['h']+8*z['F'],root_defect=z['J'],scope='Finite candidate with every inherited name preserved; native guard audit pending')
 (D/f'candidate{word}.json').write_text(json.dumps(out,indent=2)+'\n');print(word,changes,ports,out['target_H'])
