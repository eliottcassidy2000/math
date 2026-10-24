from pathlib import Path
import json
H=Path(__file__).resolve().parent
routes=json.loads((H/'packet.json').read_text())['routes']+json.loads((H/'packet_eighteen.json').read_text())['routes']+json.loads((H/'packet_eighteen37.json').read_text())['routes'];byname={r['name']:r for r in routes};charges=[]
for file in ['charges.json','charges_eighteen.json','charges_eighteen37.json']:charges.extend(json.loads((H/file).read_text()))
keys={tuple(sorted((int(k),v) for k,v in r['charge'].items())):r['name'] for r in charges};newcharges=[];variants=[]
for base,suffix,path in [('buffer_1394','min_work',[65280,116276,86071,255]),('eighteen_buffer_1394','min_work',[261632,437825,274217,279403,511]),('eighteen37_buffer_2071','min_clock',[261632,453144,428849,360765,360700,511])]:
 r=dict(byname[base]);r['name']=base+'_'+suffix;r['path']=path;r['step_sizes']=[(a^b).bit_count() for a,b in zip(path,path[1:])];names=[]
 for j,(a,b) in enumerate(zip(path,path[1:])):
  charge={s:2*((b>>i&1)-(a>>i&1)) for i,s in enumerate(r['sources']) if (a^b)>>i&1};key=tuple(sorted(charge.items()))
  if key not in keys:
   name=r['name']+f'_step_{j+1}';newcharges.append(dict(name=name,charge=charge));keys[key]=name
  names.append(keys[key])
 r['step_image_names']=names;variants.append(r)
(H/'packet_variants.json').write_text(json.dumps(dict(routes=variants),indent=2,sort_keys=True)+'\n');(H/'charges_variants.json').write_text(json.dumps(newcharges,indent=2,sort_keys=True)+'\n');print('variants',len(variants),'new image charges',len(newcharges))
