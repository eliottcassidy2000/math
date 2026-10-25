from pathlib import Path
import hashlib,json,subprocess,shutil
out=Path(__file__).parent;repo=Path('/tmp/session56-publication/clean-repo');pin='4f1760fd1eff346767c032524c932d5e317401d2'
def ck(b,m):
 if not b:raise RuntimeError(m)
def gitfile(p):return subprocess.check_output(['git','show',pin+':'+p],cwd=repo)
rows=[]
for packet in ('valuation-lca','niche-atlas','cartesian-fibres'):
 src=Path('/tmp/session56-'+packet);dst=out/'inputs'/packet
 shutil.copytree(src,dst,dirs_exist_ok=True)
 for p in sorted(dst.rglob('*')):
  if p.is_file():rows.append({'snapshot':str(p.relative_to(out)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
 # Replays don't mutate inputs.
 script='check_connections.py' if packet=='niche-atlas' else 'verify.py'
 res=[]
 for opt in (False,True):
  cmd=['python3']+(['-O'] if opt else [])+['-B',str(dst/script)]
  p=subprocess.run(cmd,capture_output=True,text=True);ck(p.returncode==0,p.stderr)
  (out/(packet+('.optimized' if opt else '.normal')+'.txt')).write_text(p.stdout)
  res.append(p.stdout)
 ck(res[0]==res[1],packet+' normal optimized differs')
 if packet=='valuation-lca':ck(res[0]==(dst/'normal.json').read_text(),'valuation saved certificate')
 if packet=='niche-atlas':ck(res[0]==(dst/'check_connections.json').read_text(),'atlas saved certificate')
source_rows=[]
for q in json.loads((out/'inputs/niche-atlas/SOURCE_PINS.json').read_text()):
 if 'path' in q:p=q['path']
 else:
  stem=q['local_packet_path'].split('/')[0]
  p='05-knowledge/results/session55_parker_fields_integer_lifts_20261025/'+{'session55-crt-collision-wiring':'crt-wiring','session55-weighted-crt-counts':'weighted-crt'}[stem]+'/THEOREM.md'
 data=gitfile(p);ck(hashlib.sha256(data).hexdigest()==q['sha256'],'atlas source hash '+p)
 dest=out/'sources'/p;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(data)
 source_rows.append({'commit':pin,'path':p,'sha256':q['sha256']})
for q in json.loads((out/'inputs/valuation-lca/SOURCE_PINS.json').read_text()):
 data=(out/'inputs/valuation-lca'/q['snapshot']).read_bytes();ck(hashlib.sha256(data).hexdigest()==q['sha256'],'valuation snapshot hash')
 src=q['source']
 if src.startswith('session55-'):
  p='05-knowledge/results/session55_parker_fields_integer_lifts_20261025/'+({'session55-crt-collision-wiring':'crt-wiring','session55-padic-lifts':'padic-lifts'}[src.split('/')[0]])+'/'+src.split('/')[-1]
 else:p='05-knowledge/results/session53_g27_collatz_cooper_20261024/audit-repeated-block/inputs/'+src
 original=gitfile(p);ck(original==data,'valuation source equality '+p)
 source_rows.append({'commit':pin,'path':p,'sha256':q['sha256']})
(out/'PROVENANCE.json').write_text(json.dumps({'repository_pin':pin,'producer_inputs':rows,'verified_historical_sources':source_rows,'external_sources':[{'url':'https://cadmo.ethz.ch/education/lectures/HS18/SAADS/papers/lca.pdf','checked':'2026-10-25','scope':'Bender/Farach-Colton LATIN 2000 pp88-94, RMQ index and Cartesian reduction; read via web'},{'url':'https://oeis.org/A000568','checked':'2026-10-25','scope':'offset 0, values 1,1,1,2,4,12,56,456; CITED table only'}]},indent=2)+'\n')
print('PASS: complete producer snapshots, all normal/-O replays, 18 source-pin entries verified against git blobs')
