#!/usr/bin/env python3
"""Replay this independent audit against an explicit frozen producer packet."""
from pathlib import Path
from hashlib import sha256
import argparse,json,subprocess,sys,tempfile
from itertools import combinations
HERE=Path(__file__).resolve().parent
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def verify(folder):
 m=json.loads((folder/'MANIFEST.json').read_text())
 files={str(p.relative_to(folder)) for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p!=folder/'MANIFEST.json'}
 need(files==set(m),'closed manifest '+str(folder))
 for name,row in m.items():
  b=(folder/name).read_bytes();need(len(b)==row['bytes'] and sha256(b).hexdigest()==row['sha256'],'manifest '+name)
def main():
 p=argparse.ArgumentParser();p.add_argument('packet',type=Path);p.add_argument('--verify-only',action='store_true');a=p.parse_args();packet=a.packet.resolve()
 verify(HERE);verify(packet)
 snap=json.loads((HERE/'PACKET_SNAPSHOT.json').read_text());need(sha256((packet/'MANIFEST.json').read_bytes()).hexdigest()==snap['manifest_sha256'],'frozen packet snapshot')
 pins=json.loads((packet/'SOURCEPINS.json').read_text());auth=json.loads((HERE/'source_authentication.json').read_text())
 for r in auth['verified_files']:
  q=pins[r['input']];need(q['upstream_commit']==auth['upstream_commit'] and q['upstream_path']==r['upstream_path'] and q['sha256']==r['sha256'],'source authentication')
 for path,r in auth['nested_source_manifests'].items():
  name=path.split('/')[-2];data=(HERE/'source-manifests'/name/'MANIFEST.json').read_bytes();need(sha256(data).hexdigest()==r['sha256'],'nested source manifest')
  m=json.loads(data)
  for rel in r['matched_pinned_files']:
   found=[v for v in pins.values() if v['upstream_path']==path.rsplit('/',1)[0]+'/'+rel]
   need(len(found)==1 and m[rel]==found[0]['sha256'],'nested source member '+rel)
 host=json.loads((packet/'inputs/host.json').read_text());frames={int(s):set(v) for s,v in host['source_frames'].items()};shared=[]
 for x,y in combinations(frames,2):
  inter=frames[x]&frames[y];expected={3,4} if {x,y}<={368,776,3139} else set();need(inter==expected,'actual helper/core frame intersections')
  if inter:shared.append({'sources':[x,y],'owners':sorted(inter)})
 need(len(set.union(*frames.values()))==host['whole_source_host_order']==652,'exact host union')
 inventory={'status':'PASS','frame_count':len(frames),'pairwise_intersections_checked':len(frames)*(len(frames)-1)//2,'host_union_order':652,'shared_frame_pairs':shared,'all_other_intersections_empty':True}
 need(inventory==json.loads((HERE/'host_inventory.json').read_text()),'host inventory receipt')
 if a.verify_only:print('PASS complete producer/audit closures, nested manifests, authenticated sources and all3321 frame intersections');return
 cmd=[sys.executable]+(['-O'] if sys.flags.optimize else [])+['-B']
 subprocess.run(cmd+[str(packet/'reproduce.py')],check=True,cwd=packet)
 with tempfile.TemporaryDirectory(prefix='relational-independent-replay-') as d:
  out=Path(d)/'receipt.json';subprocess.run(cmd+[str(HERE/'independent_check.py'),str(packet),'--out',str(out)],check=True,cwd=d)
  need(out.read_bytes()==(HERE/'independent-normal.json').read_bytes(),'independent frozen receipt')
 need((HERE/'independent-normal.json').read_bytes()==(HERE/'independent-optimized.json').read_bytes(),'normal optimized receipts identical')
 print('PASS independent audit closed reproduction')
if __name__=='__main__':main()
