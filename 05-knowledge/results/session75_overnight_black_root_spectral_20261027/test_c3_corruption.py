#!/usr/bin/env python3
"""Mutation controls must reject at the intended gate in Python and Python -O."""
import argparse,copy,json,shutil,tempfile
from pathlib import Path
from c3_transport import verify

def main(packet):
 index=json.loads((packet/'TRANSPORT.checkpoint3.json').read_text());passed=[]
 with tempfile.TemporaryDirectory(prefix='session75-c3-corruption-') as tmp:
  root=Path(tmp)/'packet';root.mkdir()
  for lane in index['lanes']:shutil.copytree(packet/lane['name'],root/lane['name'])
  verify(root,index)
  def reject(label,j,expected):
   try:verify(root,j)
   except (ValueError,OSError,KeyError) as e:
    if expected not in str(e):raise ValueError(label+' wrong failure '+str(e))
    passed.append({'control':label,'gate':expected});return
   raise ValueError(label+' unexpectedly accepted')
  gz=next(r for r in index['files'] if r['encoding']=='gzip');p=root/gz['stored_path'];b=p.read_bytes();v=bytearray(b);v[len(v)//2]^=1;p.write_bytes(v)
  reject('compressed byte flip',index,'stored byte/hash mismatch');p.write_bytes(b)
  j=copy.deepcopy(index);next(r for r in j['files'] if r['encoding']=='gzip')['logical_sha256']='0'*64
  reject('expanded hash mismatch',j,'expanded byte/hash mismatch')
  j=copy.deepcopy(index);next(r for r in j['files'] if r['encoding']=='gzip')['logical_bytes']-=1
  reject('expanded size mismatch',j,'expanded size mismatch')
  j=copy.deepcopy(index);j['files'][0]['logical_path']='../escape';reject('unsafe logical path',j,'unsafe path')
  j=copy.deepcopy(index);j['files'][0]['stored_path']='../escape';reject('unsafe stored path',j,'unsafe path')
  j=copy.deepcopy(index);j['files'].pop();reject('missing logical manifest entry',j,'native manifest membership mismatch')
  j=copy.deepcopy(index);j['files'].append(j['files'][0]);reject('duplicate logical entry',j,'duplicate logical/stored path')
  j=copy.deepcopy(index);j['lanes'][0]['manifest_sha256']='0'*64;reject('false native manifest pin',j,'native manifest pin mismatch')
  p=root/index['files'][0]['stored_path'];mode=p.stat().st_mode&0o777;p.chmod(0o755 if mode==0o644 else 0o644)
  reject('changed stored executable mode',index,'stored mode mismatch');p.chmod(mode)
  extra=root/index['lanes'][0]['name']/'unexpected.data';extra.write_text('unlisted');reject('unexpected native lane member',index,'unexpected stored lane files');extra.unlink()
  verify(root,index)
 return {'status':'PASS','controls':passed,'native_packet_unchanged':True}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--packet',type=Path,default=Path(__file__).resolve().parent);a=ap.parse_args();print(json.dumps(main(a.packet),sort_keys=True,indent=2))
