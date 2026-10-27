#!/usr/bin/env python3
import argparse,copy,json,shutil,tempfile
from pathlib import Path
from c3_provenance import check

def main(packet):
 results=[]
 with tempfile.TemporaryDirectory(prefix='session75-c3-provenance-corruption-') as tmp:
  p=Path(tmp)/'packet';shutil.copytree(packet,p,ignore=shutil.ignore_patterns('__pycache__'))
  check(p)
  def reject(label,file,mutate,expected):
   b=(p/file).read_bytes();v=json.loads(b);mutate(v);(p/file).write_text(json.dumps(v))
   try:
    try:check(p)
    except (ValueError,KeyError) as e:
     if expected not in str(e):raise ValueError(label+' wrong gate '+str(e))
     results.append(dict(control=label,gate=expected));return
    raise ValueError(label+' unexpectedly accepted')
   finally:(p/file).write_bytes(b)
  f='SOURCE_ORIGINS.checkpoint3.json'
  for key,value,gate in [('sha256','0'*64,'origin snapshot mismatch'),('bytes',0,'origin snapshot mismatch'),('git_blob_sha1','0'*40,'Git blob identity mismatch'),('commit','0'*40,'unapproved historical source pin'),('repository_path','../escape','unsafe path')]:
   reject('source-'+key,f,lambda x,k=key,v=value:x['sources'][0].__setitem__(k,v),gate)
  f='PROVENANCE.checkpoint3.json'
  reject('ledger-hash',f,lambda x:x['ledger_bindings'][0].__setitem__('sha256','0'*64),'ledger field mismatch')
  reject('audited-copy-hash',f,lambda x:x['exact_audit_copy_bindings'][0].__setitem__('sha256','0'*64),'audit copy mismatch')
  reject('historical-raw-hash',f,lambda x:x['unresolved_historical_hash_fields'][0].__setitem__('sha256','0'*64),'historical provenance field mismatch')
  check(p)
 return dict(status='PASS',controls=results,native_packet_unchanged=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--packet',type=Path,default=Path(__file__).resolve().parent);a=ap.parse_args();print(json.dumps(main(a.packet),indent=2,sort_keys=True))
