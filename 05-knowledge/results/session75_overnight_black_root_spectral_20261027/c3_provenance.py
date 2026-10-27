#!/usr/bin/env python3
"""Portable binding checks. Historical Git origins were verified at assembly;
this command authenticates their exact retained snapshots, not live GitHub.
"""
import argparse,gzip,hashlib,json
from pathlib import Path
from c3_transport import verify,read_regular,need,sha,safe

def check(packet):
 packet=Path(packet).resolve();index=json.loads(read_regular(packet,'TRANSPORT.checkpoint3.json'));verify(packet,index)
 rows={r['logical_path']:r for r in index['files']};cache={}
 def data(n):
  if n not in cache:
   if n in rows:
    r=rows[n];b=read_regular(packet,r['stored_path']);cache[n]=gzip.decompress(b) if r['encoding']=='gzip' else b
   else:cache[n]=read_regular(packet,n)
  return cache[n]
 p=json.loads(read_regular(packet,'PROVENANCE.checkpoint3.json'));need(p['schema']=='session75-c3-bindings-v1','binding schema')
 origins=json.loads(data(p['source_origins_file']));need(origins['status']=='PASS','historical origin status')
 for row in origins['sources']:
  need(row['commit'] in ('96e35642110104cc547516a63f37eebcf930f707','87d61f67673a6b10cbafff1bf9e1377f3c41e56c'),'unapproved historical source pin')
  safe(row['repository_path'])
  b=data(row['packet_logical_path']);need(len(b)==row['bytes'] and sha(b)==row['sha256'],'origin snapshot mismatch')
  need(hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()==row['git_blob_sha1'],'Git blob identity mismatch')
 for row in p['ledger_bindings']:
  x=json.loads(data(row['ledger']))
  for k in row['pointer']:x=x[k]
  need(x==row['sha256'],'ledger field mismatch')
  for n in row['targets']:need(sha(data(n))==x,'ledger target hash mismatch')
 for row in p['exact_audit_copy_bindings']:
  need(data(row['path'])==data(row['equal_to']) and sha(data(row['path']))==row['sha256'],'audit copy mismatch')
 for row in p['unresolved_historical_hash_fields']:
  x=json.loads(data(row['ledger']))
  for k in row['pointer']:x=x[k]
  need(x==row['sha256'],'historical provenance field mismatch')
 return dict(status='PASS',historical_git_origin_snapshots=len(origins['sources']),ledger_bindings=len(p['ledger_bindings']),exact_audit_copy_bindings=len(p['exact_audit_copy_bindings']),historical_hash_fields_outside_new_transport=len(p['unresolved_historical_hash_fields']),live_git_origin_check='NOT_RUN by portable verifier; exact targeted assembly receipt retained')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--packet',type=Path,default=Path(__file__).resolve().parent);a=ap.parse_args();print(json.dumps(check(a.packet),indent=2,sort_keys=True))
