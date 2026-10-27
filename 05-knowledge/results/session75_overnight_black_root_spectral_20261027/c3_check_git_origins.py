#!/usr/bin/env python3
"""Optional targeted historical Git-origin replay; never a corpus search."""
import argparse,gzip,json,subprocess
from pathlib import Path
from c3_provenance import check
from c3_transport import read_regular,need,sha

def main(packet,repo):
 check(packet)
 index=json.loads(read_regular(packet,'TRANSPORT.checkpoint3.json'));rows={r['logical_path']:r for r in index['files']}
 origins=json.loads(read_regular(packet,'SOURCE_ORIGINS.checkpoint3.json'))
 for row in origins['sources']:
  r=rows[row['packet_logical_path']];b=read_regular(packet,r['stored_path']);b=gzip.decompress(b) if r['encoding']=='gzip' else b
  command=['git','-C',str(repo),'show',row['commit']+':'+row['repository_path']]
  result=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True)
  need(result.stdout==b,'historical Git source byte mismatch: '+row['repository_path'])
 return dict(status='PASS',targeted_historical_git_snapshots=len(origins['sources']),scope='exact declared source bytes at historical commits; no whole-corpus scan or mathematical audit')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--packet',type=Path,default=Path(__file__).resolve().parent);p.add_argument('--repo',type=Path,required=True);a=p.parse_args();print(json.dumps(main(a.packet,a.repo),indent=2,sort_keys=True))
