#!/usr/bin/env python3
"""Verify the closed research packet and replay all exact suites in two modes."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib,json,os,re,subprocess,sys
from pathlib import Path
P=Path(__file__).resolve().parent
ROOT=P.parents[2]
def need(ok,msg):
 if not ok: raise RuntimeError(msg)
def sha(b): return hashlib.sha256(b).hexdigest()
def local(s):
 q=P/s
 need(not q.is_symlink() and q.resolve().is_relative_to(P),'unsafe packet path '+s)
 return q
def blob(rev,path):
 r=subprocess.run(['git','show',rev+':'+path],cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
  env={**os.environ,'GIT_NO_LAZY_FETCH':'1'},timeout=60)
 need(r.returncode==0,'source unavailable '+rev+':'+path+' '+r.stderr.decode(errors='replace'))
 return r.stdout
def run(job):
 s,opt=job;cmd=[sys.executable,'-I','-B']+(['-O'] if opt else [])+[str(local(s['program']))]
 r=subprocess.run(cmd,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=600)
 name='optimized' if opt else 'normal'
 need(r.returncode==0 and not r.stderr,'suite failure '+s['program']+' '+r.stderr.decode(errors='replace'))
 need(r.stdout==local(s[name]).read_bytes(),'receipt mismatch '+s['program']+' '+name)
 d=json.loads(r.stdout);need(d['status']=='PASS' and d[s['count_key']]==s['gates'],'receipt gate/status mismatch')
 return {'program':s['program'],'mode':name,'gates':s['gates']}
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--integrity-only',action='store_true');ap.add_argument('--live',action='store_true');args=ap.parse_args()
 m=json.loads(local('MANIFEST.json').read_text())
 actual=set()
 for f in P.rglob('*'):
  need(not f.is_symlink(),'symlink '+str(f))
  if f.is_file():actual.add(f.relative_to(P).as_posix())
 need(actual==set(m['files'])|{'MANIFEST.json'},'closed file universe mismatch')
 for f,h in m['files'].items():need(sha(local(f).read_bytes())==h,'file hash '+f)
 for b in m['review_bindings']:need(sha(local(b['path']).read_bytes())==b['sha256'],'review binding '+b['path'])
 for b in m['source_pins']:need(sha(blob(b['revision'],b['path']))==b['sha256'],'historical source '+b['path'])
 candidates={}
 for x in m['promotions']:
  text=local(x['path']).read_text();need(sha(text.encode())==x['promoted_sha256'],'promoted hash '+x['path'])
  need(text.count(x['new'])==1,'nonunique promotion '+x['path'])
  old=text.replace(x['new'],x['old'],1);need(sha(old.encode())==x['candidate_sha256'],'promotion inverse '+x['path']);candidates[x['path']]=old
 # Additional proof amendments preserve their documented earlier candidate, not old code receipts.
 for x in m['proof_amendments']:
  text=candidates[x['path']];need(text.count(x['new'])==1,'amendment match '+x['path'])
  text=text.replace(x['new'],x['old'],1);need(sha(text.encode())==x['before_sha256'],'amendment inverse '+x['path'])
 patch=json.loads(local('NAVIGATION_PATCH.json').read_text());before=blob(patch['source_revision'],patch['path'])
 need(sha(before)==patch['before_sha256'],'navigation before hash')
 after=before+patch['append'].encode();need(sha(after)==patch['after_sha256'],'navigation append replay')
 if args.live:need((ROOT/patch['path']).read_bytes()==after,'live navigation changed')
 # Ordinary local Markdown links: verify paths, avoiding external URLs and headings.
 link_count=0
 for f in P.rglob('*.md'):
  for target in re.findall(r'\]\(([^)]+)\)',f.read_text()):
   if target.startswith(('http:','https:','#','mailto:')):continue
   path=target.split('#',1)[0]
   if path:
    need((f.parent/path).exists(),'broken local link '+str(f)+': '+target);link_count+=1
 for s in m['suites']:need(local(s['normal']).read_bytes()==local(s['optimized']).read_bytes(),'mode mismatch '+s['program'])
 runs=[]
 if not args.integrity_only:
  with ThreadPoolExecutor(max_workers=3) as pool:
   runs=list(pool.map(run,[(s,o) for s in m['suites'] for o in (False,True)]))
 print(json.dumps({'status':'PASS','packet_files':len(actual),'source_pins':len(m['source_pins']),
   'promotion_inverses':len(m['promotions']),'proof_amendment_inverses':len(m['proof_amendments']),
   'local_links':link_count,'gates_single_run_per_suite':sum(s['gates'] for s in m['suites']),
   'live_navigation_checked':args.live,'runs':runs},indent=2,sort_keys=True))
if __name__=='__main__':main()
