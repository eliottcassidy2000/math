import json,re,pathlib,hashlib
R=pathlib.Path(__file__).resolve().parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
expected={line.split(': ',1)[0].split('/')[-1]:line.split(': ',1)[1] for line in (R/'window_scan.out').read_text().splitlines()}
summary={}
for tag,want in [('normal',83),('n12',45)]:
 rs=[json.loads(l) for l in (R/f'replay-{tag}.jsonl').read_text().splitlines()];need(len(rs)==want,tag+': incomplete')
 need(len({r['job'] for r in rs})==want,tag+': repeats')
 for r in rs:
  need(r['returncode']==0 and not r['stderr'],tag+': executable error')
  need('HIT' not in r['stdout'],tag+': target hit')
  if tag=='normal':need(r['stdout'].strip().split(': ',1)[1]==expected[r['job']],r['job']+': original mismatch')
 s=''.join(r['stdout'] for r in rs)
 totals={x:sum(map(int,re.findall(r' '+x+r'=(\d+)',s))) for x in ['strong','autcalls','checks','failures','hits']}
 need(totals['hits']==totals['failures']==0,tag+': failure')
 totals.update(jobs=want,total_job_seconds=sum(r['seconds'] for r in rs),source_lines_matched=want if tag=='normal' else None,sha256=hashlib.sha256((R/f'replay-{tag}.jsonl').read_bytes()).hexdigest())
 summary[tag]=totals
summary['source_commit']='87d61f67673a6b10cbafff1bf9e1377f3c41e56c'
summary['aut_control']=(R/'test-aut.out').read_text().strip()
summary['input_audit']=json.loads((R/'audit-jobs-normal.json').read_text())
print(json.dumps(summary,sort_keys=True,indent=2))
