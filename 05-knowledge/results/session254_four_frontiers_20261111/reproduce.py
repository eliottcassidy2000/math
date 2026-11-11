#!/usr/bin/env python3
"""Replay only Session254's new-claim controls, not inherited conjecture audits."""
from pathlib import Path
import json,subprocess,sys,time
root=Path(__file__).resolve().parent
rows=[]
for lane,script,expected,is_json in [
 ('chan_pak','probe.py','output.json',True),
 ('chan_pak','probe81.py','output81.json',True),
 ('hadwiger_nelson','check.py','output.json',True),
 ('jacobian','probe.py','probe.out',False),
 ('jacobian','depth_semigroup.py','depth_semigroup.out',False),
 ('lonely_runner','reproduce.py','output.json',True),
 ('niche','probe.py','probe.out',False)]:
 start=time.monotonic(); folder=root/lane
 before=(folder/'certificate.json').read_bytes() if lane=='niche' else None
 result=subprocess.run([sys.executable,str(folder/script)],cwd=folder,text=True,capture_output=True,check=True)
 golden=(folder/expected).read_text()
 if is_json: match=json.loads(result.stdout)==json.loads(golden)
 elif lane=='niche':
  def normalize(s):return [' '.join(z.split(' seconds')[0].split()) for z in s.splitlines()]
  match=normalize(result.stdout)==normalize(golden) and before==(folder/'certificate.json').read_bytes()
 else: match=result.stdout==golden
 if not match:raise RuntimeError(lane+' reproduction mismatch')
 rows.append({'lane':lane,'result':'PASS','seconds':round(time.monotonic()-start,2),'matching_artifact':expected})
 print(lane+': PASS',flush=True)
(root/'replay.json').write_text(json.dumps({'scope':'new Session254 controls only','results':rows},indent=2)+'\n')
print('PASS: all seven new-claim engines match their artifacts')
