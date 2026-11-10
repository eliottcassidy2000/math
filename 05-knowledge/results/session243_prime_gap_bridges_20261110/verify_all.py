#!/usr/bin/env python3
"""Repository-local exact Session243 replay, with intentional failing controls."""
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]
BOOTSTRAP="import runpy,sys; from pathlib import Path; p=Path(sys.argv[1]).resolve(); sys.path.insert(0,str(p.parent)); sys.argv=sys.argv[1:]; runpy.run_path(str(p),run_name='__main__')"

def require(ok,message):
    if not ok:raise RuntimeError(message)

def run(entry):
    relative,frozen,optimized,failure=entry
    command=[sys.executable,'-I','-B']+(['-O'] if optimized else [])+['-c',BOOTSTRAP,str(HERE/relative)]
    if failure:command.append('--fail')
    p=subprocess.run(command,capture_output=True,check=False)
    if failure:
        require(p.returncode!=0 and not p.stdout,'negative control did not fail: '+relative)
        require(b'deliberate' in p.stderr or b'intentional' in p.stderr,'unrelated negative failure: '+relative)
    else:
        require(p.returncode==0 and not p.stderr,'replay failed: '+relative+' '+p.stderr.decode())
        payload=json.loads(p.stdout)
        require(str(payload.get('status','')).startswith('PASS'),'missing PASS: '+relative)
        require(p.stdout==(HERE/frozen).read_bytes(),'frozen output mismatch: '+relative)
    return {'script':relative,'optimized':optimized,'negative_control':failure,'exit_code':p.returncode,'stdout_sha256':sha256(p.stdout).hexdigest()}

def main():
    arithmetic=subprocess.run([sys.executable,'-I','-B',str(HERE/'arithmetic_audit/replay.py')],capture_output=True,check=False)
    require(arithmetic.returncode==0 and not arithmetic.stderr,'arithmetic replay failed '+arithmetic.stderr.decode())
    ar=json.loads(arithmetic.stdout);require(ar['status']=='PASS','arithmetic audit failed')
    require(arithmetic.stdout==(HERE/'arithmetic_audit/REPLAY.json').read_bytes(),'arithmetic replay receipt stale')
    programs=[('constructions/verify.py','constructions/CHECK.json'),
              ('constructions/verify_gaps.py','constructions/GAPS.json'),
              ('construction_audit/verify_independent.py','construction_audit/CHECK.json'),
              ('construction_audit/verify_gaps_independent.py','construction_audit/GAPS_CHECK.json'),
              ('construction_audit/verify_arithmetic_independent.py','construction_audit/ARITHMETIC_CHECK.json')]
    jobs=[(p,o,opt,fail) for p,o in programs for opt in [False,True] for fail in [False,True]]
    with ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(run,jobs))
    for program,_ in programs:
        normal=next(r for r in results if r['script']==program and not r['optimized'] and not r['negative_control'])
        optimized=next(r for r in results if r['script']==program and r['optimized'] and not r['negative_control'])
        require(normal['stdout_sha256']==optimized['stdout_sha256'],'normal/optimized mismatch')
    pins=HERE/'construction_audit/SOURCES.json'
    source_pins=json.loads(pins.read_text())
    checked=0
    for entry in source_pins['inherited_sources']:
        raw=subprocess.run(['git','show',entry['commit']+':'+entry['path']],cwd=REPO,capture_output=True,check=True).stdout
        require(sha256(raw).hexdigest()==entry['sha256'],'inherited source pin mismatch: '+entry['path']);checked+=1
    for category in ['canon_files_reviewed','producer_files_reviewed']:
        for entry in source_pins[category]:
            require(sha256((REPO/entry['path']).read_bytes()).hexdigest()==entry['sha256'],'reviewed source pin mismatch: '+entry['path']);checked+=1
    pin_result=sha256(pins.read_bytes()).hexdigest()
    print(json.dumps({'status':'PASS','positive_replays':ar['positive_replays']+10,'negative_controls':ar['negative_replays']+10,'arithmetic_replay_sha256':sha256(arithmetic.stdout).hexdigest(),'physical_source_pins_verified':checked,'physical_source_pin_file_sha256':pin_result,'physical_replays':results,'scope':'Exact finite replays supplement the stated proofs; no global coverage or prime infinitude claim.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
