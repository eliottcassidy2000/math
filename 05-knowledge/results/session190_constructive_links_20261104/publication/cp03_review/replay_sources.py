#!/usr/bin/env python3
from pathlib import Path
from hashlib import sha256
import json,shutil,subprocess
D=Path(__file__).resolve().parent
SOURCE=D.parent
ROOT=D/'source_replay'
LOG=D/'replay_logs'
LOG.mkdir(exist_ok=True)
def need(x,msg):
    if not x:raise RuntimeError(msg)
def digest(p):return sha256(p.read_bytes()).hexdigest()
SPECS=[('session190-spectrum','RELEASE_FILES.json','spectrum/spectrum'),('session190-spectrum-audit','PUBLICATION_FILES.json','spectrum/spectrum_audit'),('session190-interface','FROZEN_MANIFEST.json','interface/session190-interface'),('session190-planar-parity-audit','FROZEN_MANIFEST.json','interface/session190-planar-parity-audit')]

def main():
    need(not ROOT.exists(),'fresh isolated replay destination required')
    originals={};packets=[]
    for name,manifest,dest in SPECS:
        src=SOURCE/name;obj=json.loads((src/manifest).read_text())
        if 'publication_files'in obj:
            names=obj['publication_files'];pins={r['path']:r['sha256']for r in obj['files']}
        else:names=list(obj['files'])+[manifest];pins=obj['files']
        target=ROOT/dest;target.mkdir(parents=True)
        for n in names:
            if n in pins:need(digest(src/n)==pins[n],'source pin '+name+'/'+n)
            shutil.copyfile(src/n,target/n)
            rel=str((target/n).relative_to(ROOT));originals[rel]=digest(src/n)
            need(digest(target/n)==originals[rel],'copied byte identity')
        packets.append(dict(source=name,manifest=manifest,manifest_sha256=digest(src/manifest),destination=dest,files=len(names)))
    checks=[]
    def run_pair(rel,script,label,reference=None):
        for optimized in [False,True]:
            suffix='optimized'if optimized else'normal';cmd=['python3']+(['-O']if optimized else[])+[str(ROOT/rel/script)]
            completed=subprocess.run(cmd,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True)
            (LOG/(label+'.'+suffix+'.stdout')).write_bytes(completed.stdout)
            (LOG/(label+'.'+suffix+'.stderr')).write_bytes(completed.stderr)
            if reference is not None:need(completed.stdout==(ROOT/rel/reference).read_bytes(),'frozen stdout '+label)
        a=(LOG/(label+'.normal.stdout')).read_bytes();b=(LOG/(label+'.optimized.stdout')).read_bytes();need(a==b,'normal optimized '+label)
        checks.append(dict(label=label,script=str(Path(rel)/script),normal_optimized_equal=True,stdout_sha256=sha256(a).hexdigest(),reference=reference))
        print('PASS '+label,flush=True)
    run_pair('spectrum/spectrum','verify.py','producer-spectrum','normal.stdout')
    run_pair('spectrum/spectrum','run_census.py','producer-census','census.stdout')
    for script,label,reference in [('independent_audit.py','independent-spectrum','audit.stdout'),('wrapper_audit.py','independent-wrapper','wrapper_audit.stdout'),('census_audit.py','independent-census','census_audit.stdout'),('bind_audit.py','independent-binding','bind_audit.stdout')]:run_pair('spectrum/spectrum_audit',script,label,reference)
    run_pair('interface/session190-interface','verify.py','producer-interface','RESULTS.json')
    run_pair('interface/session190-planar-parity-audit','audit.py','independent-interface','RESULTS.json')
    for rel,h in originals.items():need(digest(ROOT/rel)==h,'retained frozen bytes unchanged after replay '+rel)
    result=dict(status='PASS',scope='Source-bound isolated CP3 mathematical packet replays; publication bytes and immutable remote commit remain separately checked.',packets=packets,original_file_count=len(originals),source_pins=originals,checks=checks)
    (D/'SOURCE_REPLAY_RECEIPT.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS all isolated source-bound CP3 replays and final original-byte comparison',flush=True)
if __name__=='__main__':main()
