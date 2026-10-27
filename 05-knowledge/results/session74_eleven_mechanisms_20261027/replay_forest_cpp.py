#!/usr/bin/env python3
"""Replay the native C++ forest audit without modifying its source or receipts."""
import argparse,hashlib,json,subprocess,tempfile
from pathlib import Path

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--lane',type=Path,default=Path(__file__).resolve().parent/'session74-score-audit');a=ap.parse_args();lane=a.lane.resolve();source=lane/'forest_audit.cpp'
    expected=(lane/'forest_audit.out').read_bytes();need(expected==(lane/'forest_audit.sanitized.out').read_bytes(),'native C++ receipts differ')
    source_before=source.read_bytes();rows=[]
    with tempfile.TemporaryDirectory(prefix='session74-cpp-replay-') as td:
        for name,flags in [('optimized',['-O2']),('ubsan',['-O1','-g','-fsanitize=undefined','-fno-sanitize-recover=all'])]:
            binary=Path(td)/name
            c=subprocess.run(['g++','-std=c++17',*flags,'-Wall','-Wextra',str(source),'-o',str(binary)],capture_output=True,timeout=120)
            need(c.returncode==0,'C++ compile failed: '+c.stderr.decode(errors='replace'))
            q=subprocess.run([str(binary),'9'],capture_output=True,timeout=600)
            need(q.returncode==0,'C++ execution failed: '+q.stderr.decode(errors='replace'))
            need(not q.stderr,'C++ unexpected stderr: '+q.stderr.decode(errors='replace'))
            need(q.stdout==expected,'C++ output differs from native exact receipt')
            rows.append({'mode':name,'status':'PASS','stdout_sha256':hashlib.sha256(q.stdout).hexdigest()})
    need(source.read_bytes()==source_before,'C++ source changed')
    print(json.dumps({'status':'PASS','scope':'Exact native ordered-forest enumeration n=2..9; no universal theorem inferred','source_sha256':hashlib.sha256(source_before).hexdigest(),'runs':rows},indent=2,sort_keys=True))
if __name__=='__main__':main()
