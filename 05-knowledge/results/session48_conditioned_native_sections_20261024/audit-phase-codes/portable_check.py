#!/usr/bin/env python3
"""Copy only manifest-listed bytes and test from an unrelated cwd."""
from pathlib import Path
from hashlib import sha256
import argparse,json,shutil,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent

def run(script,flags,cwd,extra=()):
    p=subprocess.run([sys.executable,*flags,str(script),*map(str,extra)],cwd=cwd,capture_output=True,text=True,check=True)
    return {'stdout':p.stdout,'stderr':p.stderr,'exit_code':p.returncode}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    manifest=json.loads((HERE/'MANIFEST.json').read_text())
    with tempfile.TemporaryDirectory(prefix='independent-phase-audit-') as raw:
        temp=Path(raw);clone=temp/'packet';clone.mkdir();cwd=temp/'unrelated-cwd';cwd.mkdir()
        for name in ['MANIFEST.json',*manifest['files']]:
            dest=clone/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(HERE/name,dest)
        results={'schema':'manifest-only-portable-phase-audit-v1','copied_listed_files':len(manifest['files'])}
        results['initial_manifest']=run(clone/'verify_manifest.py',['-I','-B'],cwd)
        for label,flags in [('normal',['-I','-B']),('optimized',['-I','-O','-B'])]:
            out=temp/(label+'.json')
            results[label]=run(clone/'independent_check.py',flags,cwd,['--output',out])
            if out.read_bytes()!=(clone/(label+'.json')).read_bytes():raise RuntimeError('Portable independent receipt mismatch: '+label)
            results[label]['receipt_sha256']=sha256(out.read_bytes()).hexdigest()
        if results['normal']!=results['optimized']:raise RuntimeError('Portable independent modes differ')
        results['producer_replay']=run(clone/'replay_producer.py',['-I','-B'],cwd)
        results['final_manifest']=run(clone/'verify_manifest.py',['-I','-B'],cwd)
    if args.write:(HERE/'PORTABILITY.json').write_text(json.dumps(results,indent=2,sort_keys=True)+'\n')
    elif json.loads((HERE/'PORTABILITY.json').read_text())!=results:raise RuntimeError('Portable result differs from frozen receipt')
    print('PASS manifest-only isolated audit: independent normal/optimized, producer replay, byte-exact closure')

if __name__=='__main__':main()
