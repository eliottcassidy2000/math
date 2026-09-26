"""Reproduce four exact lanes and pin their sources, outputs and dependencies."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[2]
PREFIX='forest_20260926_'
LANES=('excursions','bernoulli','tilings','gilbreath')


def require(ok,label):
    if not ok:
        raise RuntimeError(label)


def normalized(data):
    return data.decode('utf-8').replace('\r\n','\n').encode('utf-8')


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--refresh',action='store_true')
    args=parser.parse_args()
    pins={}
    for lane in LANES:
        source=Path('04-computation/experiments')/(PREFIX+lane+'.py')
        output=Path('05-knowledge/results')/(PREFIX+lane+'.out')
        runs=[]
        for optimized in (False,True):
            cmd=[sys.executable,'-X','utf8','-B']
            if optimized:
                cmd.append('-O')
            result=subprocess.run(cmd+[str(ROOT/source)],cwd=ROOT,capture_output=True,timeout=180)
            require(result.returncode==0,lane+': '+result.stderr.decode('utf-8',errors='replace'))
            runs.append(normalized(result.stdout))
        require(runs[0]==runs[1],lane+' differs under -O')
        target=ROOT/output
        if target.exists():
            require(normalized(target.read_bytes())==runs[0],lane+' retained output mismatch')
        else:
            require(args.refresh,lane+' output missing; review and refresh')
        if args.refresh:
            target.write_bytes(runs[0])
        pins[source.as_posix()]=digest(normalized((ROOT/source).read_bytes()))
        pins[output.as_posix()]=digest(runs[0])
        print('PASS '+lane+': normal = optimized = retained')
    for rel in [
        '01-canon/theorems/THM-4507-finite-valuation-collatz-potential-obstruction.md',
        '01-canon/theorems/THM-4483-forced-charges-two-place-ranks.md',
        '05-knowledge/results/collatz_blueprint_20260921_affine.md',
        '05-knowledge/results/forest_20260926_tilings.svg',
    ]:
        pins[rel]=digest(normalized((ROOT/rel).read_bytes()))
    png='05-knowledge/results/forest_20260926_tilings.png'
    pins[png]=digest((ROOT/png).read_bytes())
    data={'hash_basis':'UTF-8 LF for text; raw bytes for PNG','sha256':pins}
    manifest=ROOT/'05-knowledge/results/forest_20260926_manifest.json'
    if args.refresh:
        manifest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    else:
        require(json.loads(manifest.read_text(encoding='utf-8'))==data,'hash manifest changed')
    print('PASS '+str(len(pins))+' source/output/dependency/figure pins')


if __name__=='__main__':
    main()
