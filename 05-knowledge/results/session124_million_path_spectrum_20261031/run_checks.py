#!/usr/bin/env python3
"""Reproduce the complete finite-spectrum packet using external working files."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parent

def require(value, message):
    if not value: raise RuntimeError(message)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--work-dir',type=Path,required=True)
    args=p.parse_args(); work=args.work_dir.resolve(); work.mkdir(parents=True,exist_ok=True)
    logs=work/'receipts'; logs.mkdir(exist_ok=True)
    started=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
    records=[]
    def run(name, command):
        t=time.monotonic()
        env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
        result=subprocess.run([str(x) for x in command],cwd=work,env=env,capture_output=True,text=True)
        (logs/(name+'.out')).write_text(result.stdout)
        if result.stderr: (logs/(name+'.stderr')).write_text(result.stderr)
        require(result.returncode==0,(name,result.returncode,result.stderr[-2000:]))
        records.append({'check':name,'returncode':result.returncode,'seconds':round(time.monotonic()-t,3)})
        print(name,'PASS',flush=True)
        return result.stdout
    def same(a,b):require(a==b,'Normal/optimized receipts differ')
    python=sys.executable
    for mode in ('normal','optimized'):
        flags=[] if mode=='normal' else ['-O']
        output=run('seeded_words.'+mode,[python,*flags,ROOT/'constructors/seeded_words.py',5,1000000,'--output-roots',work/('seed_roots.'+mode+'.json')])
        require('states 10621974 values 499980' in output,'Primary support size differs')
        require((work/('seed_roots.'+mode+'.json')).read_bytes()==(ROOT/'constructors/seed_roots.json').read_bytes(),'Seed witness data differs')
    same((logs/'seeded_words.normal.out').read_bytes(),(logs/'seeded_words.optimized.out').read_bytes())
    for filename in ('audit_seeded_words','check_target_witnesses','check_block_targets','check_mod3_and_products','check_padic_filtration'):
        a=run(filename+'.normal',[python,ROOT/'constructors'/(filename+'.py')])
        b=run(filename+'.optimized',[python,'-O',ROOT/'constructors'/(filename+'.py')]);same(a,b)
    for opt in ('O0','O2'):
        run('compile.independent.'+opt,['g++','-'+opt,'-std=c++17','-Wall','-Wextra','-pedantic',ROOT/'audit/independent.cpp','-o',work/('independent.'+opt)])
        run('independent.'+opt,[work/('independent.'+opt),work/opt])
    same((logs/'independent.O0.out').read_bytes(),(logs/'independent.O2.out').read_bytes())
    same((work/'O0/seeds.tsv').read_bytes(),(work/'O2/seeds.tsv').read_bytes())
    same((work/'O0/support.bin').read_bytes(),(work/'O2/support.bin').read_bytes())
    same((ROOT/'audit/seeds.tsv').read_bytes(),(work/'O2/seeds.tsv').read_bytes())
    a=run('independent_check.normal',[python,ROOT/'audit/check.py','--work-dir',work/'O2'])
    b=run('independent_check.optimized',[python,'-O',ROOT/'audit/check.py','--work-dir',work/'O0']);same(a,b)
    adjacency=sorted((ROOT/'audit').glob('adjacency_*.txt'))
    for opt in ('O0','O2'):
        run('compile.adjacency.'+opt,['g++','-'+opt,'-std=c++17','-Wall','-Wextra','-pedantic',ROOT/'audit/adjacency_dp.cpp','-o',work/('adjacency.'+opt)])
        run('adjacency.'+opt,[work/('adjacency.'+opt),*adjacency])
    same((logs/'adjacency.O0.out').read_bytes(),(logs/'adjacency.O2.out').read_bytes())
    for filename in ('mod3','mutation'):
        extra=[work/'adjacency.O2'] if filename=='mutation' else []
        a=run(filename+'.normal',[python,ROOT/'audit'/(filename+'.py'),*extra])
        b=run(filename+'.optimized',[python,'-O',ROOT/'audit'/(filename+'.py'),*extra]);same(a,b)
    # Reconstruct a boundary, a product, a block repair, and the top allowed value.
    for target in (1,49,34869,103061,999999):
        for mode in ('normal','optimized'):
            flags=[] if mode=='normal' else ['-O']
            run(f'construct.{target}.{mode}',[python,*flags,ROOT/'construct_witness.py',target,'--adjacency',work/f'constructed.{target}.{mode}.txt','--recipe',work/f'constructed.{target}.{mode}.json'])
        same((work/f'constructed.{target}.normal.txt').read_bytes(),(work/f'constructed.{target}.optimized.txt').read_bytes())
        same((work/f'constructed.{target}.normal.json').read_bytes(),(work/f'constructed.{target}.optimized.json').read_bytes())
    samples=sorted(work.glob('constructed.*.normal.txt'))
    for opt in ('O0','O2'):
        run('constructed_adjacency.'+opt,[work/('adjacency.'+opt),*samples])
    same((logs/'constructed_adjacency.O0.out').read_bytes(),(logs/'constructed_adjacency.O2.out').read_bytes())
    a=run('hostile_cli.normal',[python,ROOT/'audit/hostile_cli.py',work/'adjacency.O2'])
    b=run('hostile_cli.optimized',[python,'-O',ROOT/'audit/hostile_cli.py',work/'adjacency.O2']);same(a,b)
    generated={p.relative_to(work).as_posix():{'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in sorted(work.rglob('*')) if p.is_file() and p.suffix in ('.bin','.tsv','.json','.txt')}
    (logs/'verification.json').write_text(json.dumps({'status':'PASS','started_utc':started,'finished_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'python':platform.python_version(),'compiler':subprocess.check_output(['g++','--version'],text=True).splitlines()[0],'checks':records,'generated_data':generated,'source_hashes':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(ROOT.rglob('*')) if p.is_file() and p.suffix in ('.py','.cpp','.jsonl') and 'publication' not in p.relative_to(ROOT).parts and 'receipts' not in p.relative_to(ROOT).parts},'scope':'Every positive odd H<=999999 except 7,21 is attained; not an infinite completeness proof; global exclusions not reproved.'},indent=2)+'\n')
    print('ALL CHECKS PASS',flush=True)
if __name__=='__main__':main()
