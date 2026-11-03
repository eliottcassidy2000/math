#!/usr/bin/env python3
"""Reproduce complete declared censuses; independent algorithms must agree."""
from pathlib import Path
import subprocess

HERE=Path(__file__).resolve().parent


def run(exe,a,b,level,extra=()):
    return subprocess.run([exe,str(HERE),str(a),str(b),str(level),*extra],
                          capture_output=True,text=True,check=True)


if __name__=='__main__':
    for a,b,level,stem in ((2,4,1,'catalog7_2'),(3,3,1,'catalog7_3'),
                           (2,5,2,'catalog8_2'),(3,4,2,'catalog8_3'),
                           (2,6,2,'catalog9_2'),(3,5,2,'catalog9_3'),(4,4,2,'catalog9_4')):
        prod=run('/tmp/boundary-carry-catalog',a,b,level)
        check=run('/tmp/boundary-carry-audit-catalog',a,b,level)
        if prod.stderr.strip()!=check.stdout.strip():raise RuntimeError((a,b,prod.stderr,check.stdout))
        (HERE/(stem+'.jsonl')).write_text(prod.stdout)
        (HERE/(stem+'.stats')).write_text(prod.stderr)
        (HERE/(stem.replace('catalog','audit')+'.stats')).write_text(check.stdout)
        print(stem,prod.stderr.strip(),flush=True)
    for a,b in ((2,4),(3,3)):
        result=run('/tmp/boundary-carry-audit-catalog',a,b,2)
        (HERE/f'audit7targets_{a}.stats').write_text(result.stdout)
        print('lower order',a,b,result.stdout.strip(),flush=True)
    for a,b in ((3,6),(4,5)):
        result=run('/tmp/boundary-carry-catalog',a,b,2,('cover2',))
        (HERE/f'catalog10_{a}.jsonl').write_text(result.stdout)
        (HERE/f'catalog10_{a}.stats').write_text(result.stderr)
        print('partial order ten',a,b,result.stderr.strip(),flush=True)
