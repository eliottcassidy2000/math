#!/usr/bin/env python3
"""Reject malformed or unsaturated exterior context rows."""
from pathlib import Path
from tempfile import TemporaryDirectory
import subprocess
H=Path(__file__).resolve().parent

def need(x,msg):
 if not x:raise RuntimeError(msg)
base=(H/'block4.contexts').read_text().splitlines()
a=list(map(int,base[0].split())); n,k=a[:2]; edges=list(zip(a[2::2],a[3::2]));need(len(edges)==k,'source parsing')
def row(N,es): return f'{N} {len(es)} '+ ' '.join(f'{x} {y}' for x,y in es)
variants={
 'endpoint untouched': row(n,[(x,y) for x,y in edges if x!=0 and y!=0]),
 'old-old':row(n,edges+[(0,1)]),
 'loop':row(n,edges+[(9,9)]),
 'reversed duplicate':row(n,edges+[(edges[0][1],edges[0][0])]),
 'out of bounds':row(n,edges+[(0,n)]),
 'truncated':f'{n} 1 0',
 'trailing':row(n,edges)+' 17',
}
with TemporaryDirectory() as d:
 for name,v in variants.items():
  path=Path(d)/'contexts';path.write_text('\n'.join([v]+base[1:])+'\n')
  p=subprocess.run([str(H/'verify9'),'4',str(path),str(H/'block4.reps')],capture_output=True,text=True)
  need(p.returncode!=0 and 'FAIL:' in p.stderr,f'{name} accepted or crashed: {p.stderr}')
print('rejected:',', '.join(variants))
