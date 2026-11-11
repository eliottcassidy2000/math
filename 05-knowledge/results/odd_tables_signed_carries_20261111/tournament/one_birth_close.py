"""Reproduce the complete one-newborn-cut test at the named H5991 source."""
from pathlib import Path
from tempfile import TemporaryDirectory
import subprocess,json

HERE=Path(__file__).resolve().parent
rows=[162,1020,1017,1009,993,960,897,768,513,1]
root=0
for _ in range(5):
    v=len(rows)
    rows[root]|=1<<v
    rows.append(((1<<v)-1)^(1<<root))
    root=v
for u,v in ((0,10),(0,11)):
    rows[u]^=1<<v
    rows[v]^=1<<u
data=' '.join(map(str,[len(rows),root,5993]+rows))+'\n'
with TemporaryDirectory(prefix='odd-tables-unit-close-') as td:
    exe=Path(td)/'one_birth_close'
    subprocess.run(['clang++','-O3','-std=c++17',str(HERE/'one_birth_close.cpp'),'-o',str(exe)],check=True)
    result=json.loads(subprocess.check_output([str(exe)],input=data,text=True))
result.update(source_rows=rows,root=root,scope='Every one-newborn cut on the named15-vertex H5991 source; root-indegree retention is recorded, not imposed on the search.')
(HERE/'one_birth_close_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
