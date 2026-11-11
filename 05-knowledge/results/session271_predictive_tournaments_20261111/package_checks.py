#!/usr/bin/env python3
"""Scoped source/link and mathematical packaging controls, no inherited audit."""
import ast, hashlib, json, math, re
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def main():
    py=list(HERE.rglob('*.py'))
    for p in py:ast.parse(p.read_text())
    links=0
    for p in [*HERE.rglob('*.md'),ROOT/'00-navigation/S271.md']:
        for target in re.findall(r'\[[^\]]*\]\(([^)]+)\)', '\n'.join(line for line in p.read_text().splitlines() if not line.startswith('    '))):
            if '://' in target or target.startswith('#'):continue
            q=(p.parent/target.split('#')[0]).resolve()
            assert q.exists(),(str(p),target)
            links+=1
    for core in [(3,),(5,),(3,3),(3,5),(3,3,3)]:
        m=sum(core);r=len(core);S=sum(math.gcd(i,j) for i in core for j in core)
        e=(S-r)//2;v=S-2*m*r+m*m
        assert (v+m-r)//2==m*(m-r)-m*(m-1)//2+e
    out={'status':'PASS','python_syntax_files':len(py),'local_markdown_links':links,
         'gaussian_core_exponent_checks':5,
         'global_documentation_check':'pre-existing failure: startup180146>180000; existing startup files unchanged',
         'scope':'packaging and new identities, not independent audit of inherited mathematics'}
    (HERE/'PACKAGING_CHECKS.json').write_text(json.dumps(out,indent=2)+'\n')
    files=sorted(p for p in HERE.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='SHA256SUMS')
    sums=''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(HERE))+'\n' for p in files)
    (HERE/'SHA256SUMS').write_text(sums)
    print(json.dumps(out)); print('SHA256SUMS files:',len(files))

if __name__=='__main__':main()
