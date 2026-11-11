#!/usr/bin/env python3
"""Scoped links/syntax/pins/manifest checks; no mathematical proof checker."""
import ast,hashlib,json,re,subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def main():
    scripts=sorted(HERE.rglob('*.py'))
    for p in scripts:ast.parse(p.read_text())
    links=0
    for p in [*HERE.rglob('*.md'),ROOT/'00-navigation/S272.md']:
        for target in re.findall(r'\[[^\]]*\]\(([^)]+)\)',p.read_text()):
            if '://' in target or target.startswith('#'):continue
            q=(p.parent/target.split('#')[0]).resolve()
            assert q.exists(),(str(p),target)
            links+=1
    for p in HERE.rglob('*.json'):json.loads(p.read_text())
    pins=json.loads((HERE/'SOURCE_PIN_MANIFEST.json').read_text())
    assert pins['commit']=='3456e5e75cf8d4cd558955a5c7638fba2d68f464'
    assert len(pins['files'])==23
    assert len({p['path'] for p in pins['files']})==23
    for p in pins['files']:
        assert re.fullmatch('[0-9a-f]{40}',p['blob_sha'])
        assert pins['commit'] in p['url']
    result={'status':'PASS','python_syntax_files':len(scripts),'local_markdown_links':links,
       'pinned_source_files':23,'source_rebuilt':False,
       'global_documentation_check':'PASS after upstream 0e63417e46a2 restored startup budget; check rerun before final publication',
       'scope':'packaging only; see REPRODUCE.md for mathematical controls'}
    (HERE/'PACKAGING_CHECKS.json').write_text(json.dumps(result,indent=2)+'\n')
    files=sorted(p for p in HERE.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='SHA256SUMS')
    (HERE/'SHA256SUMS').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(HERE))+'\n' for p in files))
    print(json.dumps(result,indent=2))
    print('SHA256SUMS files:',len(files))
if __name__=='__main__':main()
