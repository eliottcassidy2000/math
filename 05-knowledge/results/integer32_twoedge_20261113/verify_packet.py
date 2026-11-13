#!/usr/bin/env python3
"""Verify frozen artifacts, reversible promotions, source pins and receipts.
Mathematical replay commands and finite universes are in each lane's note.
"""
from hashlib import sha256,sha1
from pathlib import Path
import json,re,subprocess

P=Path(__file__).resolve().parent
ROOT=P.parents[2]
def need(ok,message):
    if not ok:raise RuntimeError(message)
def digest(b):return sha256(b).hexdigest()
def read(name):return json.loads((P/name).read_text())
def actual():
    return {f.relative_to(P).as_posix() for f in P.rglob('*') if f.is_file()
            and '__pycache__' not in f.parts and f.name!='MANIFEST.json'}

def main():
    manifest=read('MANIFEST.json')
    need(actual()==set(manifest['files']),'closed file set differs')
    for name,row in manifest['files'].items():
        f=P/name;b=f.read_bytes()
        need(not f.is_symlink(),'symlink '+name)
        need(len(b)==row['bytes'] and digest(b)==row['sha256'],'manifest '+name)
    historical={}
    for name,row in read('promotions.json')['files'].items():
        s=(P/name).read_text()
        need(digest(s.encode())==row['after_sha256'],'promotion current '+name)
        for edit in reversed(row['edits']):
            need(s.count(edit['new'])==edit['count'],'unambiguous inverse '+name)
            s=s.replace(edit['new'],edit['old'])
        need(digest(s.encode())==row['before_sha256'],'reviewed candidate '+name)
        historical[name]=row['before_sha256']
        for review,h in row['review_files'].items():
            need(digest((P/review).read_bytes())==h,'review binding '+review)
    def bound(name,h):
        need(digest((P/name).read_bytes())==h or historical.get(name)==h,'stale binding '+name)
    for f in P.rglob('*SHA256SUMS'):
        for line in f.read_text().splitlines():
            if not line.strip():continue
            h,name=line.split('  ',1)
            bound((f.parent/name).relative_to(P).as_posix(),h)
    pins=read('SOURCE_PINS.json')['files']
    for name,row in pins.items():
        # Use the recorded object, so later edits to a current source do not
        # rewrite this historical input. All these blobs were fetched at merge.
        b=subprocess.check_output(['git','cat-file','blob',row['git_blob']],cwd=ROOT)
        need(digest(b)==row['sha256'] and len(b)==row['bytes'],'source bytes '+name)
        oid=sha1(('blob '+str(len(b))+'\0').encode()+b).hexdigest()
        need(oid==row['git_blob'],'Git object binding '+name)
    source_digests={row['sha256'] for row in pins.values()}
    total=0
    for row in manifest['receipt_pairs']:
        a=(P/row['normal']).read_bytes();b=(P/row['optimized']).read_bytes()
        need(a==b,'normal/-O mismatch '+row['normal'])
        d=json.loads(a);bound(row['producer'],d[row['producer_field']])
        total+=d.get('checks',0)
        if 'input_sha256' in d:need(d['input_sha256'] in source_digests,'missing immutable input')
        for name,h in d.get('source_pins',{}).items():
            if (P/name).is_file():bound(name,h)
            else:need(name in pins and pins[name]['sha256']==h,'receipt source '+name)
        for name,h in d.get('reviewed_sha256',{}).items():
            if (P/name).is_file():bound(name,h)
            else:
                candidates=[f.relative_to(P).as_posix() for f in P.rglob(Path(name).name) if f.is_file()]
                need(any(digest((P/c).read_bytes())==h or historical.get(c)==h for c in candidates),
                     'reviewed candidate binding '+name)
    for f in P.rglob('*.md'):
        for target in re.findall(r'\]\(([^\n)]+)\)',f.read_text()):
            if '://' in target or target.startswith('#'):continue
            name=target.split('#')[0]
            if not name:continue
            resolved=(f.parent/name).resolve()
            need(resolved.is_relative_to(ROOT),'escaped local link '+target)
            need(resolved.exists(),'broken link '+str(resolved))
    catalog=read('census/catalog.json')
    rows=catalog.get('rows',catalog.get('entries',[]))
    need(len(rows)==33,'33 integer routes')
    print('PASS: closed manifest; reviewed reversible promotions; '+str(len(pins))+
          ' historical source pins; '+str(len(manifest['receipt_pairs']))+
          ' matching receipt pairs; '+str(total)+' recorded explicit gates; local links.')

if __name__=='__main__':main()
