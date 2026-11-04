#!/usr/bin/env python3
from pathlib import Path,PurePosixPath
from hashlib import sha256,sha1
import json,gzip,subprocess
D=Path(__file__).resolve().parent;BASE=D.parent;PUB=BASE/'session190-publish';RESULT='05-knowledge/results/session190_constructive_links_20261104'
def need(x,msg):
    if not x:raise RuntimeError(msg)
def dg(data):return {'bytes':len(data),'sha256':sha256(data).hexdigest(),'git_blob_sha1':sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()}
def check(data,record,label):
    for k,v in dg(data).items():need(record[k]==v,label+' '+k)
def safe(path):
    p=PurePosixPath(path);need(not p.is_absolute()and'..'not in p.parts and str(p)==path,'safe path');return path

def read_package(folder):
    path=PUB/folder/'package.json';obj=json.loads(path.read_text());entries={}
    for item in obj['entries']:
        rel=safe(item['path']);need(rel not in entries,'distinct object path');raw=Path(item['object_file']).read_bytes();check(raw,item,rel);need(len(raw)<=48*1024,'object bound');entries[rel]=raw
    mpath=RESULT+'/ARCHIVES_'+obj['checkpoint']+'.json';manifest=json.loads(entries[mpath]);need(manifest['schema']==1,'manifest schema')
    originals={}
    for f in manifest['files']:
        rel=safe(f['path']);need(rel not in originals,'distinct original path')
        if f['storage']=='direct':raw=entries[rel]
        else:
            need(f['storage']=='gzip_chunks','known encoding');parts=[]
            for part in f['parts']:
                chunk=entries[safe(part['path'])];check(chunk,part,'chunk');parts.append(chunk)
            zipped=b''.join(parts);check(zipped,f['compressed'],'concatenated compressed bytes');raw=gzip.decompress(zipped)
        check(raw,f,rel);originals[rel]=raw
    need(len(originals)==len(obj['source_files']),'same originals count')
    for src in obj['source_files']:
        check(originals[src['path']],src,'package source metadata')
    return obj,entries,manifest,originals

def main():
    cp3,entries,manifest,originals=read_package('CP03-candidate')
    expected={}
    for source,filename,dest in [('session190-spectrum','RELEASE_FILES.json','spectrum/spectrum'),('session190-spectrum-audit','PUBLICATION_FILES.json','spectrum/spectrum_audit'),('session190-interface','FROZEN_MANIFEST.json','interface/session190-interface'),('session190-planar-parity-audit','FROZEN_MANIFEST.json','interface/session190-planar-parity-audit')]:
        p=BASE/source;m=json.loads((p/filename).read_text());names=m.get('publication_files',list(m['files'])+[filename])
        for name in names:expected[RESULT+'/'+dest+'/'+name]=(p/name).resolve()
    expected['00-navigation/S190.md']=(PUB/'S190-CP03.md').resolve();expected[RESULT+'/REPRODUCTION_CP03.md']=(PUB/'REPRODUCTION-CP03.md').resolve()
    need(set(originals)==set(expected),'CP3 exact complete approved path set')
    mapobj=json.loads((PUB/'cp03-map.json').read_text());need({f['path']:Path(f['source']).resolve()for f in mapobj['files']}==expected,'source map exact allowlist')
    for f in cp3['source_files']:need(Path(f['source']).resolve()==expected[f['path']],'expected source identity')
    for rel,source in expected.items():need(originals[rel]==source.read_bytes(),'lossless source bytes '+rel)
    root=D/'composite_repo';need(not root.exists(),'fresh composite root');root.mkdir()
    latest={};allentries={};duplicates=[];manifest_paths=[];package_pins={}
    for folder in ['CP01-final','CP02-final-v2','CP03-candidate']:
        obj,objects,m,raws=read_package(folder);package_pins[folder]=sha256((PUB/folder/'package.json').read_bytes()).hexdigest()
        for rel,raw in objects.items():allentries[rel]=raw
        for rel,raw in raws.items():
            if rel in latest:duplicates.append(rel)
            latest[rel]=raw
        manifest_paths.append(root/RESULT/('ARCHIVES_'+obj['checkpoint']+'.json'))
    need(all(p in ['00-navigation/S190.md',RESULT+'/RESTORE_ARCHIVES.py',RESULT+'/RECEIPTS.md',RESULT+'/paper_patch/REPRODUCTION.md']for p in duplicates),'only administrative records superseded')
    for rel,raw in allentries.items():
        dest=root/rel;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(raw)
    restorer=root/RESULT/'RESTORE_ARCHIVES.py'
    before={rel:sha256(raw).hexdigest()for rel,raw in allentries.items()}
    command=['python3',str(restorer),*[str(p)for p in manifest_paths],'--repo-root',str(root)]
    first=subprocess.check_output(command,text=True);(D/'composite_restore.stdout').write_text(first)
    for rel,raw in latest.items():need((root/rel).read_bytes()==raw,'independent latest original verification '+rel)
    second=subprocess.check_output(command+['--verify-only'],text=True);(D/'composite_verify.stdout').write_text(second)
    optimized=subprocess.check_output(['python3','-O']+command[1:]+['--verify-only'],text=True);need(second==optimized,'restorer normal optimized');(D/'composite_verify.optimized.stdout').write_text(optimized)
    # Idempotent second ordinary restoration creates no additional files and returns stable latest bytes.
    third=subprocess.check_output(command,text=True);need('restored 0 files'in third,'idempotent restoration')
    for rel,h in before.items():need(sha256((root/rel).read_bytes()).hexdigest()==h,'stored object was not modified')
    for rel,raw in originals.items():need((root/rel).read_bytes()==expected[rel].read_bytes(),'CP3 restored exact source '+rel)
    result=dict(status='PASS',checkpoint='CP03',original_files=len(originals),stored_objects=len(entries),source_bytes=sum(map(len,originals.values())),stored_bytes=sum(map(len,entries.values())),approved_source_count=len(expected),package_pins=package_pins,composite_current_original_count=len(latest),superseded_paths=duplicates,normal_optimized_restorer_equal=True,source_pins={p:dg(raw)for p,raw in originals.items()},entry_pins={p:dg(raw)for p,raw in entries.items()})
    (D/'PACKAGE_REVIEW_RECEIPT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items()if k not in ('source_pins','entry_pins')},indent=2))
if __name__=='__main__':main()
