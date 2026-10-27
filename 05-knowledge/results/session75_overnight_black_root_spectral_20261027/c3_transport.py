#!/usr/bin/env python3
"""Fail-closed C3 transport verification/materialization; no producer mutation.

Every original manifest entry survives as original bytes or deterministic gzip.
Expanded-byte checks happen before outputs are written. No assertions are used.
"""
import argparse,gzip,hashlib,json,os,stat,tempfile
from pathlib import Path,PurePosixPath

def need(ok,msg):
    if not ok: raise ValueError(msg)
def sha(data): return hashlib.sha256(data).hexdigest()
def safe(name):
    need(isinstance(name,str) and name and '\\' not in name,'invalid path')
    p=PurePosixPath(name)
    need(name == p.as_posix() and not p.is_absolute() and all(x not in ('..','.','') for x in name.split('/')),'unsafe path: '+name)
    return name
def read_regular(root,name):
    p=root/safe(name)
    need(p.is_file() and not p.is_symlink(),'missing/non-regular file: '+name)
    for parent in p.parents:
        if parent==root: break
        need(not parent.is_symlink(),'symlink directory: '+name)
    return p.read_bytes()
def parse_native(data,kind):
    if kind=='json':
        rows=json.loads(data)['files'];need(isinstance(rows,dict),'native files not mapping')
    else:
        rows={}
        for line in data.decode().splitlines():
            if not line.strip():continue
            h,n=line.split(maxsplit=1);n=n.lstrip('*').removeprefix('./')
            need(n not in rows,'duplicate native path');rows[n]=h
    for n,h in rows.items():
        safe(n);need(isinstance(h,str) and len(h)==64 and all(c in '0123456789abcdef' for c in h),'invalid digest')
    return rows

def verify(packet,index,output=None):
    packet=Path(packet).resolve();need(index['schema']=='session75-c3-transport-v1','wrong schema')
    if output is not None:
        output=Path(output).resolve();need(not output.exists(),'materialization target must not exist')
        need(not output.is_relative_to(packet),'materialization must be outside packet')
    seen=set();physical=set();expanded={};modes={}
    for row in index['files']:
        n=safe(row['logical_path']);p=safe(row['stored_path'])
        need(n not in seen and p not in physical,'duplicate logical/stored path')
        seen.add(n);physical.add(p)
        b=read_regular(packet,p)
        need(len(b)==row['stored_bytes'] and sha(b)==row['stored_sha256'],'stored byte/hash mismatch: '+p)
        need(stat.S_IMODE((packet/p).stat().st_mode)==int(row['stored_mode'],8),'stored mode mismatch: '+p)
        need(row['encoding'] in ('identity','gzip'),'unknown encoding')
        if row['encoding']=='gzip':
            with gzip.GzipFile(fileobj=__import__('io').BytesIO(b)) as f:
                data=f.read(row['logical_bytes']+1)
            need(len(data)==row['logical_bytes'],'expanded size mismatch: '+n)
        else:data=b
        need(len(data)==row['logical_bytes'] and sha(data)==row['logical_sha256'],'expanded byte/hash mismatch: '+n)
        mode=int(row['logical_mode'],8);need(mode in (0o644,0o755),'unexpected native mode')
        expanded[n]=data;modes[n]=mode
    for lane in index['lanes']:
        name=safe(lane['name']);mp=name+'/'+safe(lane['manifest'])
        need(mp in expanded and sha(expanded[mp])==lane['manifest_sha256'],'native manifest pin mismatch: '+name)
        rows=parse_native(expanded[mp],lane['manifest_kind'])
        expected={name+'/'+n for n in rows}|{mp}
        actual={n for n in seen if n.startswith(name+'/')}
        need(actual==expected,'native manifest membership mismatch: '+name)
        for n,h in rows.items():need(sha(expanded[name+'/'+n])==h,'native manifest byte mismatch: '+name+'/'+n)
        actual_physical={x.relative_to(packet).as_posix() for x in (packet/name).rglob('*') if x.is_file()}
        need(actual_physical=={r['stored_path'] for r in index['files'] if r['logical_path'].startswith(name+'/')},'unexpected stored lane files: '+name)
    if output is not None:
        output.mkdir(parents=True)
        for n,b in expanded.items():
            p=output/n;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b);p.chmod(modes[n])
    return dict(status='PASS',logical_files=len(seen),logical_bytes=sum(len(b) for b in expanded.values()),stored_bytes=sum(r['stored_bytes'] for r in index['files']),native_manifests=len(index['lanes']),compressed_files=sum(r['encoding']=='gzip' for r in index['files']))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--packet',type=Path,default=Path(__file__).resolve().parent);ap.add_argument('--index',default='TRANSPORT.checkpoint3.json');ap.add_argument('--materialize',type=Path);a=ap.parse_args()
    index=json.loads(read_regular(a.packet,a.index));print(json.dumps(verify(a.packet,index,a.materialize),sort_keys=True,indent=2))
if __name__=='__main__':main()
