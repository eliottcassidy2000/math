#!/usr/bin/env python3
"""Verify the closed packet and optional standard-library replay, without Z3."""
import argparse,hashlib,json,subprocess,sys
from pathlib import Path

def digest(data):return hashlib.sha256(data.replace(b'\r\n',b'\n')).hexdigest()
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--replay',action='store_true');args=parser.parse_args()
    here=Path(__file__).resolve().parent;repo=here.parents[2]
    manifest=json.loads((here/'MANIFEST.json').read_text())
    actual={p.name for p in here.iterdir() if p.is_file() and p.name!='MANIFEST.json'}
    need(actual==set(manifest['packet']),'closed packet file inventory')
    for name,pin in manifest['packet'].items():need(digest((here/name).read_bytes())==pin,('packet pin',name))
    for path,item in manifest['sources'].items():
        local=repo/path
        data=local.read_bytes() if local.exists() else subprocess.check_output(['git','show',f"{item['revision']}:{path}"],cwd=repo)
        need(digest(data)==item['sha256'],('source pin',path))
    for path,item in manifest.get('repaired_incoming',{}).get('paths',{}).items():
        local=repo/path
        data=local.read_bytes() if local.exists() else subprocess.check_output(['git','show','HEAD:'+path],cwd=repo)
        need(digest(data)==item['after_sha256'],('repaired incoming source',path))
    replayed=0
    for stem in ('geometry','ideals','port_codes','golden_tilt','repairs','reverse','incoming'):
        normal=(here/f'{stem}.normal.json').read_bytes();optimized=(here/f'{stem}.optimized.json').read_bytes()
        need(normal==optimized,('normal/optimized receipt pair',stem))
        if args.replay:
            for flags in ([],['-O']):
                output=subprocess.check_output([sys.executable,*flags,'-B',str(here/f'check_{stem}.py')],cwd=repo)
                need(output==normal,('fresh replay',stem,flags));replayed+=1
    print(json.dumps({'status':'PASS','packet_files':len(actual),'source_pins':len(manifest['sources']),
        'normal_optimized_pairs':7,'fresh_checker_runs':replayed,
        'repaired_incoming_sources':len(manifest.get('repaired_incoming',{}).get('paths',{})),
        'scope':'closed inventory and exact controls; mathematical proofs internally reviewed, not formally verified'},sort_keys=True))
if __name__=='__main__':main()
