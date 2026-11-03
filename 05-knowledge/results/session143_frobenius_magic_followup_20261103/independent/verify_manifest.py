#!/usr/bin/env python3
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent


def need(c,m):
    if not c:raise ValueError(m)


def main():
    manifest=json.loads((ROOT/'manifest.json').read_text())
    for path,wanted in manifest['files'].items():
        need(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==wanted,'manifest mismatch: '+path)
    pins=json.loads((ROOT/'source_pins.json').read_text())
    for item in pins['files']:
        data=(ROOT/item['local_path']).read_bytes()
        need(len(data)==item['bytes'] and hashlib.sha256(data).hexdigest()==item['sha256'],'source pin mismatch')
    need((ROOT/'result.json').read_bytes()==(ROOT/'result.optimized.json').read_bytes(),'optimized result mismatch')
    need((ROOT/'negative_controls.json').read_bytes()==(ROOT/'negative_controls.optimized.json').read_bytes(),'optimized negative-control mismatch')
    print(json.dumps({'status':'PASS','manifest_files':len(manifest['files']),'repository_source_pins':len(pins['files']),'normal_optimized':'byte-identical'},sort_keys=True))

if __name__=='__main__':main()
