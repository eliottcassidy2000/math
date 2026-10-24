#!/usr/bin/env python3
import json
from pathlib import Path
import independent as a
HERE = Path(__file__).resolve().parent

def main():
    manifest_hash = a.verify_manifest(HERE)
    pins = json.loads((HERE/'SOURCE_PINS.json').read_text())
    a.verify_manifest(HERE/'producer',pins['producer_manifest_sha256'])
    a.verify_manifest(HERE/'nonrigid',pins['nonrigid_manifest_sha256'])
    a.require(pins['nonrigid_manifest_sha256'] == a.NONRIGID_PIN,'hard-coded inherited freeze')
    for left,right in [('normal.json','optimized.json'),('hostile-normal.json','hostile-optimized.json'),('crosscheck-normal.json','crosscheck-optimized.json')]:
        a.require((HERE/left).read_bytes() == (HERE/right).read_bytes(),'normal optimized equality: '+left)
    output = json.loads((HERE/'normal.json').read_text())
    a.require(a.derive_closure(output['closure']['labelled_regular11'],output['nonrigid']) == output['closure'],'closure receipt arithmetic')
    a.require(json.loads((HERE/'hostile-normal.json').read_text())['semantic_rejections'] == 14,'all independent hostiles')
    print(json.dumps(dict(status='PASS',manifest_sha256=manifest_hash,producer_manifest_sha256=pins['producer_manifest_sha256'],nonrigid_manifest_sha256=pins['nonrigid_manifest_sha256']),sort_keys=True))
if __name__ == '__main__':
    main()
