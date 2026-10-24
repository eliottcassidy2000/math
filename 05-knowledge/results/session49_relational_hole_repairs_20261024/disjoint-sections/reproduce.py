#!/usr/bin/env python3
"""Actual manifest-only relocation, no workspace runtime dependencies."""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,sys,tempfile
from verify_manifest import verify
HERE=Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);a=ap.parse_args()
    verify();manifest=json.loads((HERE/'MANIFEST.json').read_text())
    with tempfile.TemporaryDirectory(prefix='session49-disjoint-relocated-') as d:
        dest=Path(d)/'packet';dest.mkdir()
        for row in manifest['files']:
            p=dest/row['path'];p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(HERE/row['path'],p)
        shutil.copyfile(HERE/'MANIFEST.json',dest/'MANIFEST.json')
        result_verify=subprocess.run([sys.executable,'-B','verify_manifest.py'],cwd=dest,check=True,capture_output=True,text=True)
        outputs=[]
        for flags in [[],['-O']]:
            z=subprocess.run([sys.executable,*flags,'-B','check.py'],cwd=dest,check=True,capture_output=True,text=True)
            outputs.append(json.loads(z.stdout))
        expected=json.loads((HERE/'normal.json').read_text())
        if outputs[0]!=outputs[1] or outputs[0]!=expected:raise RuntimeError('relocated mathematical receipt mismatch')
        receipt={'status':'PASS','actual_relocation':True,'temporary_directory':str(dest),'copied_payload':manifest['files'],'manifest_sha256_at_replay':hashlib.sha256((HERE/'MANIFEST.json').read_bytes()).hexdigest(),'manifest_verification':json.loads(result_verify.stdout),'normal_optimized_mathematical_equivalence':True,'normal_receipt_sha256':hashlib.sha256((HERE/'normal.json').read_bytes()).hexdigest(),'checks_per_interpreter':expected['checks'],'commands':[['python','-B','check.py'],['python','-O','-B','check.py']],'temporary_directory_removed_after_success':True,'scope':'Copied only the manifest-listed payload and its manifest; no network or imports from prior sessions. If this receipt is subsequently included in a final manifest, the replayed payload hashes remain authoritative.'}
    Path(a.output).write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
if __name__=='__main__':main()
