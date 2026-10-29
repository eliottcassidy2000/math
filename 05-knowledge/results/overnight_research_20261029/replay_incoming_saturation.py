"""Replay a frozen incoming certificate from exact Git blobs, preserving hashes.

The original verifier hashes source bytes. A Windows CRLF checkout correctly
fails that raw-byte guard. Extract the six original blobs to an owned temporary
directory rather than altering the verifier or weakening its source check.
"""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[3]
PIN = '1128ee8635775e21277560e0902cf615cc4f2015'
BASE = '05-knowledge/results/session97_overnight_arithmetic_carriers_20261029/halving'
FILES = ['halving_certificate.json', 'phase2/verify_extension.py',
         'phase2/extension_3567.json', 'phase2/extension_3568.json',
         'phase2/inverse_transport_3567.json', 'phase2/inverse_transport_3568.json']


def main():
    folder = Path(tempfile.mkdtemp(prefix='math-saturation-git-replay-')).resolve()
    if not folder.is_relative_to(Path(tempfile.gettempdir()).resolve()):
        raise RuntimeError('unexpected temporary location')
    hashes = {}
    for relative in FILES:
        data = subprocess.check_output(['git','show',f'{PIN}:{BASE}/{relative}'],cwd=ROOT)
        target = folder/relative
        target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes(data)
        hashes[relative] = hashlib.sha256(data).hexdigest()
    verifier = folder/'phase2/verify_extension.py'
    normal = subprocess.check_output([sys.executable,'-B',str(verifier)],cwd=folder)
    optimized = subprocess.check_output([sys.executable,'-B','-O',str(verifier)],cwd=folder)
    if normal != optimized:
        raise RuntimeError('incoming normal/optimized outputs differ')
    out = {'incoming_commit':PIN,'git_blob_sha256':hashes,
           'normal_optimized_identical':True,'verifier_result':json.loads(normal),
           'checkout_boundary':'raw CRLF bytes differ; original Git-blob guard retained unchanged'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__ == '__main__':
    main()
