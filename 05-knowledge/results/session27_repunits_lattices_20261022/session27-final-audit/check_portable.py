#!/usr/bin/env python3
"""Reconstruct omitted repair snapshots in a temporary directory, then run the original source-bound integration check. Requires a math Git checkout containing the manifest's pinned preimage blobs. Usage: python check_portable.py /path/to/math"""
from pathlib import Path
import subprocess,sys,json,tempfile,shutil
HERE=Path(__file__).resolve().parent
PACKET=HERE.parent
REPO=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(subprocess.check_output(['git','-C',str(PACKET),'rev-parse','--show-toplevel'],text=True).strip())
with tempfile.TemporaryDirectory(prefix='session27-integration-') as td:
 shared=Path(td)
 for p in PACKET.iterdir():
  if p.is_dir() and p.name.startswith('session27-') and p.name!='session27-corrections':
   (shared/p.name).symlink_to(p,target_is_directory=True)
 corrections=shared/'session27-corrections';corrections.mkdir()
 for p in (PACKET/'session27-corrections').iterdir():
  if p.is_file():shutil.copyfile(p,corrections/p.name)
 manifest=json.loads((corrections/'manifest.json').read_text())
 for f in manifest['files']:
  target=corrections/'originals'/f['path'];target.parent.mkdir(parents=True,exist_ok=True)
  target.write_bytes(subprocess.check_output(['git','-C',str(REPO),'cat-file','blob',f['git_blob_oid']]))
 shutil.copytree(corrections/'originals',corrections/'proposed')
 subprocess.run(['git','apply',str(corrections/'corrections.patch')],cwd=corrections/'proposed',check=True)
 subprocess.run([sys.executable,str(HERE/'check_integration.py'),str(PACKET),str(shared)],check=True)
