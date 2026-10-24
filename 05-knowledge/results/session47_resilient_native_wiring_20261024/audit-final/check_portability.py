#!/usr/bin/env python3
"""Verify final audit against a copied, manifest-only publication layout.
This writes temporary files only, never the source packets. No Git or network.
"""
from pathlib import Path
import argparse,json,shutil,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--root-dir',type=Path,default=HERE.parent);args=parser.parse_args();root=args.root_dir.resolve()
aliases={'session47-helper-resilience/final-v3':'helper-resilience-v3','session47-audit-helper-resilience/third-audit':'audit-helper-resilience-v3','session47-audit-suspension-obstruction':'audit-five-to-seven','session47-audit-fresh570':'fresh570-audit'}
def locate(packet):
 original=root/packet
 if original.exists():return original
 return root/aliases.get(packet,packet.removeprefix('session47-'))
def copy(source,destination):
 destination.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,destination)
reports=json.loads((HERE/'manifests.normal.json').read_text())['packets'];tests=[];copied=0
with tempfile.TemporaryDirectory(prefix='session47-final-manifest-only-') as td:
 bundle=Path(td)/'renamed-bundle';bundle.mkdir()
 for row in reports:
  packet=row['packet'];src=locate(packet);dest=bundle/aliases.get(packet,packet.removeprefix('session47-'))
  for name in row['paths']+['MANIFEST.json']:
   copy(src/name,dest/name);copied+=1
 audit=bundle/'audit-final'
 own=json.loads((HERE/'MANIFEST.json').read_text())['files']
 for name in list(own)+['MANIFEST.json']:copy(HERE/name,audit/name);copied+=1
 for name in ['SYNTHESIS.md','README.md','SHARED_OWNER_REPLAY.json']:
  src=root/'session47-root'/'final'/name
  if not src.exists():src=root/name
  copy(src,bundle/name);copied+=1
 for script,receipt in [('check.py','normal.json'),('check_manifests.py','manifests.normal.json')]:
  for optimized in [False,True]:
   command=[sys.executable,'-I']+(['-O'] if optimized else [])+['-B',str(audit/script),'--root-dir',str(bundle)]
   result=subprocess.run(command,cwd='/',capture_output=True,check=True)
   if result.stdout!=(HERE/receipt).read_bytes():raise RuntimeError('Manifest-only receipt differs: '+script)
   tests.append({'script':script,'optimized':optimized,'status':'PASS byte-identical receipt'})
 for optimized in [False,True]:
  command=[sys.executable,'-I']+(['-O'] if optimized else [])+['-B',str(audit/'verify_manifest.py'),'--root-dir',str(bundle)]
  result=subprocess.run(command,cwd='/',capture_output=True,check=True)
  if json.loads(result.stdout)['status']!='PASS':raise RuntimeError('Copied bundle integrity failed')
  tests.append({'script':'verify_manifest.py','optimized':optimized,'status':'PASS copied integrity and reviewed source hashes'})
print(json.dumps({'status':'PASS','layout':'Manifest-listed files only, copied into renamed publication siblings; isolated Python -I from /','copied_files_including_manifests':copied,'tests':tests},indent=2,sort_keys=True))
