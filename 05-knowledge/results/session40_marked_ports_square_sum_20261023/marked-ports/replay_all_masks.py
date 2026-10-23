#!/usr/bin/env python3
import argparse,gzip,hashlib,json,subprocess,tempfile
from pathlib import Path
H=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--wiring',type=Path,default=H.parent/'session38-wiring');p.add_argument('--write',action='store_true');a=p.parse_args()
with tempfile.TemporaryDirectory(prefix='session40-endpoint-') as d:
 d=Path(d);names=['classes.tsv.gz','class_of_mask.u16le.gz','edges.tsv.gz','ports.tsv.gz'];paths=[]
 for name in names:
  x=d/name.removesuffix('.gz');x.write_bytes(gzip.decompress((a.wiring/name).read_bytes()));paths.append(str(x))
 subprocess.run(['g++','-O3','-std=c++17',str(H/'replay_all_masks.cpp'),'-o',str(d/'replay')],check=True)
 run=subprocess.run([str(d/'replay'),*paths],check=True,text=True,capture_output=True)
 out=json.loads(run.stdout);out['input_sha256']={name:hashlib.sha256((a.wiring/name).read_bytes()).hexdigest() for name in names}
 out['scope']='Every fixed-path n8 target decoded independently from endpoint refinement plus shared chirality; exact full ordinary-port and merged-edge ledgers; literal within-class partner identity not reconstructed.'
 dest=H/'full_mask_receipt.json'
 if a.write:dest.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
 elif out!=json.loads(dest.read_text()):raise ValueError('frozen full-mask receipt differs')
 print(json.dumps(out,indent=2,sort_keys=True))
