#!/usr/bin/env python3
"""Replay a receipt-writing checker in an isolated copy; return exact receipt bytes."""
from pathlib import Path
import argparse,shutil,subprocess,sys,tempfile
ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--directory',required=True);p.add_argument('--script',required=True);a=p.parse_args()
rel=Path(a.directory)
if rel.is_absolute() or '..' in rel.parts:raise ValueError('unsafe lane directory')
with tempfile.TemporaryDirectory(prefix='s145-receipt-') as temp:
 d=Path(temp)/'lane';shutil.copytree(ROOT/rel,d)
 command=[sys.executable,'-B']+(['-O'] if sys.flags.optimize else [])+[str(d/a.script),'verification.replay.json']
 result=subprocess.run(command,cwd=d,capture_output=True,check=True)
 sys.stdout.buffer.write((d/'verification.replay.json').read_bytes())
