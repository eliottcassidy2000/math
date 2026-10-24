#!/usr/bin/env python3
"""Verify complete packet closure; optionally replay only listed copied files."""
from hashlib import sha256
from pathlib import Path
import argparse
import json
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent


def need(condition,message):
    if not condition:
        raise RuntimeError(message)


def validate():
    manifest = json.loads((HERE/'MANIFEST.json').read_text())
    actual = {p.relative_to(HERE).as_posix() for p in HERE.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='MANIFEST.json'}
    need(actual==set(manifest), 'manifest must list every packet file except itself')
    for name,row in manifest.items():
        relative = Path(name)
        need(not relative.is_absolute() and '..' not in relative.parts, 'safe manifest relative path')
        data = (HERE/name).read_bytes()
        need(len(data)==row['bytes'] and sha256(data).hexdigest()==row['sha256'], 'manifest integrity '+name)
    return manifest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--relocate',action='store_true')
    args = parser.parse_args()
    manifest = validate()
    if args.relocate:
        with tempfile.TemporaryDirectory(prefix='relational-hole-manifest-only-') as directory:
            root = Path(directory)/'packet'
            root.mkdir()
            for name in list(manifest)+['MANIFEST.json']:
                target = root/name
                target.parent.mkdir(parents=True,exist_ok=True)
                shutil.copyfile(HERE/name,target)
            copied = {p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file()}
            need(copied==set(manifest)|{'MANIFEST.json'}, 'relocated copy contains only manifest closure')
            for optimize in (False,True):
                command = [sys.executable]+(['-O'] if optimize else [])+['-B',str(root/'reproduce.py')]
                subprocess.run(command,cwd=root,check=True)
            print('PASS actual manifest-only relocated normal/optimized replay;',len(copied),'copied files')
    else:
        command = [sys.executable]+(['-O'] if sys.flags.optimize else [])+['-B',str(HERE/'check.py')]
        subprocess.run(command,cwd=HERE,check=True)
        print('PASS complete manifest closure;',len(manifest),'listed files')


if __name__=='__main__':
    main()
