"""Load the unchanged Session87 implementation, checking its exact source hash.

Public use finds the repository root among this file's ancestors. Staging may
pass --baseline-root PATH to any checker; this loader consumes that option.
No source snapshot is vendored and no files in the baseline tree are written.
"""
from hashlib import sha256
from pathlib import Path
import sys
import types

SOURCE_REL = Path('05-knowledge/results/session87_four_vertices_enumeration_20261028/check_generation.py')
SOURCE_SHA256 = 'f5d3da961f34dc1a32eaa0c1b69fad9bda983154578ea0187c5fd6b2e2d10fae'
SOURCE_PIN = '7a0bab9bfa7bb87c5a10c7ba239592d90cfafd51'


def locate():
    roots=[]; argv=sys.argv; i=1
    while i<len(argv):
        if argv[i]=='--baseline-root':
            if i+1>=len(argv): raise RuntimeError('--baseline-root needs a directory')
            roots.append(Path(argv[i+1]).resolve()); del argv[i:i+2]
        elif argv[i].startswith('--baseline-root='):
            roots.append(Path(argv[i].split('=',1)[1]).resolve()); del argv[i]
        else: i+=1
    if len(roots)>1: raise RuntimeError('provide --baseline-root at most once')
    if not roots: roots=list(Path(__file__).resolve().parents)
    for root in roots:
        path=root/SOURCE_REL
        if path.is_file(): return path
    raise RuntimeError('Session87 baseline not found; run within the repository or supply --baseline-root PATH')


def load_generation():
    path=locate(); data=path.read_bytes()
    if sha256(data).hexdigest()!=SOURCE_SHA256:
        raise RuntimeError('Session87 baseline hash mismatch: '+str(path))
    module=types.ModuleType('session88_verified_session87_generation')
    module.__file__=str(path)
    exec(compile(data,str(path),'exec'),module.__dict__)
    return module


generation=load_generation()
