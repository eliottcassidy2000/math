"""Reproduce both fixed ordered FRACTRAN lists; no external dependencies."""
import hashlib
import json
from pathlib import Path
from compiler import Compiler, write_artifacts
ROOT = Path(__file__).resolve().parent
for stem, metadata in [('hamilton','PROGRAM.json'),('square_shell','SQUARE_SHELL_PROGRAM.json')]:
    source = ROOT / (stem+'_source.py')
    p = Compiler().compile(source.read_text())
    info = write_artifacts(p, ROOT / stem)
    info.update(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                compiler_sha256=hashlib.sha256((ROOT/'compiler.py').read_bytes()).hexdigest())
    (ROOT/metadata).write_text(json.dumps(info,indent=2,sort_keys=True)+'\n')
    print(stem,json.dumps(info,sort_keys=True))
