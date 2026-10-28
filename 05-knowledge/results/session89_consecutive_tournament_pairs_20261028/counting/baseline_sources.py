"""Read hash-pinned baselines from existing files or one exact Git blob.

No repository writes, broad searches, source rewriting, or vendored baselines.
"""
from hashlib import sha256
from pathlib import Path
import subprocess
import types

COMMIT = '03408f84ac43fcdcfb0887e5d6bc1fc8792415cf'
SOURCES = {
 'session87': ('05-knowledge/results/session87_four_vertices_enumeration_20261028/fast_counts.py',
 '3844a4a590daaf01410cf99d3584441bfa2f8e5e95371ccaafe534fd152b0ece'),
 'session22': ('05-knowledge/results/converse_tournament_count_channels_20261021/session22-converse-channels/channels.py',
 'b27d30d8f8cc81768c8f6b382f446658515d86f2cf5f212d2b59e20859ba6792'),
}


def find_root():
    here=Path(__file__).resolve()
    for parent in here.parents:
        if (parent/'AGENTS.md').is_file() and (parent/'05-knowledge/results').is_dir():
            return parent
        candidate=parent/'session88-symmetry/verified-repo'
        if (candidate/'AGENTS.md').is_file():
            return candidate
    raise FileNotFoundError('pass --baseline-root for the repository')


def load(name, root=None):
    root=Path(root) if root else find_root()
    relative,expected=SOURCES[name]
    path=root/relative
    if path.is_file():
        source=path.read_bytes()
    else:
        source=subprocess.run(['git','-C',str(root),'show',f'{COMMIT}:{relative}'],
                              check=True,capture_output=True).stdout
    actual=sha256(source).hexdigest()
    if actual != expected:
        raise ValueError(f'{name} hash mismatch: expected {expected}, got {actual}')
    module=types.ModuleType(name)
    exec(compile(source,str(path),'exec'),module.__dict__)
    return module
