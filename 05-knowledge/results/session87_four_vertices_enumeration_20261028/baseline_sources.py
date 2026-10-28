"""Hash-verified access to existing repository sources; no vendored copies.

--baseline-root may point to a repository or a flat staging directory holding
these exact sources. The published packet resolves the repository by ancestry.
"""
from hashlib import sha256
from pathlib import Path

SOURCES = {
    'local_scale.py': ('05-knowledge/results/degree_local_tournament_enumeration_20261021/session22-enumeration/local_scale.py',
        'b76167ef82dda0b134681abce33e159d5dd80d15e92edb9e31ab9b490f6939ad'),
    'channels.py': ('05-knowledge/results/converse_tournament_count_channels_20261021/session22-converse-channels/channels.py',
        'b27d30d8f8cc81768c8f6b382f446658515d86f2cf5f212d2b59e20859ba6792'),
    'recreate_a000568.py': ('05-knowledge/results/session86_relations_attachment_carriers_20261028/recreate_a000568.py',
        'a8d85c9bc0e355ed4bf8140f6ae23b6d2c4141aad51d1b8bf3c067f302ee564f'),
    'a000568_turbo.py': ('04-computation/a000568_turbo.py',
        '584535fa7429b729b9e9c2ad30867eedcd31b3b0b46235238b0bd76ce2fd1888'),
    'a000568_gmp_enum_v4.c': ('04-computation/a000568_gmp_enum_v4.c',
        'b9d6445b8dd7833e25fecb23707155b4e3d0a104888f3a1249fce96d4efb305a'),
}
ROOT = None


def configure(root=None):
    global ROOT
    if root is not None:
        ROOT = Path(root).resolve()
        return
    for p in Path(__file__).resolve().parents:
        if (p/'AGENTS.md').is_file() and (p/'05-knowledge/results').is_dir():
            ROOT = p
            return
    raise FileNotFoundError('Cannot locate repository ancestor; provide --baseline-root')


def source_path(name):
    if ROOT is None:
        configure()
    relative,expected = SOURCES[name]
    path = ROOT/relative
    if not path.is_file():
        path = ROOT/name  # explicit flat staging override; hash still enforced
    actual = sha256(path.read_bytes()).hexdigest()
    if actual != expected:
        raise ValueError(f'Historical baseline changed: {relative}; expected {expected}, got {actual}')
    return path
