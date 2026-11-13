#!/usr/bin/env python3
"""Replay a frozen review against its exact pre-promotion document snapshots.
The live worktree is never modified. Mathematical producers run directly.
"""
from pathlib import Path
from hashlib import sha256
import argparse,json,shutil,subprocess,sys,tempfile
P=Path(__file__).resolve().parent
ROOT=P.parents[2]
parser=argparse.ArgumentParser()
parser.add_argument('review',choices=['bridges/check_reconstruction_audit.py','reconstruction/check_final_bridge_audit.py'])
parser.add_argument('--optimized',action='store_true')
args=parser.parse_args()
with tempfile.TemporaryDirectory(prefix='integer32-reviewed-') as temp:
    root=Path(temp)/'repo'
    snapshot=root/'05-knowledge/results'/P.name
    shutil.copytree(P,snapshot,ignore=shutil.ignore_patterns('__pycache__'))
    for name,row in json.loads((P/'promotions.json').read_text())['files'].items():
        f=snapshot/name;s=f.read_text()
        if sha256(s.encode()).hexdigest()!=row['after_sha256']:raise RuntimeError('current promotion '+name)
        for edit in reversed(row['edits']):
            if s.count(edit['new'])!=edit['count']:raise RuntimeError('ambiguous inverse '+name)
            s=s.replace(edit['new'],edit['old'])
        if sha256(s.encode()).hexdigest()!=row['before_sha256']:raise RuntimeError('snapshot '+name)
        f.write_text(s)
    # The inherited data are immutable inputs, not another copied census.
    (snapshot.parent/'incidence_closure_20261113').symlink_to(ROOT/'05-knowledge/results/incidence_closure_20261113',target_is_directory=True)
    command=[sys.executable,'-B']+(['-O'] if args.optimized else [])+[str(snapshot/args.review)]
    subprocess.run(command,cwd=root,check=True)
