#!/usr/bin/env python3
"""Portable launcher for the unchanged pinned independent review."""
import argparse, pathlib, shutil, subprocess, sys, tempfile
HERE=pathlib.Path(__file__).resolve().parent
p=argparse.ArgumentParser()
p.add_argument('--repo',type=pathlib.Path)
p.add_argument('--optimized',action='store_true')
p.add_argument('--negative',action='store_true')
a=p.parse_args()
repo=a.repo or next((x for x in HERE.parents if (x/'AGENTS.md').is_file() and (x/'01-canon').is_dir()),None)
if repo is None:raise SystemExit('Pass --repo /absolute/path/to/math')
repo=repo.resolve();flags=['-B']+(['-O'] if a.optimized else [])
if not a.negative:
    subprocess.run([sys.executable,*flags,str(HERE/'review/independent_audit.py'),'--repo',str(repo),'--producer',str(HERE)],check=True)
else:
    with tempfile.TemporaryDirectory(prefix='metagraph-review-layout-') as temp:
        root=pathlib.Path(temp)
        (root/'math').symlink_to(repo,target_is_directory=True)
        (root/'session20-metagraph-networks').symlink_to(HERE,target_is_directory=True)
        review=root/'review';review.mkdir()
        for name in ['negative_pin_tests.py','independent_audit.py','frozen_inputs.json']:
            shutil.copy2(HERE/'review'/name,review/name)
        subprocess.run([sys.executable,*flags,str(review/'negative_pin_tests.py')],check=True)
