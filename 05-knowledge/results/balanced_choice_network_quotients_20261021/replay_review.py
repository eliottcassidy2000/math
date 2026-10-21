#!/usr/bin/env python3
"""Portable paths for the unchanged independent choice-network audit."""
import argparse,pathlib,subprocess,sys
HERE=pathlib.Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--repo',type=pathlib.Path);p.add_argument('--optimized',action='store_true');a=p.parse_args()
repo=a.repo or next((x for x in HERE.parents if (x/'AGENTS.md').is_file() and (x/'01-canon').is_dir()),None)
if repo is None:raise SystemExit('Pass --repo /absolute/path/to/math')
flags=['-B']+(['-O'] if a.optimized else [])
subprocess.run([sys.executable,*flags,str(HERE/'review/check_review.py'),'--producer',str(HERE),'--repo',str(repo.resolve()),'--workspace',str(HERE/'review/workspace')],check=True)
