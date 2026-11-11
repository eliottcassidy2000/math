#!/usr/bin/env python3
"""Replay new-artifact controls. Historical claims are inherited, not re-audited."""
import argparse,hashlib,json,subprocess,sys,tempfile,time,importlib.util
from pathlib import Path
P=Path(__file__).resolve().parent

def run():
    parser=argparse.ArgumentParser();parser.add_argument('--full',action='store_true',help='also replay the complete atlas and independent subset-DP controls')
    args=parser.parse_args();rows=[]
    def command(label,cmd,cwd=P):
        t=time.perf_counter();r=subprocess.run(cmd,cwd=cwd,text=True,capture_output=True)
        rows.append({'check':label,'exit_code':r.returncode,'seconds':time.perf_counter()-t,'stdout_sha256':hashlib.sha256(r.stdout.encode()).hexdigest()})
        if r.returncode:raise RuntimeError(label+'\n'+r.stdout+'\n'+r.stderr)
        return r.stdout
    for lane in ('enumeration','chanpak'):
        command(lane+' checksums',['sha256sum','-c','SHA256SUMS'],P/lane)
    command('final exact enumeration controls',[sys.executable,'enumeration/verify_final.py'])
    command('new even-order residue formulas',[sys.executable,'bridges/even_mod8.py'])
    command('marked classes and labelled path mass',[sys.executable,'carrier/check.py'])
    command('positive transfer small physical controls',[sys.executable,'chanpak/spine.py'])
    command('averaged GF and tomography',[sys.executable,'chanpak/moments.py'])
    command('all-word rational rank certificate',[sys.executable,'chanpak/rank.py'])
    witness=json.loads(command('compile target999999',[sys.executable,'chanpak/compile_target.py','999999']))
    if witness['target']!=999999 or witness['vertices']>16:raise RuntimeError('witness endpoint')
    if args.full:
        with tempfile.TemporaryDirectory(prefix='s270-replay-') as d:
            binary=str(Path(d)/'atlas.bin');exe=str(Path(d)/'verify')
            command('decode full atlas',[sys.executable,'chanpak/compile_target.py','--decode',binary])
            command('compile full atlas verifier',['g++','-O3','-std=c++17','chanpak/verify_atlas.cpp','-o',exe])
            output=json.loads(command('complete atlas replay and subset checks',[exe,binary]))
            if output['all_certificate_transfer_replays']!=499998:raise RuntimeError('incomplete atlas')
    receipt={'status':'PASS','full_atlas_replayed':args.full,'checks':rows,'historical_inputs':'Inherited, not re-audited. New outputs independently checked in stated scopes.'}
    (P/'FINAL_CHECKS.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':run()
