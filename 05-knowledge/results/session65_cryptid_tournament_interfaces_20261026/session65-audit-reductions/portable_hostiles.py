#!/usr/bin/env python3
"""Package portability/mutation tests only; no new mathematics and no Git."""
from pathlib import Path
import json
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def run(packet, optimized):
    command = [sys.executable, '-I', '-B'] + (['-O'] if optimized else [])
    command += [str(packet/'session65-audit-reductions'/'portable_check.py')]
    return subprocess.run(command, cwd=packet, env={'PATH':'/nonexistent'},
                          capture_output=True, text=True, check=False)

with tempfile.TemporaryDirectory(prefix='session65-offline-package-') as temp:
    packet = Path(temp)
    shutil.copytree(ROOT.parent/'session65-reductions', packet/'session65-reductions')
    shutil.copytree(ROOT, packet/'session65-audit-reductions')
    need(not (packet/'session65-repo').exists(), 'unexpected source checkout')
    need(not any(x.name == '.git' for x in packet.rglob('*')), 'unexpected Git metadata')
    positives = [run(packet, opt) for opt in (False, True)]
    need(all(p.returncode == 0 and not p.stderr for p in positives), 'isolated replay failed')
    need(positives[0].stdout == positives[1].stdout, 'isolated normal/-O disagree')
    result = json.loads(positives[0].stdout)
    old = json.loads((ROOT/'normal.json').read_text())
    keys = ('codec','continuation_hostile','lrc_threshold_controls')
    need(all(result[key] == old[key] for key in keys), 'historical math consequence changed')
    need(result['source_integrity']['fresh_git_readback_performed'] is False,
         'offline runtime misrepresented as Git readback')
    pins = json.loads((packet/'session65-reductions'/'SOURCE_PINS.json').read_text())
    hostile = []
    for row in pins['files']:
        path = packet/'session65-reductions'/row['snapshot']
        original = path.read_bytes()
        for mutation in ('wrong_bytes','missing_file'):
            try:
                if mutation == 'wrong_bytes':
                    path.write_bytes(original+b'\nOFFLINE_HOSTILE_CORRUPTION\n')
                else:
                    path.unlink()
                outcomes=[]
                for opt in (False,True):
                    p=run(packet,opt)
                    need(p.returncode != 0 and not p.stdout, 'bad snapshot accepted: '+row['snapshot'])
                    need(path.name in p.stderr, 'failure did not identify attacked snapshot')
                    outcomes.append({'mode':'optimized' if opt else 'normal','exit_code':p.returncode,
                                     'error':p.stderr.strip().splitlines()[-1].replace(str(packet),'<isolated-package>')})
                hostile.append({'snapshot':row['snapshot'],'mutation':mutation,'outcomes':outcomes})
            finally:
                path.write_bytes(original)
    # An altered recorded blob ID must not be accepted merely because bytes match.
    provenance=packet/'session65-audit-reductions'/'OFFLINE_SOURCE_PINS.json'
    original=provenance.read_bytes()
    edited=json.loads(original); edited['files'][0]['git_blob_oid_sha1']='0'*40
    provenance.write_text(json.dumps(edited))
    outcomes=[]
    for opt in (False,True):
        p=run(packet,opt)
        need(p.returncode != 0 and not p.stdout and 'offline provenance pin manifest changed' in p.stderr,
             'tampered offline Git blob pin accepted')
        outcomes.append({'mode':'optimized' if opt else 'normal','exit_code':p.returncode})
    provenance.write_bytes(original)
    need(run(packet,False).stdout == positives[0].stdout, 'restored package failed')
    print(json.dumps({'status':'PASS_OFFLINE_PACKAGE_AND_HOSTILES',
                      'runtime_environment':{'git_checkout_present':False,'git_on_PATH':False,
                                             'ambient_PYTHONPATH_used':False,'python_isolated_mode':True},
                      'normal_optimized_byte_equal':True,
                      'historical_math_consequence_keys_equal':list(keys),
                      'source_snapshot_hostile_cases':len(hostile),
                      'source_snapshot_hostile_runs':2*len(hostile),
                      'snapshot_hostiles':hostile,'tampered_provenance_manifest_rejected':outcomes},
                     sort_keys=True,indent=2))
