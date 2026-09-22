"""Replay four odd-square lanes; bind finite certificates to written proofs.

Run from any directory. Explicit exceptions keep all gates active under -O.
No finite replay proves Collatz or substitutes for the notes' infinite proofs.
"""
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[2]
STEM='odd_square_20260921'
EXPERIMENTS=Path('04-computation/experiments')
RESULTS=Path('05-knowledge/results')
LANES=('triangles','shapes','edges','inverse')
MANIFEST=RESULTS/f'{STEM}_manifest.json'
ARCHIVE=Path('05-knowledge/reference/ODD-SQUARE-GLUED-LINE-2026-09-21-SOURCE.md')
INHERITED=(
    Path('01-canon/theorems/THM-3756-odd-square-ordinal-berggren-affine-descent.md'),
    RESULTS/'arithmetic_braids2_20260917_signed_cycles.md',
    RESULTS/'arithmetic_braids2_20260917_inverse_completion.md',
    RESULTS/'arithmetic_braids2_20260917_floor_reciprocity.md',
    RESULTS/'arithmetic_braids_20260917_divisors.md',
)


def require(condition,message):
    if not condition:
        raise RuntimeError(message)


def stamp():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def file_hash(path,lf=True):
    data=(ROOT/path).read_bytes()
    data.decode('utf-8')
    if lf:
        data=data.replace(b'\r\n',b'\n').replace(b'\r',b'\n')
    return sha256(data).hexdigest()


def json_hash(value):
    return sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def save(record):
    (ROOT/MANIFEST).write_text(json.dumps(record,indent=2,sort_keys=True)+'\n',
                               encoding='utf-8',newline='\n')


def check_marker(value,lane):
    status=value.get('status','')
    require(status=='PASS' or 'FINITE-EXACT' in status,f'Missing final status: {lane}')
    if lane=='triangles':
        require(value['checks_passed']==368573 and
                value['universe']['primitive_root_pairs']==9242,'Triangle universe changed')
    elif lane=='shapes':
        require(value['checks']==480912 and
                value['scope']['all_triples_enumerated']==159139,'Shape universe changed')
    elif lane=='edges':
        require(value['checks_passed']==813891 and
                value['universe']['positive_odd_starts']==[1,100001],'Edge universe changed')
    else:
        require(value['primitive_mixed_completions']['count']==14040 and
                value['primitive_mixed_completions']['additional_family_checks']==14040,
                'Inverse completion universe changed')


def replay(lane,record):
    script=EXPERIMENTS/f'{STEM}_{lane}.py'
    certificate=script.with_suffix('.json')
    values=[]
    for optimized in (False,True):
        (ROOT/certificate).write_text('{"status":"MASTER_REPLAY_RUNNING"}\n',
                                       encoding='utf-8',newline='\n')
        command=[sys.executable,*(['-O'] if optimized else []),str(ROOT/script)]
        shown=['python',*(['-O'] if optimized else []),script.as_posix()]
        print('RUN: '+' '.join(shown),flush=True)
        result=subprocess.run(command,cwd=ROOT,capture_output=True,text=True,
                              encoding='utf-8',errors='replace',timeout=180,check=False)
        diagnostic=result.stdout+result.stderr
        for spelling in sorted({str(ROOT),ROOT.as_posix(),json.dumps(str(ROOT))[1:-1]},
                               key=len,reverse=True):
            diagnostic=diagnostic.replace(spelling,'<repo>')
        record['commands'].append({'command':shown,'cwd':'.',
                                  'exit_code':result.returncode,'output':diagnostic.strip()})
        save(record)
        require(result.returncode==0 and 'PASS' in result.stdout,
                f'Failed or incomplete replay {lane}: {diagnostic}')
        value=json.loads((ROOT/certificate).read_text(encoding='utf-8'))
        require(isinstance(value,dict),f'Invalid certificate object: {lane}')
        check_marker(value,lane)
        fields=[]
        for field in ('source_sha256','script_sha256','source_sha256_lf'):
            if field in value:
                require(value[field]==file_hash(script,lf=field.endswith('_lf')),
                        f'Stale embedded source hash: {lane}:{field}')
                fields.append(field)
        values.append(value)
    require(values[0]==values[1],f'Normal/-O decoded JSON mismatch: {lane}')
    record['lanes'][lane]={
        'status':'PASS','normal_optimized_decoded_json_equal':True,
        'script':script.as_posix(),'certificate':certificate.as_posix(),
        'note':(RESULTS/f'{STEM}_{lane}.md').as_posix(),
        'embedded_hash_fields_checked':fields,
        'source_frozen_and_bound_by_master':True,
        'decoded_json_sha256':json_hash(values[1]),
        'finite_universe':values[1].get('universe',values[1].get('scope',{
            'specified_in_certificate':certificate.as_posix(),
            'sections':sorted(values[1]),'no_extrapolation':True})),
    }
    save(record)
    return certificate


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    record={
        'status':'RUNNING','started_utc':stamp(),'commands':[],'lanes':{},
        'scope':'Four exact finite replays with separately proved infinite statements and explicit inherited dependencies.',
        'collatz_convergence_proved':False,'all_signed_cycles_classified':False,
        'new_general_Lean_formalization':False,
        'cited_context':{
            'Krasikov_Lagarias':'https://arxiv.org/pdf/math/0205002',
            'pin':'Section2 counting definitions; Theorem6.1 and proof in section6',
            'role':'Historical inverse-orbit lower bound; no current-best or universal-coverage claim.'},
        'hash_basis':'SHA256 UTF8 bytes with CRLF and lone CR normalized to LF; embedded hashes checked in their declared basis.',
        'manifest_self_hash':'Excluded to avoid circularity; master source is pinned.',
    }
    save(record)
    try:
        inputs=[Path(__file__).relative_to(ROOT),ARCHIVE,RESULTS/f'{STEM}_synthesis.md',*INHERITED]
        for lane in LANES:
            inputs.extend([EXPERIMENTS/f'{STEM}_{lane}.py',RESULTS/f'{STEM}_{lane}.md'])
        frozen={p.as_posix():file_hash(p) for p in inputs}
        outputs=[replay(lane,record) for lane in LANES]
        for p in inputs:
            require(file_hash(p)==frozen[p.as_posix()],f'Input changed during replay: {p}')
        for lane,details in record['lanes'].items():
            value=json.loads((ROOT/details['certificate']).read_text(encoding='utf-8'))
            require(json_hash(value)==details['decoded_json_sha256'],
                    f'Certificate changed after replay: {lane}')
        record['artifacts_sha256_lf']={p.as_posix():file_hash(p) for p in sorted(set(inputs+outputs))}
        require(len(record['lanes'])==4,'Missing completed lane')
        record.update(status='PASS',finished_utc=stamp())
        save(record)
        print('PASS: four fresh normal/-O pairs agree; certificates, notes and dependencies pinned.')
    except Exception as error:
        record.update(status='FAIL',finished_utc=stamp(),error=f'{type(error).__name__}: {error}')
        save(record)
        raise


if __name__=='__main__':
    main()
