"""Freeze four proof lanes, force fresh normal/-O replays, and bind outputs."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[2]
STEM='prime_shells_20260921'
EXP=Path('04-computation/experiments')
RES=Path('05-knowledge/results')
LANES={'shells':3834,'primes':19605,'dissections':25846,'elliptic':30856}
MANIFEST=RES/f'{STEM}_manifest.json'
INHERITED=(
 '01-canon/theorems/THM-3756-odd-square-ordinal-berggren-affine-descent.md',
 '01-canon/theorems/THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit.md',
 '01-canon/theorems/THM-3991-periodic-unimodular-toric-cusp-factorial-euler-obstruction.md',
 '05-knowledge/results/odd_square_20260921_triangles.md',
 '05-knowledge/results/arithmetic_braids_20260917_summand.md',
 '05-knowledge/results/collatz_mod6_20260917_row_braid_typing.md',
 '05-knowledge/results/catalan_elliptic_20260921_elliptic.md',
 '05-knowledge/results/catalan_elliptic_20260921_catalan.md',
)


def require(condition,message):
    if not condition:
        raise RuntimeError(message)


def digest(path,lf=True):
    data=(ROOT/path).read_bytes()
    if lf:
        data=data.replace(b'\r\n',b'\n').replace(b'\r',b'\n')
    return sha256(data).hexdigest()


def jd(value):
    return sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def save(record):
    (ROOT/MANIFEST).write_text(json.dumps(record,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')


def check_certificate(lane,value):
    require(isinstance(value,dict),'certificate must be an object')
    require(value.get('status') in ('FINITE-EXACT','FINITE-EXACT PASS','PASS',
            'FINITE-EXACT; all-quantifier proofs are in the note'),f'incomplete {lane}')
    require(value.get('checks',value.get('checks_passed'))==LANES[lane],f'changed check universe {lane}')
    if lane=='shells':
        require(value['universe']['all_odd_shells']==[3,501],'shell range')
        require(value['cycles17']==[[1,15,13,9],[3,11,5,7]],'shell controls')
    elif lane=='primes':
        require(value['strict_prime_count']==155559,'prime interval')
        require(value['H4']=='25/12' and value['H6']=='49/20','harmonic controls')
    elif lane=='elliptic':
        require(value['repaired_digits']==[81,80,79],'seed sizes')
        require(value['positive_nG_through_40']==[9,17],'positivity control')
    field='source_sha256_lf' if 'source_sha256_lf' in value else 'source_sha256'
    require(value[field]==digest(EXP/f'{STEM}_{lane}.py',lf=field.endswith('_lf')),f'stale source {lane}')


def main():
    record={'status':'RUNNING','started_utc':datetime.now(timezone.utc).isoformat(timespec='seconds'),
      'commands':[],'lanes':{},'hash_basis':'UTF8 bytes with CRLF/CR normalized to LF unless an embedded raw hash is declared',
      'scope':'Four finite exact replays; infinite statements have separate written proofs; cited classifications remain external.',
      'collatz_proved':False,'new_exceptional_prime_discovered':False,'new_Lean_formalization':False,
      'manifest_self_hash':'Excluded; verifier source included.'}
    save(record)
    inputs=[Path(__file__).relative_to(ROOT),RES/f'{STEM}_synthesis.md',*[Path(p) for p in INHERITED]]
    for lane in LANES:
        inputs.extend([EXP/f'{STEM}_{lane}.py',RES/f'{STEM}_{lane}.md'])
    try:
        frozen={p.as_posix():digest(p) for p in inputs}
        for lane in LANES:
            # Positive and hostile gate controls ensure stale sentinels cannot pass.
            try:
                check_certificate(lane,{'status':'MASTER_REPLAY_RUNNING'})
            except RuntimeError:
                pass
            else:
                raise RuntimeError('sentinel accepted')
            values=[]
            for optimized in (False,True):
                script=EXP/f'{STEM}_{lane}.py'
                cert=script.with_suffix('.json')
                (ROOT/cert).write_text('{"status":"MASTER_REPLAY_RUNNING"}\n',encoding='utf-8',newline='\n')
                command=[sys.executable,*(['-O'] if optimized else []),str(ROOT/script)]
                shown=['python',*(['-O'] if optimized else []),script.as_posix()]
                print('RUN: '+' '.join(shown),flush=True)
                result=subprocess.run(command,cwd=ROOT,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=180)
                output=(result.stdout+result.stderr).strip()
                for spelling in (str(ROOT),ROOT.as_posix()):
                    output=output.replace(spelling,'<repo>')
                record['commands'].append({'command':shown,'exit_code':result.returncode,'output':output})
                save(record)
                require(result.returncode==0 and 'PASS' in result.stdout,f'replay failed {lane}: {output}')
                value=json.loads((ROOT/cert).read_text(encoding='utf-8'))
                check_certificate(lane,value)
                values.append(value)
            require(values[0]==values[1],f'normal/-O mismatch {lane}')
            record['lanes'][lane]={'status':'PASS','normal_optimized_equal':True,
                'checks':LANES[lane],'decoded_json_sha256':jd(values[-1]),
                'source':script.as_posix(),'certificate':cert.as_posix(),
                'proof_note':(RES/f'{STEM}_{lane}.md').as_posix(),
                'universe':values[-1].get('universe',values[-1].get('universes',values[-1].get('scope','See certificate and note')))}
            save(record)
        require(all(digest(p)==frozen[p.as_posix()] for p in inputs),'frozen input changed during replay')
        outputs=[EXP/f'{STEM}_{lane}.json' for lane in LANES]
        for lane in LANES:
            value=json.loads((ROOT/EXP/f'{STEM}_{lane}.json').read_text(encoding='utf-8'))
            require(jd(value)==record['lanes'][lane]['decoded_json_sha256'],'output changed after replay')
        record['artifacts_sha256_lf']={p.as_posix():digest(p) for p in sorted(inputs+outputs)}
        record.update(status='PASS',finished_utc=datetime.now(timezone.utc).isoformat(timespec='seconds'))
        save(record)
        print('PASS: four fresh normal/-O pairs; all source, proof and output hashes bound.')
    except Exception as error:
        record.update(status='FAIL',error=f'{type(error).__name__}: {error}')
        save(record)
        raise


if __name__=='__main__':
    main()
