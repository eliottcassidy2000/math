#!/usr/bin/env python3
"""Replay the pinned producer as subprocesses, compare with independent arithmetic."""
import ast,hashlib,json,os,subprocess,sys,tempfile
from pathlib import Path
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parent

def require(ok,message):
    if not ok:raise RuntimeError(message)

def main():
    pins=json.loads((ROOT/'PRODUCER_PINS.json').read_text())
    for pin in pins['files']:
        b=(ROOT/'producer'/pin['path']).read_bytes()
        require(len(b)==pin['bytes'] and hashlib.sha256(b).hexdigest()==pin['sha256'],'producer byte pin')
    independent=json.loads((ROOT/'normal.json').read_text())['atlas']
    producer=json.loads((ROOT/'producer/blind_prediction.json').read_text())
    supplied=json.loads((ROOT/'producer/eight_moments.json').read_text())
    require(supplied==independent['moment_input']['moments'],'independent moments agree')
    own=independent['prediction']
    for key in ['lower','upper','prediction','candidate_count']:
        require(producer[key]==own[key],'producer independent agreement: '+key)
    require(producer['congruence']['residue']==own['residue'] and producer['congruence']['modulus']==own['modulus'],'producer congruence agrees')
    hand=json.loads((ROOT/'producer/hand_certificate.json').read_text());short=independent['short_certificate']
    require(hand['exact_upper']==short['upper'],'short exact evaluation agrees')
    require(hand['upper_polynomial_coefficients']==short['upper_polynomial_coefficients'],'short coefficient certificate agrees')
    require(hand['prediction']==short['prediction'],'short prediction agrees')
    modular=json.loads((ROOT/'producer/modular_prediction.json').read_text())
    require(modular['prediction']==independent['generic_twelve_moment_result']['prediction'],'universal residue-only result agrees')
    tree=ast.parse((ROOT/'producer/predictor.py').read_text())
    constants=[node.value for node in ast.walk(tree) if isinstance(node,ast.Constant)]
    require(own['prediction'] not in constants,'no literal target in producer predictor')
    # More decisive than a string scan: actual blind execution without any atlas available.
    env=os.environ.copy();env.pop('PYTHONPATH',None);env['PYTHONDONTWRITEBYTECODE']='1'
    with tempfile.TemporaryDirectory(prefix='s50-producer-blind-') as temp:
        p=Path(temp)
        for name in ['predictor.py','eight_moments.json']:(p/name).write_bytes((ROOT/'producer'/name).read_bytes())
        subprocess.run([sys.executable,'-I','predictor.py','--moments','eight_moments.json','--output','answer.json'],cwd=p,env=env,check=True,capture_output=True)
        require(json.loads((p/'answer.json').read_text())==producer,'actual isolated producer prediction')
    replay=json.loads(subprocess.check_output([sys.executable,'verify_packet.py'],cwd=ROOT/'producer',env=env,text=True))
    require(replay['status']=='PASS' and replay['fresh_directory_manifest_only'],'actual producer manifest-only replay')
    print(json.dumps({'status':'PASS','producer_pinned_files':len(pins['files']),'source_manifest_payload_files':replay['manifest_files'],'producer_normal_optimized_identical':replay['normal_optimized_identical'],'producer_blind_isolated_replay':True,'independent_exact_bounds_and_modular_agreement':True,'independent_short_certificate_agreement':True,'actual_producer_manifest_only_replay':True},indent=2,sort_keys=True))
if __name__=='__main__':main()
