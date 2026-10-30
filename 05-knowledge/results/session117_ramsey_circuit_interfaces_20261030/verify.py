"""Reproduce bounded certificates; no external packages or downloaded code."""
from pathlib import Path
import hashlib, json, subprocess, sys, tempfile
ROOT=Path(__file__).resolve().parent
SCRIPTS=['circuits/check_point_defect.py','audit_circuits/independent_checks.py','affine_ramsey/check_affine_flats.py','audit_affine/check_independent.py','ramsey/check_ramsey_twins.py','audit_ramsey/check_independent.py','challenge_audit/check_potential_frontier.py','bilinear/check_bilinear.py']
def main():
    manifest=json.loads((ROOT/'manifest.json').read_text())
    for name,expected in manifest['sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=expected:
            raise RuntimeError('File hash mismatch: '+name)
    reports=[]
    for name in SCRIPTS:
        outputs=[]
        for opt in ([],['-O']):
            r=subprocess.run([sys.executable,*opt,str(ROOT/name)],capture_output=True,check=True)
            json.loads(r.stdout);outputs.append(r.stdout)
        if outputs[0]!=outputs[1]:raise RuntimeError('Optimization mismatch: '+name)
        digest=hashlib.sha256(outputs[0]).hexdigest()
        if digest!=manifest['expected_stdout_sha256'][name]:raise RuntimeError('Receipt mismatch: '+name)
        reports.append({'script':name,'normal_and_optimized':'PASS','stdout_sha256':digest})
    cpp='not requested'
    if '--cpp' in sys.argv:
        with tempfile.TemporaryDirectory() as tmp:
            exe=str(Path(tmp)/'four')
            subprocess.run(['g++','-O2','-std=c++17',str(ROOT/'circuits/check_four_inputs.cpp'),'-o',exe],check=True)
            result=subprocess.run([exe],capture_output=True,check=True)
            if result.stdout!=(ROOT/'circuits/check_four_inputs.out').read_bytes():raise RuntimeError('C++ receipt mismatch')
        cpp='PASS'
    print(json.dumps({'status':'PASS','python_checks':reports,'cpp_four_input':cpp,'scope':'finite regression certificates; see separate uniform proofs'},indent=2,sort_keys=True))
if __name__=='__main__':main()
