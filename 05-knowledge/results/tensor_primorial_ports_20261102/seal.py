"""Pin inherited inputs and verify normalized source/packet hashes."""
from pathlib import Path
import hashlib,json,subprocess,sys
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SOURCES=[
'05-knowledge/results/stability_748_guarded_encodings_20261102/STABILITY_48.md',
'05-knowledge/results/session139_incidence_moduli_and_stability_20261102/independent_concurrent48/AUDIT.md',
'05-knowledge/results/session139_incidence_moduli_and_stability_20261102/tensor_matching/PROOF.md',
'05-knowledge/results/session139_incidence_moduli_and_stability_20261102/stability/EISENSTEIN_BERNSTEIN_TENSOR_PROOF.md',
'05-knowledge/results/stability_priority_moduli_20261102/STABILITY_N276.md',
'05-knowledge/results/degree_local_tournament_enumeration_20261021/session22-enumeration/local_scale.py',
'05-knowledge/results/degree_local_tournament_enumeration_20261021/session22-enumeration/README.md',
'05-knowledge/results/regular_parent_symmetry_20261029/verify.py',
'05-knowledge/results/session114_continuation_symmetry_20261030/lrc/GENERAL_CORE.md',
'01-canon/theorems/THM-1286-D-graded-gate-tower-primorial-cascade.md',
'05-knowledge/results/session139_incidence_moduli_and_stability_20261102/all_moduli/PROOF.md',
'05-knowledge/results/free_bits_seven_fiftysix_20261028/FIFTYSIX.md',
'05-knowledge/results/aliquot_divisor_overlap_20261028/README.md',
]
def digest(p):return hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def write(p,obj):p.write_text(json.dumps(obj,sort_keys=True,indent=2)+'\n',encoding='utf-8')
def main():
    if '--check' in sys.argv:
        pins=json.loads((HERE/'SOURCE_PINS.json').read_text(encoding='utf-8'))
        for p,h in pins['sha256_LF'].items():
            if digest(ROOT/p)!=h:raise RuntimeError('changed inherited input: '+p)
        seal=json.loads((HERE/'SHA256SUMS.json').read_text(encoding='utf-8'))
        for p,h in seal.items():
            if digest(HERE/p)!=h:raise RuntimeError('changed packet artifact: '+p)
        for stem in ['verification','recurrences']:
            if (HERE/(stem+'.json')).read_bytes()!=(HERE/(stem+'.optimized.json')).read_bytes():raise RuntimeError('optimized replay differs')
        print(json.dumps({'status':'PASS','source_pins':len(pins['sha256_LF']),'sealed_artifacts':len(seal)}))
    else:
        revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
        write(HERE/'SOURCE_PINS.json',{'revision':revision,'normalization':'SHA256 after CRLF to LF','sha256_LF':{p:digest(ROOT/p) for p in SOURCES}})
        write(HERE/'SHA256SUMS.json',{p.name:digest(p) for p in sorted(HERE.iterdir()) if p.is_file() and p.name!='SHA256SUMS.json'})
        print('sealed')
if __name__=='__main__':main()
