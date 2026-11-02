"""Freeze or check this research packet and the exact inherited inputs it uses."""
from pathlib import Path
import hashlib,json,sys,subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
MANIFEST=HERE/'MANIFEST.json'
INPUTS=[
 '05-knowledge/results/collatz_source_refinement_floor_20261005.md',
 '05-knowledge/results/collatz_tail_quotient_integration_20261009.md',
 '05-knowledge/results/collatz_coherence_20260930_trichotomy_hexagon_motifs.md',
 '05-knowledge/results/regular_parent_symmetry_20261029/verify.py',
 '05-knowledge/results/hadwiger_six_seven_eight_20261025/fano_eight/RESULT.md',
 '05-knowledge/results/hadwiger_six_seven_eight_20261025/metric_spindle/RESULT.md',
 '05-knowledge/results/session62_agent_theme_bridges_20261025/holonomy/TWISTED_PAIR_LIFTS.md',
 '05-knowledge/results/session61_hadwiger_six_seven_20261025/SOURCE_STATUS.md',
 '05-knowledge/results/lucas_mirror_gap_spaces_20261028/GAP_TOPOLOGY.md',
 '05-knowledge/results/stability_priority_moduli_20261102/FRACTRAN_PRECISION.md',
 '05-knowledge/results/clocks_carriers_genus_20261102/COLLATZ_GEOMETRY.md',
 '05-knowledge/results/joint_incidence_clocked_frames_20261102/COLLATZ_LATTICE.md',
 '05-knowledge/results/session140_collatz_hexagonal_geometry_20261102/collatz_hex/PROOF.md',
 '01-canon/theorems/THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit.md',
]

def digest(path):
    data=path.read_bytes()
    if path.suffix in ('.md','.py','.json','.tsv','.svg'):data=data.replace(b'\r\n',b'\n')
    return hashlib.sha256(data).hexdigest()

def current():
    files=[p for p in HERE.iterdir() if p.is_file() and p!=MANIFEST]
    return {
        'packet':{p.name:digest(p) for p in sorted(files)},
        'inputs':{name:digest(ROOT/name) for name in INPUTS},
    }

def main():
    if sys.argv[1:]==['--write']:
        data=current()
        data['input_checkout_commit']=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
        data['hash_convention']='SHA256; CRLF normalized to LF for md/py/json/tsv/svg; other files literal bytes'
        MANIFEST.write_text(json.dumps(data,sort_keys=True,indent=2)+'\n',encoding='utf-8')
        print('WROTE manifest after explicit input review')
    elif sys.argv[1:] in ([],['--check']):
        saved=json.loads(MANIFEST.read_text(encoding='utf-8'));now=current()
        errors=[]
        for section in ('packet','inputs'):
            for name in sorted(set(saved[section])|set(now[section])):
                if saved[section].get(name)!=now[section].get(name):errors.append(section+': '+name)
        if errors:raise RuntimeError('Seal mismatch\n'+'\n'.join(errors))
        a=(HERE/'verification.json').read_bytes().replace(b'\r\n',b'\n')
        b=(HERE/'verification.optimized.json').read_bytes().replace(b'\r\n',b'\n')
        if a!=b:raise RuntimeError('normal and optimized receipts differ')
        print('PASS: packet and inherited inputs sealed; normal/optimized receipts identical')
    else:raise SystemExit('Usage: seal.py [--check|--write]')

if __name__=='__main__':main()
