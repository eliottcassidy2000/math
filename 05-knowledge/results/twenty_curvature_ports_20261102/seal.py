"""Canonical LF hashes for reproducible research inputs and packet files."""
from pathlib import Path
import argparse
import hashlib
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
INPUTS=[
    '05-knowledge/results/joint_incidence_clocked_frames_20261102/path_kernel.py',
    '05-knowledge/results/joint_incidence_clocked_frames_20261102/PATH_KERNEL.md',
    '05-knowledge/results/lower_frontier_complement_fibonacci_20261102/COMPLEMENT_CARRIER.md',
    '05-knowledge/results/regular_parent_symmetry_20261029/verify.py',
    '05-knowledge/results/collatz_hexagonal_growth_20261102/verify.py',
    '05-knowledge/results/collatz_hexagonal_growth_20261102/JOINT_PORTS.md',
    '05-knowledge/results/session140_collatz_hexagonal_geometry_20261102/curvature_control/PROOF.md',
    '05-knowledge/results/session140_collatz_hexagonal_geometry_20261102/README.md',
    '05-knowledge/results/aliquot_divisor_overlap_20261028/THEORY.md',
    '05-knowledge/results/collatz_mod6_20260917_three_adic_g_map.md',
    '05-knowledge/results/cursor_cm_phase_20261028/TRACE_CARRIER.md',
    '05-knowledge/results/cursor_cm_phase_20261028/RAMANUJAN.md',
    '05-knowledge/results/scaled_count_fft_carriers_20261030/THEORY.md',
    '05-knowledge/results/free_bits_seven_fiftysix_20261028/ATLAS.md',
    '05-knowledge/results/twin_pythagorean_124_20261026/CONNECTION_BOARD.md',
    '05-knowledge/results/session141_metric_aliquot_square_paths_20261102/curvature_triangle/PROOF.md',
    '05-knowledge/results/session141_metric_aliquot_square_paths_20261102/square_paths/PROOF.md',
    '05-knowledge/results/session141_metric_aliquot_square_paths_20261102/square_paths_content/PROOF.md',
    '05-knowledge/results/session141_metric_aliquot_square_paths_20261102/CONNECTIONS.md',
    '05-knowledge/results/session141_metric_aliquot_square_paths_20261102/aliquot_twenty/PROOFS.md',
    '05-knowledge/results/eleven_wedge_honeycomb_20261102/README.md',
    '05-knowledge/results/honeycomb_gap_aliquot_20261102/GAP_GEOMETRY.md',
    '01-canon/theorems/THM-4661-odd-degree-carriers-a-failed-six-colour-exterior-has-at-least-twenty-vertices-and-equality-testers-are-seven-colour-certificates.md',
]

def digest(p):
    raw=p.read_bytes()
    if p.suffix in {'.md','.py','.json','.svg','.tsv','.txt'}:
        raw=raw.replace(b'\r\n',b'\n')
    return hashlib.sha256(raw).hexdigest()

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    inputs={s:digest(ROOT/s) for s in INPUTS}
    files={p.name:digest(p) for p in sorted(HERE.iterdir())
           if p.is_file() and p.name!='MANIFEST.json'}
    result={'hash_convention':'SHA256; canonical LF for listed text extensions, raw bytes otherwise',
            'starting_commit':'fb7c88537a4a','inputs':inputs,'files':files,
            'primary_sources':[
                {'url':'https://arxiv.org/html/2412.11914v2','scope':'Theorem1 and Table1 finite unit distances; read 2026-11-02'},
                {'url':'https://www.pollack-math.net/pnp-4.pdf','scope':'near-perfect definition, Proposition1, counterexample650; read 2026-11-02'},
                {'url':'https://proceedings.mlr.press/v80/sala18a.html','scope':'classical tree-embedding context; read 2026-11-02'},
                {'url':'https://www.maths.dur.ac.uk/users/anna.felikson/Geometry/Geometry23/notes.pdf',
                 'scope':'standard metric formula via pinned Session140 proof; direct refetch timed out this session'}]}
    target=HERE/'MANIFEST.json'
    if args.check:
        expected=json.loads(target.read_text(encoding='utf-8'))
        if expected!=result:
            changed=[s for kind in ('inputs','files') for s,h in result[kind].items() if expected.get(kind,{}).get(s)!=h]
            raise RuntimeError('manifest differs: '+repr(changed))
        print('PASS: '+str(len(files))+' local files and '+str(len(inputs))+' inherited inputs')
    else:
        target.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8')
        print('Sealed '+str(len(files))+' files and '+str(len(inputs))+' inputs')

if __name__=='__main__':
    main()
