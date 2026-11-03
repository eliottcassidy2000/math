"""Pin consulted sources and seal LF-normalized text artifacts."""
from pathlib import Path
import hashlib,json,subprocess
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SOURCES=[
 'AGENTS.md',
 '05-knowledge/results/fibonacci_parity_precision_20261103/COMPLEMENTS.md',
 '05-knowledge/results/overnight_rank_return_20261029/DEFECT_TRANSFER.md',
 '05-knowledge/results/tournament_future_boundary_20261020/README.md',
 '05-knowledge/results/charged_recurrence_memory_20261028/COMPLEMENT.md',
 '05-knowledge/results/metagraph_primegame_20261102/HANDOFF.md',
 '05-knowledge/results/metagraph_primegame_20261102/PADDED_EARS.json',
 '05-knowledge/results/square_birth_windows_20261103/GROWTH.md',
 '05-knowledge/results/session150_hamiltonian_spectrum_20261103/transfer/TRANSFER.md',
 '05-knowledge/results/primegame_port_spectrum_20261103/FOREST_RESPONSE_RECURRENCE.md',
 '05-knowledge/results/primegame_port_spectrum_20261103/verify.py',
 '05-knowledge/results/primegame_port_spectrum_20261103/port_scout.py',
 '05-knowledge/results/primegame_port_spectrum_20261103/growth_identity.py',
 '05-knowledge/results/regular_parent_symmetry_20261029/verify.py',
 '05-knowledge/results/positive_source_jacobian_prime_clocks_20261103/PARITY_CARRIER.md',
 '05-knowledge/results/two_selector_boolean_20261102/CARRIER.md',
 '01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md',
 '01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md'
]
def record(path):
    raw=path.read_bytes().replace(b'\r\n',b'\n')
    return {'LF_normalized_sha256':hashlib.sha256(raw).hexdigest(),'LF_normalized_bytes':len(raw)}
def main():
    normal=(HERE/'verification.json').read_bytes()
    optimized=(HERE/'verification_optimized.json').read_bytes()
    if normal!=optimized:raise RuntimeError('normal and -O receipts differ')
    sources={p:record(ROOT/p) for p in SOURCES}
    pinpath=HERE/'source_pins.json'
    if pinpath.exists():
        prior=json.loads(pinpath.read_text())
        if prior['sources']!=sources:raise RuntimeError('inherited source changed: review before resealing')
    else:
        head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
        pinpath.write_text(json.dumps({'snapshot_head':head,'sources':sources},indent=2)+'\n')
    files={p.name:record(p) for p in sorted(HERE.iterdir()) if p.is_file() and p.name!='manifest.json'}
    (HERE/'manifest.json').write_text(json.dumps({'hash_basis':'UTF-8 text with CRLF normalized to LF','files':files},indent=2)+'\n')
    print(json.dumps({'status':'PASS','source_pins':len(sources),'sealed_files':len(files),'normal_optimized_identical':True}))
if __name__=='__main__':main()
