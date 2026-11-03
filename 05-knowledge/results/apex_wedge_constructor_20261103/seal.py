"""Pin inherited mechanisms, or seal this packet after successful reproduction."""
from pathlib import Path
import hashlib,json,sys
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SOURCES=[
'feedback_word_hamilton_spectrum_20261020/proofs.md',
'session134_constructive_spectrum_frontier_20261101/two_apex/THEORY.md',
'hourglass_triangular_dynamics_20261103/HOURGLASS.md',
'ternary_repair_clock_20261103/native.py',
'overnight_research_20261029/TWO_PORTS.md',
'overnight_research_20261029/MARKED_CLOCK_PORTS.md',
'mixed_returns_20261103/controllers.json',
'inverse65_zeckendorf_20261103/ZECKENDORF.md',
'session133_fractran_spectrum_refinement_20261101/path_counter/THEORY.md',
'two_port_return_induction_20261103/NATIVE_RESET_INDUCTION.md',
'session157_gaps_20261103/universal/NOTE.md']
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
if '--pins' in sys.argv:
    output=HERE/'source_pins.json'
    if output.exists():raise RuntimeError('Do not silently repin inherited sources')
    pins={'session_base':'3d6628b1804c','sha256':{
        '05-knowledge/results/'+name:digest(ROOT/'05-knowledge/results'/name) for name in SOURCES}}
    output.write_text(json.dumps(pins,indent=2)+'\n',encoding='utf8')
else:
    if not (HERE/'reproduction.json').exists():raise RuntimeError('Reproduce before sealing')
    data={p.name:digest(p) for p in sorted(HERE.iterdir()) if p.is_file() and p.name!='manifest.json'}
    (HERE/'manifest.json').write_text(json.dumps({'sha256':data},indent=2)+'\n',encoding='utf8')
print('PASS: pinned sources' if '--pins' in sys.argv else 'PASS: sealed packet')
