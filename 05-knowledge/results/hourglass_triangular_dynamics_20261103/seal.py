"""Verify receipts and local links; pin sources and hash all packet artifacts."""
from pathlib import Path
from hashlib import sha256
import json
import os
import re
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
def digest(p):return sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
for a,b in [('normal.json','optimized.json'),('audit.normal.json','audit.optimized.json'),
            ('completions.json','completions.optimized.json')]:
    if (HERE/a).read_bytes()!=(HERE/b).read_bytes():raise ValueError(('receipt mismatch',a,b))
direct=json.loads((HERE/'completion_direct.json').read_text())
recur=json.loads((HERE/'completions.json').read_text())
if direct!={'full_strong':recur['full']['strong'],'quotient_strong':recur['line_quotient']['strong']}:
    raise ValueError('completion census mismatch')
links=0
for file in HERE.glob('*.md'):
    for link in re.findall(r'\]\(([^)]+)\)',file.read_text(encoding='utf-8')):
        if '://' in link or link.startswith('#'):continue
        target=Path(os.path.normpath(str(file.parent/link.split('#')[0])))
        if not target.exists():raise ValueError(('missing link',file.name,link,str(target)))
        links+=1
deps=[
 '05-knowledge/hypotheses/HYP-3053-tournament-diagonal-layer-flow.md',
 '04-computation/tournament_diagonal_layer_flow_codex_s217.py',
 '01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md',
 '05-knowledge/results/arithmetic_braids_20260917_geometry.md',
 '05-knowledge/results/arithmetic_seams_20260921_dynamics.md',
 '05-knowledge/results/session146_carrier_primegame_20261103/prime290/README.md',
 '05-knowledge/results/session146_carrier_primegame_20261103/prime290/WHEEL_GENERAL.md',
 '05-knowledge/results/session146_carrier_primegame_20261103/machines/README.md',
 '05-knowledge/results/positive_source_jacobian_prime_clocks_20261103/ARITHMETIC_CLOCKS.md',
 '05-knowledge/results/primegame_port_spectrum_20261103/verify.py',
 '05-knowledge/results/primegame_port_spectrum_20261103/port_scout.py',
 '05-knowledge/results/regular_parent_symmetry_20261029/verify.py',
 '05-knowledge/results/hourglass_quadratic_ports_20261103/README.md',
 '05-knowledge/results/hourglass_quadratic_ports_20261103/NATIVE.md',
 '05-knowledge/results/hourglass_quadratic_ports_20261103/ARITHMETIC.md',
 '05-knowledge/results/hourglass_quadratic_ports_20261103/COMPLETIONS.md',
 '05-knowledge/results/residue_290_transport_20261102/README.md',
 '05-knowledge/results/quadratic_2916_odd_pair_triangles_n4_20261102.md',
 '05-knowledge/results/session146_carrier_primegame_20261103/carrier9/PROOF.md',
 '05-knowledge/results/session146_carrier_primegame_20261103/carrier10/PROOF.md',
 '05-knowledge/results/session146_carrier_primegame_20261103/carrier10_upper/PROOF.md',
]
(HERE/'SOURCE_PINS.json').write_text(json.dumps({p:digest(ROOT/p) for p in deps},indent=2,sort_keys=True)+'\n',encoding='utf-8')
files={p.name:digest(p) for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.json'}
(HERE/'MANIFEST.json').write_text(json.dumps({'hash_basis':'CRLF normalized to LF',
   'artifacts':files,'local_links_checked':links,'source_count':len(deps)},indent=2,sort_keys=True)+'\n',encoding='utf-8')
print(json.dumps({'artifacts':len(files),'sources':len(deps),'local_links':links,'normal_O_match':True},sort_keys=True))
