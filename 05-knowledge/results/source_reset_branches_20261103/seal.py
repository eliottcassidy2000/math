"""Pin this research packet and its exact inherited inputs; validate local links."""
from pathlib import Path
import argparse,hashlib,json,re,os
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
INPUTS=(
 '05-knowledge/results/golden_inverse_defects_20261103/inverse.py',
 '05-knowledge/results/golden_seam_recurrence_20261103/atom_atlas.json',
 '05-knowledge/results/ternary_repair_clock_20261103/native.py',
 '05-knowledge/results/regular_parent_symmetry_20261029/verify.py',
 '05-knowledge/results/golden_hamilton_transfer_20261103/RECURRENCE.md',
 '05-knowledge/results/adaptive_metagraph_reformulation_20261103/ADAPTIVE.md',
 '05-knowledge/results/adaptive_metagraph_reformulation_20261103/ANCESTRY.md',
 '05-knowledge/results/adaptive_metagraph_reformulation_20261103/METAGRAPH.md',
 '05-knowledge/results/adaptive_state_induction_20261103/BLOCK_SWITCH.md',
 '05-knowledge/results/two_port_return_induction_20261103/ROOTED_RESETS.md',
 '05-knowledge/results/two_port_return_induction_20261103/NATIVE_RESET_INDUCTION.md',
 '05-knowledge/results/spectrum31_return_20261103/PROJECTIVE.md',
 '05-knowledge/results/spectrum31_return_20261103/TRANSPORT.md',
 '05-knowledge/results/session155_spectrum_20261103/primes/README.md',
 '05-knowledge/results/fermat_wieferich_order_20261014e.md',
 '05-knowledge/results/trijection_reset_20261020/reset_composition.md',
 '05-knowledge/results/session157_gaps_20261103/prime_power/NOTE.md',
 '05-knowledge/results/session157_gaps_20261103/strong3/NOTE.md',
 '05-knowledge/results/session157_gaps_20261103/universal/NOTE.md',
 '05-knowledge/results/session158_creative_20261103/reflection/NOTE.md',
 '05-knowledge/results/exclusion_7_21_20261103/INCOMING.md',
 '05-knowledge/results/primegame_port_spectrum_20261103/verify.py',
 '05-knowledge/results/primegame_port_spectrum_20261103/port_scout.py',
 '05-knowledge/results/primegame_port_spectrum_20261103/growth_identity.py',
)
def need(b,m):
 if not b:raise RuntimeError(m)
def fullpath(p):
 p=Path(os.path.abspath(str(p)))
 return Path('\\\\?\\'+str(p)) if os.name=='nt' and not str(p).startswith('\\\\?\\') else p
def read(p):return fullpath(p).read_text(encoding='utf-8').replace('\r\n','\n')
def sha(p):return hashlib.sha256(read(p).encode()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
 owned=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.json')
 links=0
 for p in owned:
  if p.suffix!='.md':continue
  for raw in re.findall(r'\]\(([^)]+)\)',read(p)):
   if re.match(r'^[a-zA-Z][a-zA-Z0-9+.-]*:',raw) or raw.startswith('#'):continue
   target=raw.split('#')[0]
   need(fullpath(p.parent/target).exists(),('broken relative link',p.name,raw));links+=1
 record=dict(hash_convention='SHA256 of UTF-8 text with CRLF normalized to LF',
             owned={p.name:sha(p) for p in owned},inputs={p:sha(ROOT/p) for p in INPUTS},checked_local_links=links)
 text=json.dumps(record,indent=2,sort_keys=True)+'\n';manifest=HERE/'MANIFEST.json'
 if args.write:manifest.write_text(text,encoding='utf-8',newline='\n')
 else:need(read(manifest)==text,'sealed packet or input changed')
 print(f'PASS {len(owned)} owned files, {len(INPUTS)} source pins, {links} local links')
if __name__=='__main__':main()
