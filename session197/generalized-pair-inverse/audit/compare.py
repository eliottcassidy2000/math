"""Compare independent output with the reviewed source, only after generation."""
from pathlib import Path
import hashlib,json,sys
if len(sys.argv)!=4:
 raise SystemExit('usage: compare.py SOURCE_DIRECTORY INDEPENDENT_BANK_FILE INDEPENDENT_RECEIPT')
source=Path(sys.argv[1]).resolve();bank=Path(sys.argv[2]);receipt=Path(sys.argv[3])
sys.path.insert(0,str(source))
from inverse import source_bank,candidate_sources
expected='\n'.join(['C1 '+' '.join(map(str,s))for s in sorted(source_bank(1861))]+['SRC '+' '.join(map(str,s))for s in sorted(candidate_sources(1861))])+'\n'
actual=bank.read_text()
if actual!=expected:raise RuntimeError('complete source-bank mismatch')
ind=json.loads(receipt.read_text());support=json.loads((source/'SUPPORT.json').read_text())
if ind['support']!=sorted(map(int,support['family_support'])):raise RuntimeError('complete support mismatch')
print(json.dumps(dict(status='PASS',c1_quadruples=23432,candidate_sources=2000,bank_sha256=hashlib.sha256(actual.encode()).hexdigest(),support_values=883),indent=2,sort_keys=True))
