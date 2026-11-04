#!/usr/bin/env python3
"""Recheck synthesis metrics, relative links, and small spectrum exclusion."""
from pathlib import Path
import json,sys,re,hashlib
P=Path(__file__).resolve().parent; R=P.parent
sys.path.insert(0,str(R/'clone_obstruction'))
from reviewer_audit import quotient,count

def need(x,s):
 if not x:raise ArithmeticError(s)
def main():
 note=(P/'NOTE.md').read_text();links=[]
 for s in re.findall(r'\]\(([^)]+)\)',note):
  if '://' not in s:
   need((P/s.split('#')[0]).exists(),'missing '+s);links.append(s)
 receipts=[json.loads((R/'a000568'/f).read_text()) for f in ('receipt.json','receipt-optimized.json')]
 for j,expected_ratio in zip(receipts,(1.028,0.954)):
  b=next(x for x in j['benchmarks'] if x['order']==200)
  need(round(b['median_speedup'],3)==expected_ratio,'timing ratio')
  for mode,states,transitions in [('projective',95709,2303050),('ungauged_fixed',96836,2770773)]:
   need(b['stats'][mode]['peak_states']==states and b['stats'][mode]['transitions']==transitions,'state metrics')
 memory=json.loads((R/'a000568/memory-receipt.json').read_text())
 need({x['projective']:x['traced_peak_bytes'] for x in memory['rows']}=={True:10061760,False:9719012},'allocation metrics')
 small={}
 for n in range(1,7):
  vals=sorted({count(quotient(n,k)) for k in range(1<<(n*(n-1)//2))})
  need(7 not in vals and 21 not in vals,'small global exclusions')
  small[n]={'labelled_tournaments':1<<(n*(n-1)//2),'H_values':vals}
 print(json.dumps({'status':'PASS','note_sha256':hashlib.sha256((P/'NOTE.md').read_bytes()).hexdigest(),'relative_links':links,'metrics':'PASS','small_spectra':small},indent=2,sort_keys=True))
if __name__=='__main__':main()
