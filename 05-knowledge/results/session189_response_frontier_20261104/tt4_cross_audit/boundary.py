#!/usr/bin/env python3
from pathlib import Path
import json,hashlib,subprocess
from algebra import tt,need
from audit import fullbanks,coef
P=Path(__file__).resolve().parent

def main():
 source=tt([0,1,2,3])+[0]
 for i in range(4):source[i]|=1<<4
 targets=[(1,0,3,2),(1,3,0,2),(2,0,3,1),(2,3,0,1)];combos=[[(0,1,3,1),(0,1,2,-1)],[(0,1,3,2),(0,1,2,-1)],[(0,1,3,2),(0,1,2,3)],[(0,2,3,1),(0,2,1,-1)]];graphs=[source]
 for order in targets:graphs.append(tt(order)+[15])
 bank=fullbanks(graphs,'q2_boundary_banks');outcomes=[]
 for j,terms in enumerate(combos):
  value=sum(w*(coef(bank[0],p,q,r,4)-coef(bank[j+1],p,q,r,4))for p,q,r,w in terms);outcomes.append(value)
 need(outcomes==[4,4,-4,2],'quotient order two positive boundary')
 result={'status':'PASS','Q':[2,0],'source':source,'targets':[b['rows']for b in bank[1:]],'positive_sym_capacity':2,'constant_gammas':[2,2,-2,1],'actual_grouped_quartic_gaps':outcomes,'full_named_graphs':len(bank),'full_named_coefficients':sum(len(b['coefficients'])for b in bank),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'bank_sha256':hashlib.sha256((P/'q2_boundary_banks.json').read_bytes()).hexdigest()};(P/'boundary_receipt.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(json.dumps(result,sort_keys=True))
if __name__=='__main__':main()
