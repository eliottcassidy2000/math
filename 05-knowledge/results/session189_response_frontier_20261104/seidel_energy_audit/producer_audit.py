#!/usr/bin/env python3
"""Independent full-bank comparison plus clean replay of the frozen producer."""
from pathlib import Path
from fractions import Fraction
import json,hashlib,subprocess,tempfile,shutil,sys
import audit as A
P=Path(__file__).resolve().parent

def main():
    source=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else P.parent/'session189-seidel-energy'
    names=['PROOF.md','verify.py','rational_banks.json','small_energy_bank.json','carousel_controls.json','normal.json','optimized.json']
    pins={f:hashlib.sha256((source/f).read_bytes()).hexdigest() for f in names};A.require((source/'normal.json').read_bytes()==(source/'optimized.json').read_bytes(),'producer normal optimized identity')
    raw=json.loads((source/'rational_banks.json').read_text());positions=0;quartics=0
    for record in raw['records']:
        q=tuple(record['rows']);e,o,E,_,_=A.banks(q,False);ec,oc,_,_,_=A.banks(A.opposite(q),False);bs=record['half_masks'];A.require(set(bs)==set(A.basis(len(q))),'complete block basis');A.require(E==Fraction(*record['energy']),'producer energy')
        for idx,(a,c) in enumerate((a,c) for a in bs for c in bs):
            A.require(e[a,c]+o[a,c]==Fraction(*record['source_coefficients'][idx]),'producer source coefficient');A.require(ec[a,c]+oc[a,c]==Fraction(*record['converse_coefficients'][idx]),'producer target coefficient');A.require(o[a,c]==Fraction(*record['common_odd_subset_coefficients'][idx]),'producer odd coefficient');positions+=1
            if a.bit_count()+c.bit_count()==4:quartics+=1
    small=json.loads((source/'small_energy_bank.json').read_text())
    for record in small:
        q=tuple(record['rows']);e,o,E,_,_=A.banks(q,False);ec,oc,_,_,_=A.banks(A.opposite(q),False);A.require(E==Fraction(*record['energy']),'small energy')
        for i in range(len(q)):
            for j in range(len(q)):
                if i!=j:A.require(e[(1<<i)|(1<<j),1<<i]+o[(1<<i)|(1<<j),1<<i]-ec[(1<<i)|(1<<j),1<<i]-oc[(1<<i)|(1<<j),1<<i]==Fraction(*record['shared_cubics'][i][j]),'small shared cubic')
    with tempfile.TemporaryDirectory(prefix='seidel-producer-review-') as td:
        D=Path(td)
        for f in names:shutil.copyfile(source/f,D/f)
        for opt,label in [([], 'normal'),(['-O'],'optimized')]:
            result=subprocess.check_output([sys.executable,*opt,str(D/'verify.py')]);A.require(result==(source/'normal.json').read_bytes(),'producer exact replay '+label);(P/('producer.'+label+'.json')).write_bytes(result)
            for f in ['rational_banks.json','small_energy_bank.json','carousel_controls.json']:A.require(hashlib.sha256((D/f).read_bytes()).hexdigest()==pins[f],'producer regenerated data '+f)
    A.require(pins=={f:hashlib.sha256((source/f).read_bytes()).hexdigest() for f in names},'producer sources changed during review')
    receipt={'status':'PASS','complete_named_producer_banks':len(raw['records']),'named_positions_per_side':positions,'quartic_positions_per_side':quartics,'small_labeled_quotients':len(small),'producer_source_pins':pins,'sources':{f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in ['producer_audit.py','audit.py']},'artifacts':{f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in ['producer.normal.json','producer.optimized.json']}}
    (P/'producer_receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':main()
