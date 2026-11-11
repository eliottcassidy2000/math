"""Literal full-state controls for the parameterized H13 difference."""
from pathlib import Path
import importlib.util,json

def main():
    path=Path(__file__).resolve().parents[1]/'modular_cycles/probe.py'
    spec=importlib.util.spec_from_file_location('native',path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    controls=0
    for A in range(8):
        for B in range(4):
            states=[]
            for cut in (6,3):
                rows=m.c.construct((18,8,3,5,14),cut,A,B)
                states.append(m.state(rows,len(rows)-1))
            b=B.bit_length();r=int(f'{B:0{b}b}'[::-1],2) if B else 0
            side=(13 if B else 7,(1<<b)*(104*A+45)-13*r,13)
            delta=3*A-(A&1)+2
            if states[0][1:]!=side or states[1][1:]!=side or states[1][0]-states[0][0]!=delta:
                raise RuntimeError((A,B,states,side,delta))
            controls+=1
    prefix_controls=0
    for A in (0,1):
        for exponent in (1,2):
            for cut in (6,3):
                original=m.c.construct((18,8,3,5,14),cut,A,0)
                st=m.state(original,len(original)-1)
                rows=list(original)
                for _ in range(exponent):
                    # C3 -> old source, with the retained root shifted by3.
                    oldsize=len(rows)
                    rows=[r|(((1<<oldsize)-1)<<3) for r in (2,4,1)]+[r<<3 for r in rows]
                actual=m.state(rows,len(rows)-1)
                if actual!=tuple((3**exponent)*v for v in st):raise RuntimeError(('prefix',A,exponent,cut,st,actual))
                prefix_controls+=1
    st=(1,1,1,1);seq=[]
    for n in range(40):
        h,e,u,d=st;seq.append(h);st=(h+u-e,h-e,h+u,h)
    if seq[:3]!=[1,1,3] or any(seq[n+3]!=sum(seq[n:n+3]) for n in range(37)):
        raise RuntimeError(('tribonacci',seq))
    report={'status':'FINITE-EXACT controls on all-address proof in PORTABLE_GAPS.md',
            'literal_pair_controls':controls,'universe':'A0..7,B0..3, both specified H13 cuts, every rooted coordinate by subset DP',
            'J_singleton_counts_first_16':seq[:16],'J_recurrence_controls':37,
            'literal_ordinal_prefix_controls':prefix_controls,
            'positive_controls':'A0 unit width; A1 width4; nonempty left/right flanks',
            'boundary':'multiples of6 use an actual C3 ordinal prefix; absolute baseline compatibility remains open'}
    Path(__file__).with_name('portable_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
