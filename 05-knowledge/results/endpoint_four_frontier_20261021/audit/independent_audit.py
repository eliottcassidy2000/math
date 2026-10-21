#!/usr/bin/env python3
"""Independent guarded affine-family, printed-word and literal-hostile audit."""
from pathlib import Path
from math import gcd
from collections import Counter
import re,json,sys,hashlib
from audit_core import *
HERE=Path(__file__).resolve().parent
PACKET=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else HERE.parent

def printed_words(text,n,ell=None):
    formulas=re.findall(r'Q[123]=\(([^)]+)\)',text)
    need(len(formulas)==6,'six labelled printed theorem words')
    values={'a':[0],'b':[1],'x':[2],'y':[n-3],'c':[n-2],'z':[n-1]}
    if ell is not None:
        values['D']=list(range(3,ell+3));values['U']=list(range(ell+3,n-1));formulas=formulas[:3]
    else:values['D']=list(range(3,n-3));formulas=formulas[3:]
    return [sum((values[s.strip()] for s in f.split(',')),[]) for f in formulas]

def guards(n,h,ell):
    return ell>0 and ell%2 and n>=max(15,2*(ell+3)) and h not in {2,ell+2} and (ell==1 or h not in {4,ell-1})

def odd_proof_controls(n,ell,text):
    q1,q2,q3=printed_words(text,n,ell);r1,r2,r3=map(reversal,[q1,q2,q3]);L=ell+3
    need(r1==reversal(list(range(1,n))+[0]),'cyclic-cut reversal identity')
    if ell==1:need(r2[0]==n-1,'singleton-prefix pairs a,z')
    else:
        need(r2[4]==n-1 and r2[0]==n-ell,'actual second-word images')
        need(ell+3<=n-ell<=n-2 and n-ell not in {0,2,n-2},'n-ell actual U position')
    q=q3;inv={v:i for i,v in enumerate(q)}
    need(q==list(range(3,L))+[0,1,2]+list(range(L,n)),'printed final position rotation')
    g={i:r3[n-1-i] for i in range(n)}
    need(g=={i:q[n-1-inv[n-1-i]] for i in range(n)},'commutator expression')
    cycles=permutation_cycles(g);actual=sorted(map(len,cycles));expected=sorted([L//gcd(L,3)]*(2*gcd(L,3))+[1]*(n-2*L))
    need(actual==expected and L//gcd(L,3)%2==0,'exact disjoint block cycles')

def main():
    text=(PACKET/'README.md').read_text();counts=Counter()
    for n in range(15,51):
        for ell in range(1,n//2-2,2):
            qs=printed_words(text,n,ell)
            for h in range(2,n-3):
                if not guards(n,h,ell):continue
                A=symbolic_source(n,{h,n-3,n-2},[(ell+3,2,1)])
                replay_family(A,qs);counts['odd_prefix_affine_cells']+=1
                S=k_move(A,qs[0]);k=n-1-h
                need(all(S[v][n-1]==int(v in {0,n-2,2,k}) for v in range(n-1)),'exact S(z) incoming row')
                if ell==1:
                    need(all(S[0][v] in {0,1} for v in range(1,n)) and sum(S[0])==5,'singleton-prefix score5')
                else:
                    need(S[4][0]==1 and S[n-ell][n-1]==0,'actual Q2 opposite incidences')
    for n in range(15,105):
        qs=printed_words(text,n);A=symbolic_source(n,{2,n-3,n-2},[(n-2,1,1),(2,n-2,1)])
        replay_family(A,qs);counts['complete_h_equals_x_affine_cells']+=1
        r1,r2=map(reversal,qs[:2]);need((r1[n-1],r1[0])==(2,n-2),'h=x first anti-pair')
        need((r2[1],r2[0])==(2,n-2),'h=x second anti-pair')
    for n in range(15,301):
        for ell in range(1,n//2-2,2):odd_proof_controls(n,ell,text);counts['exact_permutation_and_witness_controls']+=1
    for line in ['(a,x,3,4,5,6,...,c,z,b)','(3,4,5,a,x,6,...,c,z,b)','(3,4,5,a,b,x,6,...,c,z)']:need(line in text,'printed cut-six literal formula')
    hostiles=[];expected=json.loads((PACKET/'receipt-normal.json').read_text())['repaired_same_labelled_sources']
    frozen=json.loads((PACKET/'old-template-hostiles.json').read_text())
    need(len(frozen)==len(expected)==4 and [(w['n'],w['h']) for w in frozen]==[(15,3),(15,6),(18,3),(18,9)],'exact four hostile sources')
    for w,record in zip(frozen,expected):
        n,h=w['n'],w['h'];A=literal_rows(w['rows']);old=check_literal(A,{h,n-3,n-2},w['words'],[False,False,True])
        need(A[6][2]==1 and guards(n,h,3),'hostile cut-six guards')
        qs=printed_words(text,n,3);need(record['repair_words']==qs,'receipt differs from printed new words')
        new=check_literal(A,{h,n-3,n-2},qs,[False,False,False])
        hostiles.append({'n':n,'h':h,'old_blue_steps':old,'new_blue_steps':new,'cut':6,'same_literal_source_and_target':True})
    # Exact logical scope subtraction in formal branch coordinates, not a tournament census.
    for n in range(15,41):
        for h in range(2,n-3):
            for xc in (False,True):
                for xy in (False,True):
                    prior=(xc and h not in {2,n-4}) or (not xc and xy and h not in {2,n-4})
                    added=(h==2 and xc)
                    remainder=(h==2 and not xc) or h==n-4 or (not xc and not xy)
                    need((not (prior or added))==bool(remainder),'remaining-frontier branch logic')
    need('must fail every applicable odd-prefix guard' in text,'missing failed-guard frontier condition')
    print(json.dumps({'status':'PASS','method':'independent GF(2) elimination over complete affine source families; no producer imports','odd_prefix_all_h_orders':[15,50],'h_equals_x_orders':[15,104],'permutation_and_direct_witness_orders':[15,300],'counts':dict(counts),'same_source_hostile_repairs':hostiles,'scope':'Guarded odd prefixes and the closed h=x,x->c branch. Remaining sources must also fail every applicable odd-prefix guard.'},sort_keys=True,indent=2))
if __name__=='__main__':main()
