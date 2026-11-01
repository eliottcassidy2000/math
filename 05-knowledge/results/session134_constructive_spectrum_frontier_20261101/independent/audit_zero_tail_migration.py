#!/usr/bin/env python3
"""Independent arbitrary-length zero-tail migration formula controls.
Uses only auditor graph/state functions; no producer import.
"""
import json
from audit_migration_compiler import state,literal,tournament,build

def need(c,msg):
    if not c:raise RuntimeError(msg)
def local_tail(r):
    # local0 old apex; fresh local1..r point to0 and increase chronologically.
    return [0]+[(1|sum(1<<j for j in range(i+1,r+1))) for i in range(1,r+1)]
def expected(old,r):
    h,F,q,g=old
    if r==1:return F,F+q+g,q+h,h-g
    m=r-1;a=2**m;b=2*3**m-2**(m+1)+1
    return F+(a-1)*q+g,F+(2*a-1)*q+g,a*h+(a-1)*F+b*q+(a-2)*g,F+(a-2)*q+g

def main():
    cases=literal_cases=0;hostile=None
    for n in range(1,4):
        for code in range(1<<(n*(n-1)//2)):
            A=tournament(n,code)
            for cut in range(1<<n):
                old=state(A,cut)
                for r in range(1,7):
                    B,S=build(A,cut,local_tail(r),1);got=state(B,S)
                    need(got==expected(old,r),f'zero-tail formula mismatch {n,code,cut,r,old,got}')
                    if r<=3:
                        need(literal(B,S)==got,'literal arbitrary-tail control');literal_cases+=1
                    if (n,code,cut,r)==(1,0,1,3):
                        h,F,q,g=old;a=4
                        wrong=a*h+(a-1)*F+9*q+(a-2)*g
                        need(got[2]==20 and wrong==18,'r3 coefficient hostile')
                        hostile={'old_h_F_q_g':old,'r':r,'actual_new_h_F_q_g':got,'wrong_q_if_coefficient_9':wrong,'correct_q_coefficient':11}
                    cases+=1
    print(json.dumps({'status':'PASS','old_marked_parent_cases':74,'tail_lengths':[1,6],'actual_graph_checks':cases,'literal_permutation_cover_checks':literal_cases,'discarded_r3_coefficient_hostile':hostile,'scope':'All labelled old parents orders1..3/all cuts, tails1..6. All-order formula requires the stated local cover enumeration.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
