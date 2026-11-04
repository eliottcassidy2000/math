"""Exact finite certificate in the proved nine-dimensional module spectrum.

This is not an unconstrained sequence fit: the two-bin module-run bound is
n_exterior+bins=7+2=9, proved in RESPONSE.md, and supplies the finite space.
"""
from pathlib import Path
from fractions import Fraction as Q
import json
from repair import dial,native,need,endpoint
HERE=Path(__file__).resolve().parent

def main():
    records=json.loads((HERE/'repair.json').read_text())['sharp_entry_bridges']
    samples={r['m']:[value for pair in r['pairs'] for value in pair['fields']] for r in records}
    bulk_transport=[]
    for m in range(9):
        T=tuple(records[m]['target_rows']);U=native.extend(T,{0,1,2})
        # Keep every old physical name; only the canonical role indices move.
        order=list(range(m+4))+[len(T)]+list(range(m+4,m+7))
        relabelled=[sum(1<<j for j,b in enumerate(order) if U[a]>>b&1) for a in order]
        need(relabelled==records[m+1]['target_rows'],'Actual bulk birth with full role transport')
        bulk_transport.append(dict(m=m,new_out_mask=7,canonical_to_actual_indices=order))
    A=[[Q(j**m) for j in range(1,10)]+list(map(Q,samples[m])) for m in range(1,10)]
    for c in range(9):
        p=next(r for r in range(c,9) if A[r][c]);A[c],A[p]=A[p],A[c]
        d=A[c][c];A[c]=[x/d for x in A[c]]
        for r in range(9):
            if r!=c:
                d=A[r][c];A[r]=[x-d*y for x,y in zip(A[r],A[c])]
    coefficients=[[A[j][9+i] for j in range(9)] for i in range(48)]
    need(all(x.denominator==1 for row in coefficients for x in row),'Integral power coefficients')
    coefficients=[[int(x) for x in row] for row in coefficients]
    def evaluate(m):return [sum(c*j**m for j,c in enumerate(row,1)) for row in coefficients]
    need(all(evaluate(m)==samples[m] for m in range(1,10)),'All independent DP samples')
    # An actual larger graph counted by a different algorithm, not used to fit.
    from next_returns import gap
    D=(2,4,9,3)
    for _ in range(20):D=native.extend(D,{0,1,2})
    X=dial(D,3);fb=(0,3,len(D));actual=[]
    for v in range(3):actual+=gap.pair_fields(X,fb,v,3)
    need(actual==evaluate(20),'Independent larger actual gap-count holdout')
    holdout_order=len(X)
    need([evaluate(1)[i] for i in (5,21,37)]==[411,393,183],'Three named F fields')
    for m in range(1,10):
        v=evaluate(m)
        need([v[i] for i in (5,21,37)]==[76*3**m+96*2**m-9,76*3**m+96*2**m-27,96*2**m-9],'Closed three-mark responses')
        need(all(v[i]%3==0 for i in (0,5,21,37)),'Ternary-divisible wrapper boundary')
    active=sorted({j for row in coefficients for j,c in enumerate(row,1) if c})
    need(all(coefficients[16*p+i][2]==0 for p in range(3) for i in (0,2,6,7,11)),'All four golden marks omit the ternary mode')
    need(len({tuple(x%7 for x in evaluate(m)) for m in range(1,7)})==6,'Six distinct realized modulo7 responses')
    C=[row[:3] for row in coefficients]
    need(active==[1,2,3] and all(row[1]%6==0 and row[2]%38==0 for row in C),'Integral response lattice containment')
    need(C[3]==[1,0,0] and [C[2][j]-5*C[3][j] for j in range(3)]==[0,6,0],'First two lattice generators')
    need([C[1][j]-3*C[5][j]-16*C[2][j]+107*C[3][j] for j in range(3)]==[0,0,38],'Third lattice generator')
    # Realization guards: on this orbit only, a hidden 19-coordinate is free.
    need(pow(2,13,19)==3 and len({pow(2,m,19) for m in range(1,19)})==18,'Actual modulo19 orbit')
    for m in range(1,37):
        need(pow(pow(2,m,19),13,19)==pow(3,m,19),'Realized third coordinate recovered modulo19')
        need(pow(3,m,7)==((-1)**m*pow(2,m,7)**2)%7,'Necessary parity phase modulo7')
    # Two literal sources have equal golden field charts but unequal continuations.
    hostile=[]
    from check import phase
    for m in (1,4):
        r=records[m];X=tuple(r['target_rows']);charts={}
        for pair in r['pairs']:
            v=pair['fields']
            for p,off in ((pair['ports'][0],2),(3,7)):
                golden=[v[0],(v[off+4]-v[off])//2,(v[off+4]+v[off])//2]
                charts[str(p)]=phase(golden)
        U=native.extend(X,set(range(len(X)))-{0});H=endpoint(U,0,3)[0]
        need(H==samples[m][0]+samples[m][5]-1,'Literal one-birth continuation')
        hostile.append(dict(m=m,rows=list(X),charts=charts,E2=samples[m][5],birth_out_mask=U[-1],birth_H=H,birth_H_mod7=H%7))
    need(hostile[0]['charts']==hostile[1]['charts'] and hostile[0]['birth_H_mod7']!=hostile[1]['birth_H_mod7'],'Both field charts miss a legal next-step parity bit')
    out=dict(status='EXACT_CERTIFICATE_IN_PROVED_FINITE_RESPONSE_SPACE',exterior_vertices=7,bins=2,run_bound=9,
             sample_m=list(range(1,10)),ports=[[0,3],[1,3],[2,3]],fields_per_pair=16,
             field_order=['H','P2','E1_first','S1_first','q_first','E2_first','S2_first','E1_second','S1_second','q_second','E2_second','S2_second','W2_00','W2_01','W2_10','W2_11'],
             bases=list(range(1,10)),coefficients=coefficients,active_bases=active,
             actual_bulk_transport=bulk_transport,
             integral_row_lattice_diagonal=[1,6,38],smith_diagonal=[1,2,114],
             boundary_m0=samples[0],positive_formula_at_m0=evaluate(0),
             holdout=dict(m=20,order=holdout_order,feedback=list(fb),actual_fields=actual),
             phase_mod7_period_divides=6,realized_mod19_period=18,phase_hostile=hostile)
    (HERE/'response_recurrence.json').write_text(json.dumps(out,indent=2)+'\n')
    print('48 joint fields; proved space dimension9; active bases',active,'independent holdout order',holdout_order,'phase hostile',[z['birth_H'] for z in hostile])

if __name__=='__main__':main()
