"""Exact five-mode interpolation, extra-size controls, and both field charts."""
from pathlib import Path
from fractions import Fraction
from math import comb,factorial
import json
from audit_moving import need,chart,bridge
HERE=Path(__file__).resolve().parent
def fit(y):
    A=[[Fraction(j**k) for j in range(5)]+[Fraction(y[k])] for k in range(5)]
    for i in range(5):
        p=next(j for j in range(i,5) if A[j][i]);A[i],A[p]=A[p],A[i];a=A[i][i];A[i]=[x/a for x in A[i]]
        for j in range(5):
            if j!=i:
                b=A[j][i];A[j]=[x-b*z for x,z in zip(A[j],A[i])]
    v=[r[-1] for r in A];need(all(x.denominator==1 for x in v),'integer spectral coefficients');return list(map(int,v))
def main():
    states=json.loads((HERE/'feedback_ports.json').read_text())['states'];source={x['m']:x for x in states if x['stage']=='source'};target={x['m']:x for x in states if x['stage']=='second_port'};series={}
    def obs(x):return [x['H'],*x['joint_pairs'][0]['port_data'][0]]
    for j,name in enumerate(['H','E1','S1','deletion_H','E2','S2']):
        delta=[obs(target[m])[j]-3*obs(source[m])[j] for m in range(1,8)];v=fit(delta)
        wanted=[0,0,-135276,135276,0] if name=='S2' else [0]*5
        need(v==wanted,'complete exact defect coefficients');need(all(sum(v[a]*(a**k) for a in range(5))==y for k,y in enumerate(delta)),'two extra sizes outside interpolation');series[name]=dict(coefficients=v,samples=delta)
    source_coefficients={name:fit([source[m]['golden'][j] for m in range(1,8)]) for j,name in enumerate(('H','X','Y'))}
    need(source_coefficients=={'H':[0,-1061,2486,0,0],'X':[0,-877,-74696,99409,0],'Y':[0,-876,-74696,99409,0]},'source modal coefficients')
    field=[];direction=chart([0,1,1])['F343_vector_codec'];qdir=bridge.five_charts([0,1,1])
    for m in range(1,8):
        k=m-1;d=67638*(3**k-2**k);g=source[m]['golden'];h=target[m]['golden'];need(h==[3*g[0],3*g[1]+d,3*g[2]+d],'exact golden affine transport')
        sc=source[m]['charts'];tc=target[m]['charts'];need(tc['F49_charge']==[(3*x+d*y)%7 for x,y in zip(sc['F49_charge'],(-3,1,-1))],'F49 and charge transport')
        A,rem=divmod(g[0]+1061,2486);need(rem==0,'dyadic realization guard');B,rem=divmod(g[1]+74696*A+877,99409);need(rem==0 and d==67638*(B-A),'ternary realization guard and recovered correction')
        need(d%7==(g[0]+2*g[1]+1)%7,'phase correction free from current golden vector')
        a,b,c=sc['F49_charge'];need(tc['F49_charge']==[(4*a+4*b+4)%7,(2*a+4*b+1)%7,(4*a+5)%7],'affine F49 pair transport with source-family charge')
        for a,b,e in zip(sc['F343_vector_codec'],tc['F343_vector_codec'],direction):need(list(b)==[(3*x+d*y)%7 for x,y in zip(a,e)],'all three F343 Frobenius frames')
        need((d%7==0)==(m%6==1),'same source clock, six-step correction period')
        for name,a in source[m]['five_Q4'].items():
            b=target[m]['five_Q4'][name];need(b['golden']==[3*x+d*y for x,y in zip(a['golden'],qdir[name]['golden'])],'all five actual Q4 linear transports')
        field.append(dict(m=m,delta=d,delta_mod7=d%7,source_codec=sc,target_codec=tc))
    # Integer factorization V=S diag(r!) B^t gives the exact sample lattice.
    S=[[0]*5 for _ in range(5)];S[0][0]=1
    for k in range(1,5):
        for r in range(1,5):S[k][r]=S[k-1][r-1]+r*S[k-1][r]
    need(all(j**k==sum(S[k][r]*factorial(r)*comb(j,r) for r in range(j+1)) for k in range(5) for j in range(5)),'unimodular Stirling/Pascal factorization')
    for m in range(1,8):
        need(obs(source[m])[4]-source[m]['H']==(0 if m==1 else 1),'source readiness boundary')
        need(obs(target[m])[4]-target[m]['H']==(0 if m==1 else 3),'target readiness boundary')
    result=dict(status='FINITE_EXACT_COEFFICIENTS_WITH_PROVED_FIVE_MODE_REDUCTION',m_range=[1,7],independent_interpolation_samples=5,extra_size_controls=2,
        spectral_roots=[0,1,2,3,4],annihilator_ascending=[0,24,-50,35,-10,1],sample_smith_factors=[1,1,2,6,24],sample_lattice_index=288,
        delta_response_coefficients=series,source_golden_coefficients=source_coefficients,field_transport=field,
        F49_affine_matrix=[[4,4],[2,4]],F49_affine_translation=[4,1],F49_affine_determinant_mod7=1,
        all_size_claim='Requires the ordered two-bin/two-feedback skeleton proof in TRANSPORT.md, not polynomial fitting alone.')
    (HERE/'feedback_transport.json').write_text(json.dumps(result,indent=2)+'\n');print('PASS exact3-ready transport, hidden3-mode, and characteristic7 phase')
if __name__=='__main__':main()
