"""Finite certificates used in the all-length arithmetic obstruction proofs."""
from common import *
from collections import Counter
from math import prod

def det(M):
    return M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])-M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])+M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0])

def run():
    barriers=[]
    for r in read('even_valuations.json.gz')['records']:
        if r['type']=='bounded_constant':continue
        g=r['odd_content_valuation'];f={int(b):c for b,c in r['normalized_odd_modes'].items()}
        if r['type']=='bounded_parity':
            p=next(p for p,j in enumerate(r['eventual_valuations_by_parity']) if j>=g+2);period=2;residue=p;j=r['eventual_valuations_by_parity'][p]
        else:
            p=r['unbounded_parities'][0];t=next(t for t in [0,1] if sum(c*pow(b,2*t+p,16) for b,c in f.items())%16==8)
            period=4;residue=2*t+p;j=g+3
        first=residue
        while first<j+1:first+=period
        for k in range(64):need(v2(value(r['exponential_coefficients'],first+period*k))==j,'literal exact-valuation AP control')
        bad=[b for b,c in enumerate(r['exponential_coefficients'],1) if b%2 and c and v2(c)<j-1]
        need(bad,'nonintegral odd-mode coefficient in lowest head')
        barriers.append({'birth_id':r['birth_id'],'port':r['port'],'type':r['type'],'odd_content_valuation':g,'exact_valuation':j,'a_residue':residue,'a_period':period,'first_control_a':first,'nonintegral_odd_bases':bad})
    need(len(barriers)==126,'all60 parity and66 unbounded families')
    sys.path.insert(0,str(HERE.parent/'inverse65_zeckendorf_20261103'))
    from bank import MATS
    determinants={name:det(M) for name,M in MATS.items()}
    need(all(abs(d)==1 or abs(d)>=5 for d in determinants.values()) and all(abs(d)!=6 for d in determinants.values()),'determinant6 absent from the generated multiplicative monoid')
    smooth=[]
    for name,d in determinants.items():
        d=abs(d);a=0;b=0
        while d%2==0:d//=2;a+=1
        while d%3==0:d//=3;b+=1
        if d==1:
            need(b>=3*a,'2,3-smooth determinant cone');smooth.append({'letter':name,'v2':a,'v3':b})
    periodic=[];byid={(r['birth_id'],r['port']):r for r in barriers}
    for r in read('even_valuations.json.gz')['records']:
        e=r['exponential_coefficients'];j=r['eventual_valuations_by_parity'][0] if r['type']=='bounded_constant' else byid[r['birth_id'],r['port']]['exact_valuation'];d=j-1
        support=[b for b,c in enumerate(e,1) if c+(2**d if b==1 else 0)]
        reason='scalar recurrence rank exceeds3' if len(support)>3 else 'determinant cone excludes6^L' if support==[1,2,3] else None
        if reason:periodic.append({'birth_id':r['birth_id'],'port':r['port'],'lowest_head_division':d,'mode_support':support,'obstruction':reason})
    need(len(periodic)==470,'all but eight pure3 families in periodic-word barrier')
    out={'scope':'Exact source coefficients for the affine TT-chart odd-mode obstruction and the periodic65-seam rank/determinant obstruction, including every arithmetic-progression reparameterization. Neither excludes general recursive tournament construction.','affine_TT_barrier_families':len(barriers),'barrier_type_histogram':dict(Counter(r['type'] for r in barriers)),'exact_valuation_AP_controls':len(barriers)*64,'records':barriers,'periodic65_determinants':determinants,'periodic65_absolute_determinant_histogram':dict(Counter(abs(x) for x in determinants.values())),'smooth_determinant_cone':smooth,'periodic_barrier_families':len(periodic),'periodic_barrier_histogram':dict(Counter(r['obstruction'] for r in periodic)),'periodic_barriers':periodic}
    save('repair_obstructions.json.gz',out);print(json.dumps({k:v for k,v in out.items() if k not in ['records','periodic65_determinants','periodic_barriers','smooth_determinant_cone']}),flush=True)
if __name__=='__main__':run()
