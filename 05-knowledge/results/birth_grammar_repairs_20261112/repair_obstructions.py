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
    out={'scope':'Exact source coefficients for the affine TT-chart odd-mode obstruction and the periodic65-seam determinant obstruction. Neither excludes general recursive tournament construction.','affine_TT_barrier_families':len(barriers),'barrier_type_histogram':dict(Counter(r['type'] for r in barriers)),'exact_valuation_AP_controls':len(barriers)*64,'records':barriers,'periodic65_determinants':determinants,'periodic65_absolute_determinant_histogram':dict(Counter(abs(x) for x in determinants.values()))}
    save('repair_obstructions.json.gz',out);print(json.dumps({k:v for k,v in out.items() if k not in ['records','periodic65_determinants']}),flush=True)
if __name__=='__main__':run()
