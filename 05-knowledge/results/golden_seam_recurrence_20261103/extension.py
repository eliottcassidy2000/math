"""Compile all strong five-vertex marked atoms, with exact seam guards."""
import json
from pathlib import Path
from transfer import (STATES, PROJECTION, POSITIVE, COST, chart, complement,
                      reversed_rows, advance, need)
from search import explore, primes

HERE=Path(__file__).resolve().parent


def strong(rows):
    n=len(rows)
    for u in range(n):
        seen=1<<u
        while True:
            more=seen
            for v,row in enumerate(rows):
                if seen>>v&1: more|=row
            if more==seen: break
            seen=more
        if seen!=(1<<n)-1: return False
    return True


def compile_atom(atom):
    R=reversed_rows(complement(atom))
    M=[advance(state,R) for state in STATES]
    Z=PROJECTION
    W=[[sum(M[i][k]*Z[k][j] for k in range(6)) for j in range(3)] for i in range(6)]
    need(all(v%2==0 for v in W[0]), 'integral first projection row')
    need(all((W[3][j]-W[0][j])%2==0 for j in range(3)), 'integral third projection row')
    p0=tuple(v//2 for v in W[0])
    p2=tuple((W[3][j]-W[0][j])//2 for j in range(3))
    p1=tuple(2*p0[j]+p2[j]-W[1][j] for j in range(3))
    P=(p0,p1,p2)
    need(all(v>=0 for row in P for v in row),'positive compiler image')
    need(all(sum(Z[i][k]*P[k][j] for k in range(3))==W[i][j]
             for i in range(6) for j in range(3)), 'all-word intertwining')
    initial=advance((),R,False)
    need(all(sum(initial[k]*Z[k][j] for k in range(6))==2*p0[j] for j in range(3)),
         'all-word initialization')
    return P


def atlas():
    return {f'r{mask}':{'source_rows':chart(5,mask),'matrix':compile_atom(chart(5,mask))}
            for mask in range(64) if strong(chart(5,mask))}


def main():
    bank=atlas()
    need(len(bank)==50,'all strong five-vertex charts')
    all_mats=POSITIVE|{k:v['matrix'] for k,v in bank.items()}
    all_costs=COST|{k:5 for k in bank}
    data={'atlas':bank,'distinct_five_matrices':len(set(v['matrix'] for v in bank.values()))}
    for name,mats,costs,bound in [
            ('golden_plus_r13',POSITIVE|{'r13':bank['r13']['matrix']},COST|{'r13':5},100000),
            ('all_strong_through_five',all_mats,all_costs,10000)]:
        d=explore(bound,mats,costs,8)
        present={r['H'] for r in d['values']}
        p=primes(bound); p7=primes(bound//7)
        d['prime_lane']={'eligible':len([x for x in p if x!=7]),'missing':[x for x in p if x!=7 and x not in present]}
        d['seven_prime_lane']={'eligible':len([x for x in p7 if x!=3]),'missing_p':[x for x in p7 if x!=3 and 7*x not in present]}
        data[name]=d
        print(json.dumps({'family':name,'bound':bound,'states':d['states'],'values':len(present),
                          'missing_prime_count':len(d['prime_lane']['missing']),
                          'first_missing_primes':d['prime_lane']['missing'][:12],
                          'missing_7p_count':len(d['seven_prime_lane']['missing_p']),
                          'first_missing_7p_primes':d['seven_prime_lane']['missing_p'][:12]}))
    (HERE/'extensions.json').write_text(json.dumps(data,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':main()
