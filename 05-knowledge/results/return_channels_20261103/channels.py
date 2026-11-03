"""Exact finite error transport, common ternary family, and parity obstruction."""
from collections import Counter
from itertools import permutations
from math import gcd,isqrt
from pathlib import Path
import json,sympy as sp
from common import *
from charges import modular_clock
T=(2,-1,0);R=(3,-1,-1);Q=(2,-2,1)
BRIDGE=('r43','a','r13')
SEEDS=(('r5','a','a'),('a','a','q1','a','a'),('q0','a','a','a','a','a'))
LOOPS={'L2':('r24','a','a','a'),'L6':('r24','q4'),
       'L45':('r19','r32','r48','q4'),'L13':('r19','r19','a','a','r18')}

def main():
    for l,m in MATS.items():
        need(m[0][0]%2==1 and m[1][0]%2==m[2][0]%2==0,('universal first-coordinate parity',l))
    need(step(T,matrix(BRIDGE))==tuple(2*x for x in R),'exact channel bridge')
    need(step(T,matrix(('r24','a','a')))==tuple(2*x for x in Q),'Q4 channel bridge')
    need(step(Q,MATS['a'])==tuple(-x for x in T),'return edge')
    states={T:()};gains={};census=[]
    for depth in range(1,6):
        nxt={}
        for v,w in states.items():
            for l,m in MATS.items():
                z=step(v,m)
                if min(z)>=0 or max(z)<=0:continue
                if z[2]==0 and z[0]==-2*z[1]:gains.setdefault(-z[1],w+(l,))
                nxt.setdefault(z,w+(l,))
        states=nxt;census.append({'length':depth,'mixed_error_states':len(states)})
    need(set(gains)=={-45,-36,-6,-2,-1,1,2,6,13,36,45},'complete five-letter gain set')
    loops={};X=sp.symbols('X')
    for name,w in LOOPS.items():
        p=sp.Matrix(matrix(w));z=step(T,matrix(w));gain=-z[1]
        need(z==tuple(gain*x for x in T),'loop gain')
        poly=p.charpoly(X).as_expr();quad=sp.Poly(sp.cancel(poly/(X-gain)),X)
        loops[name]={'word':w,'vertices':sum(COSTS[l] for l in w),'gain':gain,
                     'matrix':[[int(x) for x in row] for row in p.tolist()],'bulk_trace':-int(quad.nth(1)),'bulk_norm':int(quad.nth(0)),
                     'characteristic_polynomial':str(sp.factor(poly))}
    need([loops[x]['gain'] for x in LOOPS]==[-2,6,45,13],'four native loop gains')
    noncommuting=[]
    for order in (('L2','L6'),('L6','L2')):
        word=SEEDS[0]+LOOPS[order[0]]+LOOPS[order[1]]
        noncommuting.append({'macro_order':order,'word':word,'H':vector(word)[0],
                            'state':vector(word),'vertices':sum(COSTS[l] for l in word)})
    need(noncommuting[0]['H']!=noncommuting[1]['H'],'gain and clock do not determine bulk H')
    for entry in noncommuting:
        need(literal_frontier(actual_target(entry['word']))[0]==entry['H'],'literal noncommuting clock control')
    # Independent graph DP is validated on every atom and the original three seeds.
    small=[];mod=old_graph_module()
    for l in MATS:
        rows=actual_target((l,));h,width,peak=literal_frontier(rows)
        literal=sum(all(rows[u]>>v&1 for u,v in zip(p,p[1:])) for p in permutations(range(len(rows))))
        need(h==literal==vector((l,))[0],('all atom literal frontier',l))
    for seed in SEEDS:
        rows=actual_target(seed);need(literal_frontier(rows)[0]==mod.subset_h(rows)==vector(seed)[0],'old seed literal check')
    families=[];base0=vector(SEEDS[0]+BRIDGE);J=2*(base0[1]+base0[2])-3*base0[0]
    need(base0==(119051,70956,98913) and J==-17415,'common seed and charge')
    for k in range(4):
        records=[]
        for seed in SEEDS:
            word=seed+BRIDGE+('r21',)*k+('b',);rows=actual_target(word)
            h,width,peak=literal_frontier(rows)
            need(h==vector(word)[0],('literal long common family',k,seed))
            records.append({'processing_word':word,'state':vector(word),'H':h,'vertices':len(rows),
                            'physical_port':len(rows)-1,'target_rows':rows,'backward_width':width,'frontier_peak_states':peak})
        m=records[0]['H'];need(records[1]['H']==3*m-10 and records[2]['H']==3*m-20,'joint negative offsets')
        need(records[2]['H']%7==0,'all-k seven-divisible branch check')
        families.append({'k':k,'M':m,'seven_quotient':records[2]['H']//7,'sources':records})
    for k in (0,1):
        need(families[k+2]['M']==40*families[k+1]['M']-13*families[k]['M']-10*J,'M recurrence')
        need(families[k+2]['seven_quotient']==40*families[k+1]['seven_quotient']-13*families[k]['seven_quotient']+74710,'quotient recurrence')
    for key,numerator in (('M',(798527,-928700,304323)),('seven_quotient',(342223,-397900,130387))):
        seq=[f[key] for f in families];d=(1,-41,53,-13)
        need(tuple(sum(d[j]*seq[k-j] for j in range(k+1)) for k in range(3))==numerator,'rational series numerator')
    primes={}
    for n in (798527,21461063791):
        bound=isqrt(n);need(all(n%d for d in range(2,bound+1)),('trial prime',n))
        primes[str(n)]={'trial_division_through':bound,'prime':True}
    need(families[0]['seven_quotient']==7*48889,'first quotient composite hostile')
    q0,q1=(f['seven_quotient'] for f in families[:2]);r=-7471
    U0,U1=13*q0-5*r,13*q1-5*r
    need(U0%13!=0 and (27*U0-U1)%13==0 and ((27*U0-U1)//13)%13!=0,'exact norm13 clock hypotheses')
    clocks=[]
    for modulus in (7,13,169,2197,91,1183):
        c=modular_clock(q0,q1,r,modulus)
        if modulus==7:need(c['tail']==0 and c['period']==6 and c['zero_indices_in_stored_tail_and_cycle']==[0],'prime7 index exclusion')
        if modulus==13:need(c['tail']==1 and c['period']==13 and c['zero_indices_in_stored_tail_and_cycle']==[1],'prime13 index exclusion')
        if modulus in (13,169,2197):
            e={13:1,169:2,2197:3}[modulus]
            need(c['tail']==e and c['period']==modulus,'power13 new source clock')
        if modulus in (91,1183):need((c['tail'],c['period'])==((1,78) if modulus==91 else (2,1014)),'compatible new CRT clocks')
        clocks.append({k:v for k,v in c.items() if k!='states'})
    forbidden=[k for k in range(78) if k%6==0 or k%13==1]
    need(len(forbidden)==18 and 3 not in forbidden,'complete two-prime periodic exclusion')
    need(2 not in forbidden and families[2]['seven_quotient']==29*1327*14057,'passing residue guards does not certify primality')
    return_failures=[]
    for v in ((287,122,225),(289,121,225),(293,119,225),(295,118,225)):
        z=step(v,matrix(BRIDGE));need(decode_state(z) is None,'bridge does not repair second-step state')
        return_failures.append({'before':v,'after':z})
    result={'status':'PASS','bridge_word':BRIDGE,'bridge_matrix':matrix(BRIDGE),'common_initial_state':base0,
            'common_charge':J,'families':families,'prime_certificates':primes,'failed_next_affine_states':return_failures,
            'seven_quotient_clocks':clocks,'forbidden_prime_indices_mod78':forbidden,
            'new_clock_homogeneous_seed':[U0,U1],
            'mixed_state_census':census,'all_return_gains_through_five_letters':[
                {'gain':g,'word':w} for g,w in sorted(gains.items())],
            'loop_generators':loops,'noncommuting_clock_control':noncommuting,
            'scope':'Actual common-source replacement families and exact loop submonoid; no globally closed ternary induction.'}
    (HERE/'channels.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','M':[f['M'] for f in families],
                      'seven_quotients':[f['seven_quotient'] for f in families],
                      'gains':sorted(gains),'literal_max_vertices':max(r['vertices'] for f in families for r in f['sources']),
                      'loop_characteristics':{k:v['characteristic_polynomial'] for k,v in loops.items()},
                      'noncommuting_H':[r['H'] for r in noncommuting]}))
if __name__=='__main__':main()
