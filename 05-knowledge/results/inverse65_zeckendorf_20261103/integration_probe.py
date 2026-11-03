"""Dual error/charge projectors and actual norm-13 suffix replacements."""
from pathlib import Path
from fractions import Fraction
import importlib.util,json,sys
import sympy as sp
from bank import *
from charges import primitive
from inverse import decode_state
HERE=Path(__file__).resolve().parent

def main():
    projectors=[]
    for l in ('b','q1','q2','q4','q5','r21','r38'):
        p=sp.Matrix(MATS[l]);sg=next(s for s in (1,-1) if (p-s*sp.eye(3)).det()==0)
        ell=sp.Matrix(primitive((p-sg*sp.eye(3)).nullspace()[0]))
        tau=sp.Matrix(1,3,primitive((p.T-sg*sp.eye(3)).nullspace()[0]))
        c=int((tau*ell)[0]);pi=ell*tau/c
        need(pi*pi==pi and p*pi==sg*pi and pi*p==sg*pi,'dual spectral projector')
        need(sp.ilcm(*[v.q for v in pi])==abs(c),'minimal integral projector denominator')
        projectors.append({'letter':l,'sign':sg,'charge':list(map(int,ell)),
                           'error_direction':list(map(int,tau)),'pairing':c,'gluing_index':abs(c)})
    # Independent replay of the incoming full-state transport and its failed closure.
    U=('r5','a','a');U2=('a','a','q1','a','a');U4=('q0','a','a','a','a','a')
    t=(2,-1,0);base=vector(U)
    need(base==(31,14,25),'incoming actual seed convention')
    need(vector(U2)==tuple(3*x+y for x,y in zip(base,t)),'incoming +2 source')
    need(vector(U4)==tuple(3*x+2*y for x,y in zip(base,t)),'incoming +4 source')
    need(step(t,MATS['b'])==tuple(-x for x in t) and step(t,MATS['q2'])==t,'incoming signed transport')
    failed=((287,122,225),(289,121,225),(293,119,225),(295,118,225))
    need(all(decode_state(v) is None for v in failed),'incoming four closure failures independently rejected')
    states,parents,ids=closure(2000);tau=(3,-1,-1);hits={2:[],4:[]};both=[]
    for i,v in enumerate(states):
        if v==(1,0,0):continue
        matched=[]
        for k in (2,4):
            target=tuple(3*a+k*b for a,b in zip(v,tau))
            if target in ids:
                matched.append(k);hits[k].append({'state':v,'word':recover(i,parents),
                    'target_state':target,'target_word':recover(ids[target],parents)})
        if len(matched)==2:both.append(v)
    need(len(states)==20703 and [len(hits[k]) for k in (2,4)]==[2,2] and not both,'bounded dual-lift census')
    old=HERE.parent/'golden_hamilton_transfer_20261103'
    spec=importlib.util.spec_from_file_location('independent_graph_dp',old/'transfer.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    selected=[(2,('r4','a','a'),('a','a','r27','a')),
              (4,('a','r12'),('a','a','a','r10'))]
    families=[]
    for m,u,v in selected:
        need(vector(v)==tuple(3*a+m*b for a,b in zip(vector(u),tau)),'actual ternary lift')
        records=[];J=2*sum(vector(u)[1:])-3*vector(u)[0]
        for k in range(4):
            left=u+('r21',)*k+('b',);right=v+('r21',)*k+('b',)
            a,b=vector(left)[0],vector(right)[0]
            need(b==3*a-5*m,'constant-offset transport')
            if k<=1:
                need(mod.subset_h(actual_target(left))==a and mod.subset_h(actual_target(right))==b,'literal offset graphs')
            records.append({'k':k,'base_word':left,'target_word':right,'H':a,'target_H':b,
                            'base_vertices':sum(COSTS[l] for l in left),'target_vertices':sum(COSTS[l] for l in right)})
        need(all(records[k+2]['H']==40*records[k+1]['H']-13*records[k]['H']-10*J for k in range(2)),'new family recurrence')
        families.append({'multiple':m,'charge_J':J,'offset':-5*m,'terms':records})
    result={'status':'PASS','projectors':projectors,'incoming_rejected_states':failed,
            'complete_reachable_states_through2000':len(states),'individual_lifts':hits,'joint_lifts':both,
            'actual_negative_offset_families':families,
            'scope':'Infinite actual suffix replacements; no jointly closed ternary bank.'}
    (HERE/'integration.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','projector_indices':[r['gluing_index'] for r in projectors],
                      'states':len(states),'lift_counts':[len(hits[k]) for k in (2,4)],
                      'initial_offsets':[(f['terms'][0]['H'],f['terms'][0]['target_H']) for f in families]}))
if __name__=='__main__':main()
