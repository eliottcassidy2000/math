"""Independent forward endpoint audit of the sibling two-birth formulas."""
from pathlib import Path
from itertools import combinations
import importlib.util,json,hashlib
from check_coverage import profile,birth,need,strong,brute_h
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'construction/algebra.py'
spec=importlib.util.spec_from_file_location('candidate_algebra',SOURCE);alg=importlib.util.module_from_spec(spec);spec.loader.exec_module(alg)
def c2(T,a):
    n=len(T);return birth(birth(T,((1<<n)-1)^(1<<a)),(1<<n)-1)
def main():
    inputs=[]
    for n in range(1,4):
        edges=list(combinations(range(n),2))
        for mask in range(1<<len(edges)):
            T=[0]*n
            for k,(a,b) in enumerate(edges):
                if mask>>k&1:T[a]|=1<<b
                else:T[b]|=1<<a
            inputs.append(tuple(T))
    inputs.extend([(10,4,1,6),(14,12,8,0),(20,1,2,7,14)])
    count=0;strong_count=0
    for T in inputs:
        n=len(T)
        for w in range(n):
            for a in range(n):
                z=alg.operator(T,w,a);U=c2(T,a);q=profile(U);V=birth(birth(T,((1<<n)-1)^(1<<a)),1<<w);v=profile(V)
                need(tuple(z['C2'])==U and tuple(z['C4'])==V,'independent literal construction')
                need(tuple(z['predicted2'])==(q['H'],q['E2'][w]),'C2 exact polynomial output')
                need(tuple(z['predicted4'])==(v['H'],v['E2'][w]),'C4 exact polynomial output')
                need(q['S2'][n]==q['H'],'transported auxiliary x is start-ready')
                need(q['E1'][n]==z['h'],'new x end count remembers ancestor H')
                need(q['S2'][n+1]-q['H']==z['h']-z['e_a'],'new y start defect is ancestor h minus end count')
                z2=alg.operator(U,w,n);need(z2['L']+z2['K']==q['E2'][w],'full-root coordinate transport invariant')
                W=c2(U,n);r=profile(W);need((r['H'],r['E2'][w])==(3*q['H'],3*q['E2'][w]),'second macro scales actual h/F by3')
                Z=c2(U,n+1);zz=profile(Z)
                need(zz['H']==3*q['H']+2*(z['h']-z['e_a']) and zz['E2'][w]==3*q['E2'][w]+2*(z['F']-z['J']),'alternative newborn carries both ancestor debts')
                need(z['F']-z['J']>=profile(T)['E1'][w],'endpoint empty-slot positive debt')
                keep=[i for i in range(n) if i!=a];D=tuple(sum(((T[i]>>j)&1)<<k for k,j in enumerate(keep)) for i in keep);dq=profile(D)
                R=0 if a==w else dq['E2'][keep.index(w)];old_a=profile(c2(U,a))
                need(q['S2'][a]==4*z['S_a_2']-dq['H'] and old_a['E2'][w]==q['E2'][w]+2*(4*(z['L']+z['K'])-R),'retained old-a response from actual second macro')
                need(z['F']-R>=profile(T)['E1'][w],'deletion-insertion positive old-a debt')
                if strong(T) and len(T)>=4:need(all(f-q['H'] not in (1,2) for f in q['S2']),'scoped inherited all-marking obstruction')
                if strong(T):need(strong(U) and strong(W),'strong source retained after both macros');strong_count+=1
                count+=1
    T=(2,4,1);U=c2(T,1);q=profile(U);good=c2(U,3);bad=c2(U,4)
    need(brute_h(U)==q['H'] and brute_h(good)==profile(good)['H'] and brute_h(bad)==profile(bad)['H'],'independent permutation selected clocks')
    need(profile(good)['H']==3*q['H'] and profile(bad)['H']!=3*q['H'],'wrong-newborn auxiliary is a hostile control')
    # Every C2 image has a start-ready vertex. The strong11/13 exclusion is inherited;
    # this bounded independent census checks its small-order core boundary.
    start_ready_counts=set();strongs=0
    for n in range(1,6):
        edges=list(combinations(range(n),2))
        for mask in range(1<<len(edges)):
            T=[0]*n
            for k,(a,b) in enumerate(edges):
                if mask>>k&1:T[a]|=1<<b
                else:T[b]|=1<<a
            T=tuple(T)
            if not strong(T):continue
            strongs+=1;q=profile(T)
            if any(f==q['H'] for f in q['S2']):start_ready_counts.add(q['H'])
    need(not ({11,13}&start_ready_counts),'small-core start-ready exclusion')
    old=HERE.parent.parent/'zero_ready_strong_switch_20261103/literal/realized_generation2_seed4_depth3_all_prefixes.json'
    old=json.loads(old.read_text())[1];q=profile(tuple(old['suffix_rows']));names=old['suffix_names'];start=[f-q['H'] for f in q['S2']];end=[f-q['H'] for f in q['E2']]
    need(q['H']==80319 and min(start)==616 and min(end)==2960 and all(v>0 for v in start+end),'exact escaped suffix has neither polarity zero-ready')
    barrier=dict(outer_H=1537279,carrier='induced named B14',rows=old['suffix_rows'],names=names,H=q['H'],start_defects=dict(zip(map(str,names),start)),end_defects=dict(zip(map(str,names),end)),outdegree_one_names=[n for n,row in zip(names,old['suffix_rows']) if row.bit_count()==1],minimum_start_defect=616,minimum_C2_H=3*q['H']+2*min(start),scope='No auxiliary re-marking gives a zero-offset C2 scalar triple on this exact B14; this does not rule out another actual entry or macro.')
    result=dict(status='FINITE_EXACT_INDEPENDENT_C2_C4_AND_AUXILIARY_TRANSPORT_AUDIT',universe='all labelled inputs orders1..3 and every w,a, plus displayedTT4/Q4/strongH9 controls',formula_cases=count,strong_cases=strong_count,source_sha256=hashlib.sha256(SOURCE.read_bytes().replace(b'\r\n',b'\n')).hexdigest(),auxiliary_hostile=dict(source_rows=U,source_H=profile(U)['H'],correct_auxiliary=3,correct_H=profile(good)['H'],wrong_auxiliary=4,wrong_H=profile(bad)['H']),strong_small_core_control=dict(universe='all labelled strong tournaments orders1..5',sources=strongs,start_ready_counts=sorted(start_ready_counts)),escaped_suffix=barrier,scope='Independent finite formula/transport controls. General proof is the sibling construction theorem. Strong11/13 all-order exclusion also requires inherited strong-order lower bound, not this finite census alone. The no1/2-start-ready C2-output obstruction requires strong INPUT order>=4; C3 input is an explicit exception.')
    (HERE/'c2_audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: '+str(count)+' exact C2/C4 cases, transported-auxiliary strong tripling, wrong-auxiliary hostile, and all strong small-core11/13 controls')
    print('auxiliary hostile',result['auxiliary_hostile'])
if __name__=='__main__':main()
