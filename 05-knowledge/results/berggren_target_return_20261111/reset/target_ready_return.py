"""All-k target-preserving ready returns after the actual (+4,+2) path."""
from pathlib import Path
import json,sys
HERE=Path(__file__).resolve().parent
RESULTS=HERE.parent.parent
sys.path.insert(0,str(RESULTS/'uniform_return_metagraph_20261111'/'reset'))
from classify_modules import H,strong,need,inherited as a,two_birth as b,converse
from atlas_two_birth import construct
from ready_rebase import full_state

P221=(34,244,249,243,225,192,129,1)
P665=(162,1020,1017,1009,993,960,897,768,513,1)

def parent(k,stage=0):return b.inflate(P221 if stage==0 else P665,5 if stage==0 else 7,k)
def family(k,stage=0):return a.birth_beaten_only_by(parent(k,stage),0)

def expand_map(kept,mapping,source_module,target_module,k):
    out_keep=[];out_map=[]
    for u,v in zip(kept,mapping):
        if u==source_module:
            need(v==target_module,'matched whole transitive module')
            out_keep.extend(range(u,u+k));out_map.extend(range(v,v+k))
        else:
            out_keep.append(u+(k-1 if u>source_module else 0))
            out_map.append(v+(k-1 if v>target_module else 0))
    return out_keep,out_map

def mark_runs(P,module):
    n=len(P);root=0;R=a.birth_beaten_only_by(P,root)
    separator=tuple(row|(1<<n) for row in P)+((1<<n)-1,)
    force_start=tuple(row&~(1<<root) for row in separator)
    force_end=tuple(0 if i==root else row for i,row in enumerate(P))
    start_paths=tuple(row&~(1<<root) for row in P)
    return {name:b.module_runs(T,module)[1:] for name,T in
            [('parent',P),('child',R),('two_covers',separator),
             ('start_covers',force_start),('end_paths',force_end),('start_paths',start_paths)]}

def falling_runs(k,count):
    S=[[0]*(count+1) for _ in range(k+1)];S[0][0]=1
    for n in range(1,k+1):
        for j in range(1,min(n,count)+1):S[n][j]=j*S[n-1][j]+S[n-1][j-1]
    result=[];fact=1
    for j in range(1,count+1):fact*=j;result.append(fact*S[k][j])
    return result

def evaluate(co,k):return sum(x*y for x,y in zip(co,falling_runs(k,len(co))))

def state_from_runs(runs,k):
    h=evaluate(runs['child'],k);q=evaluate(runs['parent'],k)
    K=evaluate(runs['two_covers'],k);U0=evaluate(runs['start_covers'],k)
    return [h,evaluate(runs['end_paths'],k),q-1,K-U0,h+1,q,K]

def main():
    runs=[mark_runs(P221,5),mark_runs(P665,7)]
    for j,(h,A) in enumerate(((221,99),(665,297))):
        need(runs[j]['child'][:2]==[h,A] and not any(runs[j]['child'][2:]),'complete target polynomial')
        need(runs[j]['start_paths'][0]==1 and not any(runs[j]['start_paths'][1:]),'all-k unique-start parent')
    controls=[]
    for k in range(1,7):
        D=1<<k;S=construct((1,1,2),0,4,14,k)[-1]
        R=family(k,0);Q=family(k,1)
        need(H(S)==H(R)==99*D+23 and H(Q)==297*D+71,'literal independent target counts')
        for stage,T in enumerate((R,Q)):
            need(strong(parent(k,stage)) and strong(T),'strong ready parents and children')
            need(full_state(T,len(T)-1)==state_from_runs(runs[stage],k),'full seven-coordinate independent count')
        # Direct same-gauge return. Old root and the two newborn tail vertices survive.
        kept,mapping=expand_map([0,5,6,7],[5,1,6,7],0,5,k)
        need(a.induced(S,kept)==a.induced(R,mapping),'actual retained source to ready target, same gauge')
        # Better-retained return in the converse gauge.
        kept_c,map_c=expand_map([0,2,4,5,6,7],[5,1,0,6,4,2],0,5,k)
        map_c[:k]=reversed(map_c[:k])
        need(a.induced(S,kept_c)==a.induced(converse(R),map_c),'actual six-of-eight converse return')
        first=a.birth_beaten_only_by(R,len(R)-1)
        next_target=a.birth_beaten_only_by(first,len(first)-1)
        kept2,map2=expand_map([0,2,3,4,5,6,7,8,9],[0,3,4,6,7,8,9,1,2],5,7,k)
        need(a.induced(next_target,kept2)==a.induced(Q,map2),'second target-preserving return retains all but two')
        need(H(next_target)==H(Q)==3*H(R)+2,'exact absolute target preserved by second return')
        following=a.two_birth_reset(Q,len(Q)-1)
        need(H(following)==891*D+215,'next exact target before readiness obstruction')
        controls.append(dict(k=k,source=S,ready221=R,ready665=Q,
            first_kept=kept,first_map=mapping,converse_kept=kept_c,converse_map=map_c,
            second_kept=kept2,second_map=map2,
            trace=[11*D+1,33*D+7,H(S),H(R),H(next_target),H(Q),H(following)],
            states=[state_from_runs(r,k) for r in runs]))
    result=dict(status='PASS',scope='Complete fixed-module run identities; literal independent controls k1..6; retained block maps',
                parent_templates=[P221,P665],module_vertices=[5,7],marked_runs=runs,controls=controls)
    (HERE/'target_ready_return.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: target-preserving all-k returns99D+23 and297D+71; exact retained maps and full states')
    print('MARKED_RUNS',runs)

if __name__=='__main__':main()
