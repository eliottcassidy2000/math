"""Least positive defect at1997 and the exact low-branching exclusions."""
from pathlib import Path
import json,sys
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parent.parent
sys.path.insert(0,str(RESULTS/'uniform_return_metagraph_20261111/reset'))
from classify_modules import H,strong,need,inherited as a,two_birth as b
from target_ready_return import family
from finite_bridge import replace
from gap_return_probe import BASE

def power_hit(N,A,B,minimum=1):
    q,r=divmod(N-B,A)
    return r==0 and q>=1<<minimum and q&(q-1)==0

def main():
    old=json.loads((RESULTS/'uniform_return_metagraph_20261111/arithmetic/tail_squares.json').read_text())
    need(len(old['rays'])==52 and len(old['isolated'])==8,'inherited complete low-defect inventory')
    queries=[]
    for N in (1997,5993,6013):
        hits=[r for r in old['rays'] if power_hit(N,r['E'],r['B'])]
        isolates=[r for r in old['isolated'] if r['H']==N]
        three=[(A,B) for A,B in ((7,-5),(9,-3)) if power_hit(N,A,B)]
        triangle=[(A,B) for A,B in ((36,-3),(72,3),(60,-3),(156,9),(168,9)) if power_hit(N,A,B,0)]
        need(not(hits or isolates or three or triangle),'complete branching defect<=6 exclusion')
        queries.append(dict(H=N,transitive_hits=hits,isolates=isolates,three_inport_hits=three,triangle_hits=triangle))
    runs=b.module_runs(BASE,1)
    need(runs[1:3]==[1997,1134] and not any(runs[3:]),'complete constant3-defect chart polynomial')
    controls=[]
    for m in range(1,6):
        T=b.inflate(BASE,1,m);p=a.root_profiles(T)[0]
        need(strong(T) and p['h']==1134*(1<<m)-271 and p['gap']==3,'literal infinite defect3 family')
        need(sum(row&1 for row in T)==1,'indegree-one escape')
        controls.append(dict(m=m,rows=T,profile=p))
    Q=family(1,1);S=a.birth_beaten_only_by(a.birth_beaten_only_by(Q,len(Q)-1),len(Q))
    kept=[0,1,2,3,4,6,7,9,10,11];mapping=[11,2,5,6,7,8,9,10,0,4]
    T,names,rewrite=replace(S,[f'old{i}' for i in range(13)],BASE,kept,mapping,['new1','new3'])
    need(names[0]=='old10' and names[11]=='old0','pre-excursion root and in-neighbour persist')
    need(H(S)==H(T)==1997 and a.root_profiles(T)[0]['gap']==3,'actual same-target least-defect return')
    ready1=json.loads((HERE/'funnel_six_state.json').read_text())['hits']
    ready2=json.loads((HERE/'ready_two_inverse.json').read_text())['hits']
    need(not ready1['1997'] and not ready2['1997'],'complete inherited ready-image exclusions pay minimality')
    result=dict(status='PASS',scope='Complete strong branching strata b<=6; actual indegree-one least-defect3 return at1997',
                low_branching_queries=queries,module_runs=runs,family_controls=controls,actual_return=rewrite,
                ready_5995_counts=[len(ready1['5995']),len(ready2['5995'])],
                ready_6013_counts=[len(ready1['6013']),len(ready2['6013'])])
    (HERE/'low_defect_escape.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS:1997minimumpositive defect3;10/13retainedreturn;no branchingdefect<=6at1997/5993/6013')
    print('REBASE_TRACE',rewrite['trace'])

if __name__=='__main__':main()
