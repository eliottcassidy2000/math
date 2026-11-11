"""Actual alternating23->73->221->667->2003 path with persistent root/port."""
from pathlib import Path
import json,sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent.parent/'uniform_return_metagraph_20261111'/'reset'))
from classify_modules import H,strong,need,inherited as a,two_birth as b
from atlas_two_birth import construct

B221=(34,108,105,113,103,64,129,63)
B667=(6,196,200,211,231,207,129,257,127)
B2003=(514,284,505,497,481,451,387,259,1,510)

def replace(S,names,T,kept,mapping,new_names):
    need(a.induced(S,kept)==a.induced(T,mapping),'every retained actual arc')
    labels=list(mapping);current=a.induced(S,kept);trace=[H(S),H(current)]
    target_names={v:names[u] for u,v in zip(kept,mapping)}
    missing=[v for v in range(len(T)) if v not in labels]
    need(len(missing)==len(new_names),'exact replacement names')
    for v,name in zip(missing,new_names):
        out=sum((T[v]>>w&1)<<j for j,w in enumerate(labels))
        current=a.attach(current,out);labels.append(v);target_names[v]=name;trace.append(H(current))
    need(current==a.induced(T,labels),'literal sequence of declared arbitrary births')
    result_names=[target_names[v] for v in range(len(T))]
    need(len(set(result_names))==len(result_names),'physical names stay distinct')
    return T,result_names,dict(kept=kept,map=mapping,source_names=names,target_names=result_names,
                               insertion_target_labels=missing,trace=trace)

def main():
    T23=a.ready_source((0,1),1);T73=a.ready_source((1,1,2),1)
    S,ids,P,X,V=construct((0,1),0,4,8,1)
    need(S==T23 and a.induced(T73,[0,1,4,5,6,2,3])==V,'inherited actual first internal+4 map')
    S,ids2,P2,X2,T221=construct((1,1,2),0,4,14,1)
    need(S==T73,'same actual intermediate source for second digit')
    # Original names are P,old00,old10,a,b,r. Newborns on the two inherited edges are x,y,u,v.
    names=['P','x','y','a','b','r','u','v']
    ready,names,first=replace(T221,names,B221,[0,1,3,5,6,7],[4,2,6,7,0,1],['f','g'])
    need(names[7]=='r' and names[6]=='a','first return retains original root and original port')
    need(a.root_profiles(ready)[7]['gap']==2,'actual2-ready221 preparation')
    born1=a.birth_beaten_only_by(ready,7);born2=a.birth_beaten_only_by(born1,8)
    need(H(born2)==667,'standard+4 at prepared221')
    names+=['c','d']
    ready667,names,second=replace(born2,names,B667,[2,3,4,5,6,7,8],[1,2,3,6,7,8,5],['e','h'])
    need(names[8]=='r' and names[7]=='a','667 return retains original root and port')
    need(a.root_profiles(ready667)[8]['gap']==2,'actual2-ready667 preparation')
    last,names,third=replace(ready667,names,B2003,[1,2,4,5,6,7,8],[4,5,2,3,6,0,9],['j','l','m'])
    need(names[9]=='r' and names[0]=='a','direct667to2003 retains original root and port')
    need(a.root_profiles(last)[9]['gap']==1 and H(last)==3*667+2,'direct+2 withready output despite667not1-ready')
    for T in(T23,T73,T221,B221,born1,born2,B667,B2003):need(strong(T),'all named completed sources are strong')
    sr=[b.module_runs(B667,v) for v in range(len(B667))]
    tr=[b.module_runs(B2003,v) for v in range(len(B2003))]
    lift_matches=[]
    for v,co in enumerate(sr):
        for w,ct in enumerate(tr):
            if all((ct[j] if j<len(ct) else 0)==3*(co[j] if j<len(co) else 0)
                   for j in range(2,max(len(co),len(ct)))):lift_matches.append((v,w))
    need(not lift_matches,'complete all-single-module uniform lift obstruction for this bridge')
    result=dict(status='PASS',scope='Literal finite alternating path; all retained comparisons and actual births; all90 single-module pairs',
        target_counts=[23,73,221,667,2003],original_persistent_names=['a','r'],
        initial_plus4_trace=[H(T) for T in(T23,P,X,V)],
        initial_plus2_trace=[H(T) for T in(T73,P2,X2,T221)],
        ready221_rows=B221,ready667_rows=B667,ready2003_rows=B2003,
        rewrites=[first,second,third],source_module_runs=sr,target_module_runs=tr,lift_matches=lift_matches)
    (HERE/'finite_bridge.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: actual23->73->221->667->2003; same originalrootr andporta survive all returns')
    print('REWRITE_TRACES',[x['trace'] for x in result['rewrites']])

if __name__=='__main__':main()
