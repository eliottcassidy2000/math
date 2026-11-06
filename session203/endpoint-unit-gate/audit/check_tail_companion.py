#!/usr/bin/env python3
"""Independent tail controls and family-C rooted tensor audit."""
from itertools import permutations, combinations, product
from collections import Counter
import json
from check import need, valid, literal, context, induced, tournaments

def source(bits):
    rows=[0]*(len(bits)+1)
    for i in range(1,len(rows)):
        for j in range(i+1,len(rows)):rows[i]|=1<<j
        if bits[i-1]:rows[0]|=1<<i
        else:rows[i]|=1
    return rows

def rooted_tensor(core,mask,q):
    # Enumerate chosen occupied gaps, rather than multiplying gap polynomials.
    answer=Counter();n=len(core);orders=0
    for word in permutations(range(n)):
        orders+=1
        optional=[];mandatory=[];good=True
        for g in range(n+1):
            left=word[g-1] if g else None
            right=word[g] if g<n else None
            empty=left is None or right is None or bool(core[left]&(1<<right))
            filled=(left is None or left==q or not mask&(1<<left)) and (right is None or right==q or bool(mask&(1<<right)))
            if not empty and not filled:good=False;break
            if filled:
                item=(int(right==q),int(left==q))
                if empty:optional.append(item)
                else:mandatory.append(item)
        if not good:continue
        for choices in range(1<<len(optional)):
            occupied=mandatory+[item for i,item in enumerate(optional) if choices&(1<<i)]
            key=(len(occupied),sum(v[0] for v in occupied),sum(v[1] for v in occupied))
            need(key[1]<=1 and key[2]<=1,'double root adjacency')
            answer[key]+=1
    return answer,orders

def main():
    tails=[];calls=0
    for t in range(1,8):
        prefix=source([1]+[0]*t)
        ph,ps,pc,w=literal(prefix);calls+=w
        need(ph==2**t+1,'prefix formula')
        hs=[]
        for a,b in product((0,1),repeat=2):
            U=source([1]+[0]*t+[a,b,0]);valid(U)
            h,s,c,w=literal(U);calls+=w
            need(s==1 and c>0 and h>=ph,'tail physical source')
            hs.append(h)
        need(hs[3]-hs[2]-hs[1]+hs[0]==-2,'tail twin curvature')
        tails.append({'t':t,'prefix_H':ph,'four_H':hs})
    P=[2,4,41,3,15,91,31];C=[0,5,41,131,15,219,159,23]
    pt,po=rooted_tensor(P,27,5);ct,co=rooted_tensor(C,155,5)
    delta={key:ct[key]-2*pt[key] for key in set(pt)|set(ct) if ct[key]!=2*pt[key]}
    need(delta=={(0,0,0):1,(1,0,1):1,(1,1,0):-1,(2,1,1):-1},'family-C tensor')
    # Relabel copied root5 as0 to use the independently implemented physical context.
    pp=[5]+[i for i in range(len(P)) if i!=5]
    cc=[5]+[i for i in range(len(C)) if i!=5]
    kp=induced(P,pp);kc=induced(C,cc)
    mp=sum(1<<i for i,v in enumerate(pp) if 27&(1<<v))
    mc=sum(1<<i for i,v in enumerate(cc) if 155&(1<<v))
    controls=0;singleton=None;triangle=[]
    for n in range(1,5):
        for U in tournaments(n):
            h,s,_,w=literal(U);calls+=w
            p=context(kp,mp,U);c=context(kc,mc,U)
            hp,_,_,w=literal(p);calls+=w
            hc,_,_,w=literal(c);calls+=w
            need(hc-2*hp==2*s-h,'family-C literal identity')
            # Original core vertex0 became new label1 and remains a sink.
            need(c[1]==0,'family-C universal sink')
            red=induced(c,[v for v in range(len(c)) if v!=1])
            hr,_,_,w=literal(red);calls+=w
            need(hr==hc and len(red)==len(p),'family-C replacement')
            controls+=1
            if n==1:singleton={'parent_H':hp,'child_H':hc,'correction':hc-2*hp}
            if n==3 and h==3:triangle.append({'parent_H':hp,'child_H':hc,'correction':hc-2*hp})
    need(controls==75,'rooted family-C input census')
    print(json.dumps({'status':'PASS','tail_controls':tails,
                      'parent_tensor':{str(k):v for k,v in sorted(pt.items())},
                      'child_tensor':{str(k):v for k,v in sorted(ct.items())},
                      'difference_tensor':{str(k):v for k,v in sorted(delta.items())},
                      'core_orders_checked':po+co,'family_C_root_zero_inputs':controls,
                      'singleton_control':singleton,'triangle_controls':triangle,
                      'literal_DFS_calls':calls},indent=2,sort_keys=True))

if __name__=='__main__':main()
