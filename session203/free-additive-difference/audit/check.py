#!/usr/bin/env python3
"""Independent literal path enumeration for freely prescribed positive gaps."""
from itertools import combinations
import json

P_CORE=[2,0,11,19,7]
C_CORE=[2,8,3,5]

def need(b,s):
    if not b:raise RuntimeError(s)

def valid(a):
    n=len(a)
    for i,row in enumerate(a):
        need(0<=row<1<<n and not row&(1<<i),'row')
        for j in range(i):need(((row>>j)&1)+((a[j]>>i)&1)==1,'pair')

def tournaments(n):
    pairs=list(combinations(range(n),2))
    for code in range(1<<len(pairs)):
        a=[0]*n
        for bit,(i,j) in enumerate(pairs):
            if code&(1<<bit):a[i]|=1<<j
            else:a[j]|=1<<i
        yield a

def glue(core,U):
    n=len(core);a=core[:]+[0]*(len(U)-1)
    image=[0]+list(range(n,len(a)))
    for i,j in combinations(range(len(U)),2):
        x,y=image[i],image[j]
        if U[i]&(1<<j):a[x]|=1<<y
        else:a[y]|=1<<x
    for core_v in range(1,n):
        for old_v in range(n,len(a)):a[core_v]|=1<<old_v
    valid(a);return a

def enumerate_paths(a,profiles=False):
    n=len(a);full=(1<<n)-1;H=S=C=calls=0
    h=[0]*(1<<n);end0=[0]*(1<<n);start0=[0]*(1<<n)
    h[0]=1
    def visit(last,mask,root):
        nonlocal H,S,C,calls
        calls+=1
        if profiles:
            h[mask]+=1
            if last==0:end0[mask]+=1
            if root==0:start0[mask]+=1
        if mask==full:
            H+=1
            if root==0:
                S+=1;C+=bool(a[last]&1)
            return
        todo=a[last]&~mask
        while todo:
            bit=todo&-todo;todo-=bit
            visit(bit.bit_length()-1,mask|bit,root)
    for v in range(n):visit(v,1<<v,v)
    result={'H':H,'start':S,'cycles':C,'calls':calls}
    if profiles:
        L=sum(end0[mask]*h[full^mask] for mask in range(1,full+1))
        R=sum(start0[mask]*h[full^mask] for mask in range(1,full+1))
        d=h[full^1]
        beta=sum(h[mask]*h[(full^1)^mask] for mask in range(full+1) if not mask&1)
        result.update({'L':L,'R':R,'end':end0[full],'deletion_H':d,'deletion_B2':beta})
    return result

def cycle_chain(m):
    # p0,...,p(m-1),x,y, with chain -> x -> y -> chain.
    a=[0]*(m+2);x=m;y=m+1
    for i in range(m):
        a[i]=sum(1<<j for j in range(i+1,m))|(1<<x)
    a[x]=1<<y;a[y]=(1<<m)-1
    valid(a);return a

def endpoint_source(s):
    m=s.bit_length();W=cycle_chain(m)
    a=[0]+[row<<1 for row in W]
    for i in range(m+2):
        beaten=i<m and bool(s&(1<<(m-1-i)))
        if beaten:a[0]|=1<<(i+1)
        else:a[i+1]|=1
    valid(a);return a

def starts_all(a):
    n=len(a);full=(1<<n)-1;out=[0]*n
    def visit(last,mask,root):
        if mask==full:out[root]+=1;return
        todo=a[last]&~mask
        while todo:
            bit=todo&-todo;todo-=bit
            visit(bit.bit_length()-1,mask|bit,root)
    for root in range(n):visit(root,1<<root,root)
    return out

def main():
    core_profiles=[];calls=0;graphs=0
    for core in (P_CORE,C_CORE):
        valid(core);z=enumerate_paths(core,True);calls+=z.pop('calls');graphs+=1
        core_profiles.append(z)
    need([(p['deletion_H'],p['R'],p['H']-p['end']) for p in core_profiles]==[(3,6,3)]*2,'core equalities')
    need([p['end'] for p in core_profiles]==[0,2],'core endpoint difference')
    need([p['deletion_B2'] for p in core_profiles]==[24,12],'higher profile mismatch')
    ordinary=0
    for n in range(1,5):
        for U in tournaments(n):
            z=enumerate_paths(U,True);calls+=z['calls'];graphs+=1
            h,s,L=z['H'],z['start'],z['L']
            outputs=[]
            for core in ([2,0],P_CORE,C_CORE):
                a=glue(core,U);w=enumerate_paths(a);calls+=w['calls'];graphs+=1;outputs.append(w)
            edge,p,c=outputs
            need(edge['H']==h+L-s and p['H']==3*edge['H'] and c['H']==p['H']+2*s,'all-input identity')
            need(len(glue(P_CORE,U))==n+4 and len(glue(C_CORE,U))==n+3,'all-input orders')
            if z['cycles']>0 or n==1:need(c['cycles']>0,'strong child')
            ordinary+=1
    chain_controls=[]
    for m in range(1,10):
        W=cycle_chain(m);starts=starts_all(W)
        need(starts==[2**(m-i-1) for i in range(m)]+[1,1],'binary start weights')
        need(sum(starts)==2**m+1,'W path count')
        chain_controls.append({'m':m,'start_counts':starts,'H':sum(starts)})
    freely=0;largest=0;samples=[]
    for s in range(1,128):
        U=endpoint_source(s);z=enumerate_paths(U,True);calls+=z['calls'];graphs+=1
        need(z['start']==s and z['cycles']>0,'free strong endpoint source')
        need(len(U)==s.bit_length()+3,'source order')
        p=glue(P_CORE,U);c=glue(C_CORE,U)
        zp=enumerate_paths(p);zc=enumerate_paths(c);calls+=zp['calls']+zc['calls'];graphs+=2
        need(zp['H']==3*(z['H']+z['L']-s) and zc['H']-zp['H']==2*s,'free additive identity')
        need(zp['H']%2==zc['H']%2==1 and zc['cycles']>0,'odd counts and strong child')
        need(len(p)==s.bit_length()+7 and len(c)==s.bit_length()+6,'output orders')
        need(zp['start']==0 and zc['start']==2**s.bit_length()+1,'retained endpoint coupling control')
        largest=max(largest,zc['H']);freely+=1
        if s in (1,2,3,4,5,7,8,15,16,31,32,63,64,85,127):
            samples.append({'s':s,'source_H':z['H'],'source_L':z['L'],'parent_H':zp['H'],'child_H':zc['H'],
                            'parent_order':len(p),'child_order':len(c),'parent_root_start':zp['start'],'child_root_start':zc['start']})
    need(ordinary==75 and freely==127 and largest==21773,'author control totals')
    print(json.dumps({'status':'PASS','core_profiles':core_profiles,'ordinary_inputs':ordinary,
                      'free_positive_endpoint_targets':freely,'largest_child_H':largest,
                      'literal_graphs':graphs,'literal_DFS_visits':calls,'chain_weight_controls':chain_controls,
                      'free_source_samples':samples,
                      'scope':'Positive result in sections1-2 only. Exploratory section3 is outside this audit. No freely selectable base count or closed full-profile iteration.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
