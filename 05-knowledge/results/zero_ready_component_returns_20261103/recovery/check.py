#!/usr/bin/env python3
"""Recovered marked-component routes and exact hostile probes, with bounded controls."""
from pathlib import Path
from itertools import permutations
import hashlib,json
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]

def need(v,msg):
    if not v:raise RuntimeError(msg)

def tables(rows):
    n=len(rows);full=(1<<n)-1;end=[[0]*n for _ in range(full+1)];start=[[0]*n for _ in range(full+1)]
    H=[0]*(full+1);H[0]=1
    incoming=[sum(1<<u for u in range(n) if rows[u]>>v&1) for v in range(n)]
    for mask in range(1,full+1):
        for v in range(n):
            if not mask>>v&1:continue
            old=mask^(1<<v)
            if not old:end[mask][v]=start[mask][v]=1;continue
            bits=old&incoming[v];e=0
            while bits:
                b=bits&-bits;bits-=b;e+=end[old][b.bit_length()-1]
            end[mask][v]=e;bits=old&rows[v];s=0
            while bits:
                b=bits&-bits;bits-=b;s+=start[old][b.bit_length()-1]
            start[mask][v]=s
        H[mask]=sum(end[mask])
    return H,end,start

def fields(rows,r):
    H,end,start=tables(rows);full=(1<<len(rows))-1
    E2=sum(end[S][r]*H[full^S] for S in range(1,full+1) if S>>r&1)
    S2=sum(start[S][r]*H[full^S] for S in range(1,full+1) if S>>r&1)
    P2=sum(H[S]*H[full^S] for S in range(full+1))
    return {'h':H[full],'s':start[full][r],'q':H[full^(1<<r)],'F':E2,'e':end[full][r],'S2':S2,'P2':P2}

def born(rows,r):
    n=len(rows);out=list(rows)+[0]
    for v in range(n):
        if v==r:out[v]|=1<<n
        else:out[n]|=1<<v
    return out

def transfer(z,cmd):
    h,s,q,F=z
    return (h+F-s,q,q,2*F-q) if cmd=='R' else (h+F-s,h-s,h,h+F)

def seam(X,a,Y,b):
    n=len(X);m=len(Y);c=n+m;out=[0]*(c+1)
    for u,row in enumerate(X):out[u]=row
    for u in range(m):
        out[n+u]=(1<<n)-1
        for v in range(m):
            if Y[v]>>u&1:out[n+u]|=1<<(n+v)
    for u in range(n):
        if u==a:out[u]|=1<<c
        else:out[c]|=1<<u
    for u in range(m):
        if u==b:out[c]|=1<<(n+u)
        else:out[n+u]|=1<<c
    return out,c

def literal_h(rows):
    return sum(all(rows[p[i]]>>p[i+1]&1 for i in range(len(rows)-1)) for p in permutations(range(len(rows))))

def small_seam_controls():
    inputs=[([0],0),([2,0],0),([2,0],1),([6,4,0],0),([6,4,0],1),([6,4,0],2),([2,4,1],0),([2,4,1],1),([2,4,1],2)]
    ncheck=0
    for X,a in inputs:
        x=fields(X,a)
        for Y,b in inputs:
            y=fields(Y,b);Z,c=seam(X,a,Y,b);z=fields(Z,c)
            eh=x['F']*y['F']+x['h']*y['F']+x['F']*y['h']+(x['h']-x['s'])*(y['h']-y['s'])
            eF=x['F']*y['P2']+x['h']*(y['P2']-y['S2'])
            need((z['h'],z['F'])==(eh,eF),'seam count/formula')
            need(literal_h(Z)==z['h'],'seam literal Hamilton check')
            for s in range(1,4):
                t=1<<s
                A=y['F']+y['h']-y['s']+t*(y['P2']-y['S2'])
                B=y['F']+y['h']+t*y['P2']
                C=y['h']-y['s']
                need(eh+t*eF==A*x['h']+B*x['F']-C*x['s'],'marked seam linear functional')
                ncheck+=1
    return {'rooted_pairs':81,'outer_functional_checks':ncheck,'literal_hamilton_checks':81}

def deletion(rows,w):
    V=[i for i in range(len(rows)) if i!=w]
    return [sum((rows[u]>>v&1)<<j for j,v in enumerate(V)) for u in V]

def p2(rows):
    H,_,_=tables(rows);full=(1<<len(rows))-1
    return sum(H[S]*H[full^S] for S in range(full+1))

def six_fields(rows,w):
    f=fields(rows,w)
    return [f['h'],f['e'],f['S2'],f['F'],p2(deletion(rows,w)),f['q']]

def append_root(rows,w,bit):
    out=list(rows)+[0];n=len(rows)
    for v in range(n):
        if bit==0 and v==w:out[n]|=1<<v
        else:out[v]|=1<<n
    return out

def six_transfer(z,bit):
    h,e,t,F,p,q=z
    return [h,0,2*t,F,2*p,q] if bit else [h+t-e,q,2*t-q,F+p,2*p,q]

def six_controls():
    count=0
    for n in range(1,5):
        for mask in range(1<<(n*(n-1)//2)):
            out=[0]*n;k=0
            for a in range(n):
                for b in range(a+1,n):
                    if mask>>k&1:out[a]|=1<<b
                    else:out[b]|=1<<a
                    k+=1
            for w in range(n):
                z=six_fields(out,w)
                for bit in (0,1):
                    need(six_fields(append_root(out,w,bit),w)==six_transfer(z,bit),'retained-mark six-scalar birth law')
                    count+=1
    return count

def joint_lower_bound_hostile():
    rows=[8,1,3,6];w=1;H,E,S=tables(rows);full=15
    cycle=[0,3,2,1,0]
    need(all(rows[a]>>b&1 for a,b in zip(cycle,cycle[1:])), 'strong hostile Hamilton cycle')
    joint=[]
    for v in range(4):
        def rec(mask,last):
            if last==w:return H[full^mask]
            return sum(rec(mask|1<<u,u) for u in range(4) if rows[last]>>u&1 and not mask>>u&1)
        joint.append(rec(1<<v,v))
    need(S[full][w]==1 and H[full^(1<<w)]==3 and joint==[2,3,2,2], 'strong unique-start joint lower-bound hostile')
    # The only nontrivial strong tournaments below order4 are the two labelled C3 orientations.
    small=[]
    for tri in ([2,4,1],[4,1,2]):
        for root in range(3):
            z=three_birth_dial(tri,root,list(range(3)))
            need(min(x[1] for x in z['old_joint_same_bin_starts_ends_w'])>=1,'C3 minimality control')
            small.append([tri,root])
    return dict(status='REFUTED_PROPOSED_SIMPLIFICATION_ONLY',proposal='For strong unique-start(D,w), every J_(v,w)>=H(D-w).',rows=rows,root=w,strong_cycle=cycle,h=H[full],q=3,joint=joint,hostile_vertex=0,minimal_nontrivial_strong_order=4,smaller_controls=len(small))

def three_birth_dial(rows,w,names):
    H0,e0,s0=tables(rows);n=len(rows);full=(1<<n)-1
    f=fields(rows,w);h,s,F=f['h'],f['s'],f['F']
    joint=[]
    for v in range(n):
        def rec(mask,last):
            if last==w:return H0[full^mask]
            out=0;more=rows[last]&~mask
            while more:
                bit=more&-more;more-=bit;out+=rec(mask|bit,bit.bit_length()-1)
            return out
        joint.append(rec(1<<v,v))
    X=born(rows,w);X=born(X,n);X=born(X,n+1)
    HX,eX,sX=tables(X);fx=(1<<len(X))-1
    profiles=[];observations=[]
    for root in (n,n+1,n+2):
        Froot=sum(eX[S][root]*HX[fx^S] for S in range(1,fx+1) if S>>root&1)
        profiles.append([HX[-1],sX[-1][root],HX[fx^(1<<root)],Froot])
        S2=sum(sX[S][root]*HX[fx^S] for S in range(1,fx+1) if S>>root&1)
        E1=eX[-1][root];need((S2-E1)%2==0,'golden parity')
        gold=[HX[-1],(S2-E1)//2,(S2+E1)//2]
        hh,xx,yy=gold
        observations.append({'name':['dial_a','dial_b','dial_c'][root-n],
          'E1':E1,'S2':S2,'golden':gold,'F49_charge':[(4*(hh+2*xx-yy))%7,yy%7,(hh-yy)%7]})
    hX=3*h+3*F-s
    expected=[[hX,h,h,3*h+4*F],[hX,h+F-s,h+F-s,hX],[hX,h+F,h+2*F,3*h+4*F-s]]
    need(profiles==expected,'three-root dial profiles')
    need(sX[-1][:n]==joint,'all old starts equal marked same-bin cover counts')
    need(HX[fx^(1<<w)]==3*f['q'],'old-mark deletion tripling')
    wrapper_pairs=[]
    for clock in (1,2):
        t=1<<clock
        c=hX+t*profiles[2][3];a=hX+t*profiles[0][3]
        need(a-c==t*s,'contextual dial jump')
        wrapper_pairs.append({'clock':clock,'newest_root_value':c,'first_root_value':a,'difference':a-c})
    return {'rows':X,'names':names+['dial_a','dial_b','dial_c'],
      'new_root_profiles':{name:z for name,z in zip(['dial_a','dial_b','dial_c'],profiles)},
      'root_observations':observations,
      'old_joint_same_bin_starts_ends_w':list(zip(names,joint)),
      'all_start_counts':sX[-1],'minimum_start_count':min(sX[-1]),
      'old_w_deletion':HX[fx^(1<<w)],'wrapper_pairs':wrapper_pairs}

def main():
    path='05-knowledge/results/triangle_return_170807_20261103/zero_ready.json'
    src=json.loads((ROOT/path).read_text())[1];rows=src['suffix_rows'];names=src['suffix_names'];r=names.index(31)
    f=fields(rows,r)
    need(f=={'h':2195,'s':1,'q':53,'F':42153,'e':53,'S2':2673,'P2':672920},'current marked suffix')
    need(f['h']+4*f['F']==170807,'zero-ready component count')
    targets=[512423,512425];limit=max(targets);z0=(f['h'],f['s'],f['q'],f['F'])
    todo=[('',z0,rows,r,names)];states=[];hits=[];pruned=0;literal=[]
    while todo:
        word,z,g,root,ns=todo.pop(0)
        h,s,q,F=z
        if h+2*F>limit:pruned+=1;continue
        clocks=[];k=1
        while h+(1<<k)*F<=limit:
            value=h+(1<<k)*F;clocks.append([k,value])
            if value in targets:hits.append({'word':word,'profile':z,'clock':k,'H':value})
            k+=1
        states.append({'word':word,'profile':list(z),'root_name':ns[root],'component_names':ns,'clocks':clocks})
        if len(word)<=1 or word=='NR':
            lf=fields(g,root)
            need(tuple(lf[k] for k in ('h','s','q','F'))==z,'literal R/N profile')
            literal.append(word)
        for cmd in ('R','N'):
            nxt=transfer(z,cmd);need(nxt[0]+2*nxt[3]>h+2*F,'positive pruning functional')
            new=born(g,root);n=len(g)
            todo.append((word+cmd,nxt,new,root if cmd=='R' else n,ns+['component_birth_'+str(n)]))
    need(not hits,'unexpected ternary return in R/N menu')
    # All-Y polarized-seam exclusion uses a proved lower bound for |Y|>=3.
    b=f['h']+2*f['F'];large_bound=f['h']+18*f['F'];need(large_bound>limit,'large-seam lower bound')
    lower_second=2*f['P2']-f['S2'];need(2*lower_second>limit,'reversed seam lower bound')
    clock_families=[]
    for kind,a,c in [('singleton',b,b),('TT2 first root',b,2*b),('TT2 second root',3*(f['h']+f['F'])-f['s'],4*f['F']+3*f['h'])]:
        vals=[];k=1
        while True:
            value=a+(1<<k)*c;vals.append([k,value]);need(value not in targets,'small-seam clock hit')
            if value>limit:break
            k+=1
        clock_families.append({'Y':kind,'constant':a,'dyadic_coefficient':c,'clocks_through_first_overshoot':vals})
    need(f['P2']-f['F']==630767 and 2*(f['P2']-f['F'])>limit,'one-start funnel obstruction')
    dial=three_birth_dial(rows,r,names)
    need(dial['minimum_start_count']==53,'actual dial loses all unique-start marks')
    six=six_fields(rows,r);need(six==[2195,53,2673,42153,52480,53],'actual six-scalar seed')
    base_dial=three_birth_dial([0],0,['base'])
    need([base_dial['wrapper_pairs'][0][k] for k in ['newest_root_value','first_root_value']]==[17,19],'singleton direct-return boundary')
    refs=[path,'05-knowledge/results/triangle_return_170807_20261103/ZERO_READY.md',
      '05-knowledge/results/polarized_path_calculus_20261103/CALCULUS.md',
      '05-knowledge/results/polarized_path_calculus_20261103/FUNNELS.md',
      '05-knowledge/results/adaptive_metagraph_reformulation_20261103/ADAPTIVE.md',
      '05-knowledge/results/apex_wedge_constructor_20261103/PROFILE.md',
      '05-knowledge/results/ordered_response_ten_20261103/NEUTRAL_CARRIES.md',
      '05-knowledge/results/session161_source_matching_20261103/matched_sources/CHAIN.md',
      '05-knowledge/results/session161_source_matching_20261103/matched_sources/FAN.md']
    print(json.dumps({'status':'PASS','scope':'One actual marked suffix; exact complete R/N plus transitive-prefix clock menu below512425; all-size one-seam and one-funnel exclusions use proved positive bounds. No arbitrary component-rewrite exclusion.',
      'component_rows':rows,'component_names':names,'physical_root':31,'component_fields':f,'targets':targets,
      'fixed_clock_inverse':[{'target':N,'equation':'h_prime+4F_prime=N','scalar_triple_h':3*f['h'],'required_F_numerator':N-3*f['h'],'required_F_denominator':4} for N in targets],
      'three_birth_unit_dial':dial,'singleton_direct_return_control':base_dial,
      'joint_lower_bound_hostile':joint_lower_bound_hostile(),
      'retained_mark_six_scalar_seed':six,
      'retained_mark_six_scalar_controls':six_controls(),
      'RN_states':states,'RN_pruned_children':pruned,'RN_hits':hits,'RN_literal_profile_recounts':literal,
      'seam_small_clock_families':clock_families,'seam_large_Y_lower_bound':large_bound,
      'seam_B_second_F_lower_bound':lower_second,'one_start_funnel_F_lower_bound':630767,
      'seam_controls':small_seam_controls(),
      'sources':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in refs},
      'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2,sort_keys=True))
if __name__=='__main__':main()
