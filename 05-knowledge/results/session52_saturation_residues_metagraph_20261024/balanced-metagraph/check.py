#!/usr/bin/env python3
"""Independent finite controls of degree-coefficient balanced port identities."""
from collections import Counter
from itertools import combinations,permutations
from pathlib import Path
from hashlib import sha256
import argparse,json
from engine import need,summary,coefficient_polynomial,eulerian_anticycle_dp
HERE=Path(__file__).resolve().parent

def rows(n,word,pairs):
    T=[0]*n
    for i in range(n-1):T[i]|=1<<(i+1)
    for k,(i,j) in enumerate(pairs):
        a,b=(i,j) if word>>k&1 else (j,i);T[a]|=1<<b
    return tuple(T)

def normalize(T,path):
    return tuple(sum(((T[u]>>v)&1)<<j for j,v in enumerate(path)) for u in path)

def literal_cube(n):
    p=n//2;pairs=[(i,j) for i in range(n) for j in range(i+2,n)];m=len(pairs);mask=(1<<m)-1
    pairindex={ij:k for k,ij in enumerate(pairs)}
    reflected=[pairindex[n-1-j,n-1-i] for i,j in pairs]
    endpoint=Counter();blue=Counter();count=Counter();digest=sha256()
    d=list(map(int,[i<n-1 for i in range(n)]))
    for i,j in pairs:d[j]+=1
    oldword=0;balanced_words=[];degrees={};rotation_counts=Counter()
    for index in range(1<<m):
        word=index^(index>>1)
        if index:
            changed=word^oldword;k=changed.bit_length()-1;i,j=pairs[k]
            sign=1 if word>>k&1 else -1;d[i]+=sign;d[j]-=sign
        oldword=word
        if any(x not in (p-1,p) for x in d):continue
        tag=('H' if d[0]==p else 'L')+('H' if d[-1]==p else 'L')
        endpoint[tag]+=1;T=rows(n,word,pairs);U=rows(n,word^mask,pairs)
        need([row.bit_count() for row in T]==d,'Gray scores independently equal actual matrix scores')
        dC=[row.bit_count() for row in U]
        need(dC==[n-1-d[v]+int(v==0)-int(v==n-1) for v in range(n)],'actual C endpoint identity')
        internal=all(x in (p-1,p) for x in dC)
        need(internal==(tag=='HL'),'actual C stays balanced iff downhill')
        R=sum(((word>>k)&1)<<reflected[k] for k in range(m))
        isblue=(word==R)
        q=tuple(reversed(range(n)));revconverse=tuple(((1<<n)-1)^r^(1<<i) for i,r in enumerate(normalize(T,q)))
        need(rows(n,R,pairs)==revconverse,'R really is reversed converse')
        if isblue:
            blue[tag]+=1;need(tag in ('HL','LH'),'blue endpoints are opposite')
            need(rows(n,(word^mask),pairs)==normalize(tuple(((1<<n)-1)^r^(1<<i) for i,r in enumerate(U)),q),'C preserves R-fixedness')
        if internal:
            count['internal_ports']+=1
            if isblue:count['internal_blue_ports']+=1
        else:
            count['boundary_ports']+=1
        if n>=4:need(word!=(word^mask),'C free except n2')
        # Unique regular extension and fixed-cycle transport.
        Rext=list(T)+[0];z=n
        for v,a in enumerate(d):
            if a==p:Rext[z]|=1<<v
            else:Rext[v]|=1<<z
        need(all(r.bit_count()==p for r in Rext),'unique extension regular')
        if internal:
            need(Rext[z]&1 and Rext[n-1]>>z&1,'downhill path closes at extension vertex')
            Uext=list(U)+[0]
            for v,a in enumerate(dC):
                if a==p:Uext[z]|=1<<v
                else:Uext[v]|=1<<z
            Cext=[((1<<(n+1))-1)^r^(1<<i) for i,r in enumerate(Rext)]
            cyc=(z,)+tuple(range(n))+(z,)
            for a,b in zip(cyc,cyc[1:]):Cext[a]|=1<<b;Cext[b]&=~(1<<a)
            need(Cext==Uext,'lift intertwines C with noncycle reversal')
            for shift in range(1,n+1):
                if (n+1)%shift:continue
                q=tuple((v+shift)%(n+1) for v in range(n+1))
                if all(((Rext[a]>>b)&1)==((Rext[q[a]]>>q[b])&1) for a in range(n+1) for b in range(n+1)):
                    rotation_counts[shift]+=1
                    if shift==1:need(isblue,'every one-step circulant fixed orientation is blue')
        digest.update(f'{word}:{tag}:{int(isblue)}:{int(internal)}\n'.encode())
        if n<=6:balanced_words.append((word,T));degrees[word]=(tag,isblue)
    result={'n':n,'cube_words':1<<m,'balanced_ports':sum(endpoint.values()),'endpoint_ports':{k:endpoint[k] for k in ('LL','LH','HL','HH')},
            'balanced_blue_ports':sum(blue.values()),'blue_endpoint_ports':dict(sorted(blue.items())),**count,'balanced_word_digest':digest.hexdigest(),'proper_rotation_fixed_counts':dict(sorted(rotation_counts.items()))}
    if n<=6:
        allperms=list(permutations(range(n)))
        canoncache={}
        def canonical(T):
            if T not in canoncache:canoncache[T]=min(normalize(T,p) for p in allperms)
            return canoncache[T]
        classes={word:canonical(T) for word,T in balanced_words}
        source=Counter(classes.values());intr=Counter();bl=Counter();loops=Counter();ordinary_lines=set();merged_lines=set();merged_class={}
        full=(1<<n)-1
        for key in source:
            op=tuple(full^r^(1<<i) for i,r in enumerate(key))
            merged_class[key]=min(key,canonical(op))
        for word,T in balanced_words:
            tag,isblue=degrees[word];key=classes[word]
            if tag=='HL':
                intr[key]+=1;bl[key]+=int(isblue)
                if n>=4:
                    line=min(word,word^mask)
                    if key==classes[word^mask]:ordinary_lines.add(line)
                    if merged_class[key]==merged_class[classes[word^mask]]:merged_lines.add(line)
        classreports=[]
        for key,mass in sorted(source.items()):
            aut=sum(normalize(key,p)==key for p in allperms)
            hp=[p for p in allperms if all(key[a]>>b&1 for a,b in zip(p,p[1:]))]
            hdown=sum(key[p[0]].bit_count()==n//2 and key[p[-1]].bit_count()==n//2-1 for p in hp)
            need(len(hp)==mass*aut,'actual H/Aut agrees with class word mass')
            need(hdown==intr[key]*aut,'downhill actual H/Aut agrees with internal class word degree')
            ext=list(key)+[0]
            for v in range(n):
                if key[v].bit_count()==n//2:ext[n]|=1<<v
                else:ext[v]|=1<<n
            cycles=sum(all(ext[a]>>b&1 for a,b in zip((n,)+p,p+(n,))) for p in allperms)
            need(cycles==hdown,'actual regular-parent Hamilton cycles match downhill paths')
            classreports.append({'rows':key,'H':len(hp),'aut':aut,'word_mass':mass,'downhill_H':hdown,'internal_word_degree':intr[key],'internal_blue_degree':bl[key],'regular_parent_Hamilton_cycles':cycles})
        result['class_control']={'ordinary_classes':len(source),'merged_classes':len(set(merged_class.values())),'ordinary_internal_loop_lines':len(ordinary_lines),'merged_internal_loop_lines':len(merged_lines),'classes':classreports}
    return result

def phi(n):
    return sum(__import__('math').gcd(i,n)==1 for i in range(1,n+1))

def rotation_fixed(n,d):
    N=n+1;need(N%d==0,'rotation orbit divisor');L=N//d
    if d==1:return 1<<(n//2-1)
    target=(d-1)*L//2-1
    factors=[]
    for i in range(d):
        for j in range(i+1,d):
            multiplicity=L-int(j==i+1 or (i==0 and j==d-1))
            factors += [(i,j)]*multiplicity
    poly={(0,)*d:1}
    for i,j in factors:
        nxt=Counter()
        for a,c in poly.items():
            for v in (i,j):
                if a[v]<target:
                    b=list(a);b[v]+=1;nxt[tuple(b)]+=c
        poly=nxt
    return (1<<(d*(L-1)//2))*poly.get((target,)*d,0)

def main():
    a=argparse.ArgumentParser();a.add_argument('--output',default=str(HERE/'normal.json'));a.add_argument('--max-order',type=int,default=12);args=a.parse_args()
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
    for pin in pins['sources']:need(sha256((HERE/pin['path']).read_bytes()).hexdigest()==pin['sha256'],'source pin '+pin['path'])
    values=[];literal=[]
    for n in range(2,args.max_order+1,2):
        val=summary(n)
        eo,states=eulerian_anticycle_dp(n);need(eo==val['internal_ports'],'independent anticycle coefficient DP')
        val['anticycle_dp_states']=states
        if n<=8:
            product,peak=coefficient_polynomial(n);need(product==val['endpoint_ports'],'independent complete degree product')
            val['ordinary_polynomial_peak_states']=peak
            lit=literal_cube(n);literal.append(lit)
            for key in ('balanced_ports','endpoint_ports','balanced_blue_ports','internal_ports','internal_blue_ports','boundary_ports'):need(val[key]==lit.get(key,0),'literal cube agreement '+key)
        proper={d:rotation_fixed(n,d) for d in range(1,n+1) if (n+1)%d==0}
        numerator=eo+sum(phi((n+1)//d)*v for d,v in proper.items())
        need(numerator%(n+1)==0,'rotation Burnside integrality')
        val['proper_rotation_fixed_counts']=proper;val['cyclic_orbits_internal_words']=numerator//(n+1)
        if n<=8:need(proper==literal[-1]['proper_rotation_fixed_counts'],'literal rotation fixed-count polynomial')
        need(val['internal_blue_ports']>=1<<(n//2-1),'circulant blue lower bound')
        need(val['endpoint_ports']['LL']==val['endpoint_ports']['HH'],'reversed-converse endpoint symmetry')
        if n>=4:need((val['internal_ports']-val['internal_blue_ports'])%4==0,'free internal black Klein orbits')
        values.append(val)
    inherited=json.loads((HERE/'sources/regular-eleven-counts.json').read_text())
    M=inherited['almost_regular_parent_classes'];SC=inherited['self_converse_parents']
    need((M,SC)==(13333,229),'pinned inherited ordinary/self-converse balanced n10 counts')
    inherited_table={'source':'sources/regular-eleven-counts.json','ordinary_balanced_classes':M,'self_converse_balanced_classes':SC,'merged_balanced_nodes':(M+SC)//2,'non_self_converse_merged_pairs':(M-SC)//2,'status':'INHERITED finite-exact, not inferred from port polynomials'}
    result={'inherited_n10_class_counts':inherited_table,'status':'PASS','scope':'n-only DP even orders through max_order; independent anticycle DP; literal full fixed-path cubes and direct score polynomial through n8; H/Aut and regular-parent cycle controls through n6','max_order':args.max_order,'counts':values,'literal_controls':literal}
    Path(args.output).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','max_order':args.max_order,'counts':values},indent=2,sort_keys=True))
if __name__=='__main__':main()
