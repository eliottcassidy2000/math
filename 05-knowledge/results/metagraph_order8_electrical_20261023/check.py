"""Sparse, exact dyadic filtration of a colored line graph.

Primary graph input comes from independently checked actual tiling exports.
Small controls use the earlier independently enumerated colored graphs.
"""
from pathlib import Path
from hashlib import sha256
import csv
import gzip
import json
import argparse

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
CHECKS=0


def need(ok, why):
    global CHECKS
    CHECKS+=1
    if not ok:raise ValueError(why)


def kernel(rows, n):
    piv={}
    for row in rows:
        while row:
            p=(row&-row).bit_length()-1
            if p in piv:row^=piv[p]
            else:
                piv[p]=row
                break
    basis=[]
    for f in range(n):
        if f in piv:continue
        x=1<<f
        for p in sorted(piv,reverse=True):
            if (piv[p]&x).bit_count()%2:x^=1<<p
        need(all((row&x).bit_count()%2==0 for row in rows),'binary kernel')
        basis.append(x)
    return basis,len(piv)


def solver(rows,n):
    piv={};zero=[]
    for i,row in enumerate(rows):
        move=1<<i
        while row:
            p=(row&-row).bit_length()-1
            if p in piv:
                row^=piv[p][0];move^=piv[p][1]
            else:
                piv[p]=(row,move)
                break
        if not row:zero.append(move)
    def solve(rhs):
        need(all((z&rhs).bit_count()%2==0 for z in zero),'consistent lift RHS')
        x=0
        for p in sorted(piv,reverse=True):
            row,move=piv[p]
            if ((move&rhs).bit_count()+(row&x).bit_count())%2:x^=1<<p
        need(all(((row&x).bit_count()%2)==((rhs>>i)&1) for i,row in enumerate(rows)),'binary lift solution')
        return x
    return solve


def components(n,edges):
    parent=list(range(n))
    def root(x):
        while parent[x]!=x:
            parent[x]=parent[parent[x]];x=parent[x]
        return x
    for i,j,w in edges:
        if w:parent[root(i)]=root(j)
    return len({root(i) for i in range(n)})


def apply(n,edges,x):
    y=[0]*n
    for i,j,w in edges:
        d=w*(x[i]-x[j]);y[i]+=d;y[j]-=d
    return y


def calculate(n,L,B):
    """L blue weights, B paired-black weights; loops harmless and retained."""
    G=L+[(i,j,2*w) for i,j,w in B]
    need(components(n,G)==1,'connected full graph')
    active=sorted({i for e in L if e[0]!=e[1] and e[2] for i in e[:2]})
    need(bool(active),'nonempty nonloop blue support')
    root=active[-1];nonroot=active[:-1];ix={v:i for i,v in enumerate(nonroot)}
    isolated=sorted(set(range(n))-set(active));a=len(nonroot)
    blue=[[0]*a for _ in range(a)]
    for i,j,w in L:
        if i==j:continue
        if i!=root:blue[ix[i]][ix[i]]+=w
        if j!=root:blue[ix[j]][ix[j]]+=w
        if i!=root and j!=root:
            blue[ix[i]][ix[j]]-=w;blue[ix[j]][ix[i]]-=w
    binary=[sum((v%2)<<j for j,v in enumerate(row)) for row in blue]
    kb,rankblue=kernel(binary,a);solve=solver(binary,a)
    base=len(isolated);d=base+len(kb)
    labels=[0]*n
    for i,v in enumerate(isolated):labels[v]=1<<i
    for t,b in enumerate(kb):
        for k,v in enumerate(nonroot):
            if b>>k&1:labels[v]|=1<<(base+t)
    Q=[0]*d
    for i,j,w in B:
        if w%2:
            z=labels[i]^labels[j];scan=z
            while scan:
                bit=scan&-scan;scan-=bit;Q[bit.bit_length()-1]^=z
    # The blue contribution is an INTEGER half-energy, not a binary diagonal.
    for t,u in enumerate(kb):
        for s,v in enumerate(kb):
            energy=sum(blue[i][j] for i in range(a) if u>>i&1 for j in range(a) if v>>j&1)
            need(energy%2==0,'blue harmonic energy even')
            if (energy//2)%2:Q[base+t]^=1<<(base+s)
    need(all(((Q[i]>>j)&1)==((Q[j]>>i)&1) for i in range(d) for j in range(i)), 'symmetric half-energy')
    rad,rankQ=kernel(Q,d)
    # Lift each radical vector to a literal integer potential with Lap*x=0 mod4.
    lifts=[];images=[]
    for z in rad:
        x=[(label&z).bit_count()%2 for label in labels]
        image=apply(n,G,x)
        need(all(v%2==0 for v in image),'first harmonic lift')
        need(all((image[v]//2)%2==0 for v in isolated),'radical at isolated rows')
        rhs=sum(((image[v]//2)%2)<<i for i,v in enumerate(nonroot))
        correction=solve(rhs)
        for i,v in enumerate(nonroot):x[v]+=2*((correction>>i)&1)
        image=apply(n,G,x)
        need(x[root]==0 and all(v%4==0 for v in image),'literal lift modulo4')
        lifts.append(x);images.append(image)
    R=[]
    for x in lifts:
        row=0
        for j,y in enumerate(images):
            value=sum(a*b for a,b in zip(x,y))
            need(value%4==0,'quarter energy integral')
            row|=((value//4)%2)<<j
        R.append(row)
    need(all(((R[i]>>j)&1)==((R[j]>>i)&1) for i in range(len(R)) for j in range(i)), 'symmetric quarter-energy')
    radr,rankR=kernel(R,len(R))
    # A radical quarter-mode can need corrections from TWO earlier layers.
    # First use the full half-form, then solve the remaining blue mod-two system.
    eighth_lifts=[];eighth_images=[]
    if radr:
        solveQ=solver(Q,d)
        for z in radr:
            x=[sum(lifts[t][i] for t in range(len(lifts)) if z>>t&1)%4 for i in range(n)]
            image=apply(n,G,x)
            need(all(v%4==0 for v in image),'quarter radical base lift')
            rhs=0
            for i,v in enumerate(image):
                if (v//4)%2:rhs^=labels[i]
            coefficients=solveQ(rhs)
            for i in range(n):x[i]+=2*((labels[i]&coefficients).bit_count()%2)
            image=apply(n,G,x)
            need(all(v%4==0 for v in image),'half-form correction retains mod4')
            need(all((image[v]//4)%2==0 for v in isolated),'second lift at isolated rows')
            rhs=sum(((image[v]//4)%2)<<i for i,v in enumerate(nonroot))
            correction=solve(rhs)
            for i,v in enumerate(nonroot):x[v]+=4*((correction>>i)&1)
            image=apply(n,G,x)
            need(x[root]==0 and all(v%8==0 for v in image),'literal lift modulo8')
            eighth_lifts.append(x);eighth_images.append(image)
    S=[]
    for x in eighth_lifts:
        row=0
        for j,y in enumerate(eighth_images):
            value=sum(a*b for a,b in zip(x,y));need(value%8==0,'eighth energy integral')
            row|=((value//8)%2)<<j
        S.append(row)
    rads,rankS=kernel(S,len(S))
    forms=[Q,R,S];bases=[None,lifts,eighth_lifts]
    layer_solvers={1:solveQ} if radr else {}
    higher=[];current_rad=rads
    for level in range(3,13):
        if not current_rad:break
        modulus=1<<level;new_lifts=[];new_images=[]
        for z in current_rad:
            x=[sum(bases[-1][t][i] for t in range(len(bases[-1])) if z>>t&1)%modulus for i in range(n)]
            for k in range(level-1,0,-1):
                image=apply(n,G,x)
                need(all(v%modulus==0 for v in image),'earlier-layer correction preserves divisibility')
                b=[v//modulus for v in image]
                if k==1:
                    rhs=0
                    for i,v in enumerate(b):
                        if v%2:rhs^=labels[i]
                else:
                    rhs=sum((sum(a*v for a,v in zip(u,b))%2)<<t for t,u in enumerate(bases[k-1]))
                if k not in layer_solvers:layer_solvers[k]=solver(forms[k-1],len(forms[k-1]))
                co=layer_solvers[k](rhs)
                for i in range(n):
                    correction=((labels[i]&co).bit_count()%2 if k==1 else
                                sum(bases[k-1][t][i] for t in range(len(bases[k-1])) if co>>t&1))
                    x[i]+=(1<<(level-k))*correction
            image=apply(n,G,x)
            need(all(v%modulus==0 for v in image),'all previous form corrections')
            need(all((image[v]//modulus)%2==0 for v in isolated),'last layer isolated rows')
            rhs=sum(((image[v]//modulus)%2)<<i for i,v in enumerate(nonroot))
            correction=solve(rhs)
            for i,v in enumerate(nonroot):x[v]+=modulus*((correction>>i)&1)
            x=[v%(2*modulus) for v in x];image=apply(n,G,x)
            need(x[root]==0 and all(v%(2*modulus)==0 for v in image),'higher literal lift')
            new_lifts.append(x);new_images.append(image)
        newform=[]
        for x in new_lifts:
            row=0
            for t,y in enumerate(new_images):
                value=sum(a*b for a,b in zip(x,y));need(value%(2*modulus)==0,'higher form integral')
                row|=((value//(2*modulus))%2)<<t
            newform.append(row)
        newrad,newrank=kernel(newform,len(newform))
        higher.append({'level':level+1,'dimension':len(newform),'rank':newrank,
                       'next_dimension':len(newrad),'form_rows':newform,'lifts':new_lifts})
        forms.append(newform);bases.append(new_lifts);current_rad=newrad
    # Direct pairing controls avoid a self-confirming symbolic construction.
    trials=range(min(d,20))
    for t in trials:
        u=[(z>>t)&1 for z in labels]
        image=apply(n,G,u)
        for s in trials:
            value=sum(y*((labels[i]>>s)&1) for i,y in enumerate(image))
            need(value%2==0 and (value//2)%2==((Q[t]>>s)&1),'direct half-energy control')
    return {'vertices':n,'root':root,'blue_components':components(n,L),
            'blue_active_vertices':len(active),'blue_reduced_rank_mod2':rankblue,
            'blue_bicycle_dimension':len(kb)-(components(n,L)-len(isolated)-1),
            'critical_two_rank':d,'half_form_rank':rankQ,'four_divisible_factors':len(rad),
            'quarter_form_rank':rankR,'eight_divisible_factors':len(radr),
            'half_form_rows_sha256':sha256(json.dumps(Q).encode()).hexdigest(),
            'quarter_form_rows':R,'radical_lifts_mod4':lifts,
            'eighth_form_rows':S,'eighth_radical_lifts_mod8':eighth_lifts,
            'sixteen_divisible_factors':len(rads),'higher_layers':higher,'complete_two_primary':not current_rad,
            'exact_valuation_one':rankQ,'exact_valuation_two':rankR,'exact_valuation_three':rankS}


def raw_read(path):
    if path.exists():raw=path.read_bytes()
    else:raw=gzip.decompress(Path(str(path)+'.gz').read_bytes())
    return raw.replace(b'\r\n',b'\n')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    old=ROOT/'05-knowledge/results/session36_information_curl_contexts_20261023/metagraph/dyadic_network_receipt.json'
    oldraw=raw_read(old);controls=[]
    for case in json.loads(oldraw):
        if case['n']>7:continue
        out=calculate(case['nodes'],case['blue_edges'],case['paired_black_edges'])
        need(out['critical_two_rank']==case['critical_two_rank'],'inherited 2-rank control')
        need(out['eight_divisible_factors']==0,'known small exponent control')
        controls.append({'n':case['n'],**out})
    folder=HERE.parent/'metagraph_order8_20261023'
    nr=raw_read(folder/'merged_nodes.tsv');er=raw_read(folder/'merged_edges.tsv')
    nodes=list(csv.DictReader(nr.decode().splitlines(),delimiter='\t'));L=[];B=[]
    for row in csv.DictReader(er.decode().splitlines(),delimiter='\t'):
        i,j,blue,black=(int(row[k]) for k in ('source','target','blue','black'))
        need(0<=i<=j<len(nodes) and blue>=0 and black>=0 and black%2==0,'literal colored edge')
        if blue:L.append((i,j,blue))
        if black:B.append((i,j,black//2))
    out=calculate(len(nodes),L,B)
    result={'status':'FINITE-EXACT graph export consumer; general forms separately proved',
            'input_sha256_lf':{'nodes':sha256(nr).hexdigest(),'edges':sha256(er).hexdigest(),'small_graphs':sha256(oldraw).hexdigest()},
            'controls':controls,'order8':out,'checks':CHECKS}
    raw=(json.dumps(result,indent=2,sort_keys=True)+'\n').encode()
    if args.write:(HERE/'receipt.json').write_bytes(raw)
    else:need(raw_read(HERE/'receipt.json')==raw,'saved receipt')
    summary={k:v for k,v in out.items() if k not in ('radical_lifts_mod4','quarter_form_rows','eighth_radical_lifts_mod8','higher_layers')}
    summary['higher_layers']=[{k:v for k,v in row.items() if k!='lifts'} for row in out['higher_layers']]
    print(json.dumps({'checks':CHECKS,'order8':summary,'receipt_sha256_lf':sha256(raw).hexdigest()}))


if __name__=='__main__':main()
