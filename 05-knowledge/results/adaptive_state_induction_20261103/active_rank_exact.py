#!/usr/bin/env python3
"""Exact rational certificates for cap-two reachable query charts.
Modular pivots select candidates only. Every closure identity is then checked
over Q against every coordinate of the proved path-cover carrier.
"""
from pathlib import Path
import importlib.util,json
from fractions import Fraction
from math import gcd,lcm

path=Path(__file__).with_name('active_rank_probe.py')
spec=importlib.util.spec_from_file_location('rank_probe',path)
rp=importlib.util.module_from_spec(spec);spec.loader.exec_module(rp)


def need(ok,msg):
    if not ok:raise RuntimeError(msg)


def dot(a,b):return sum(x*y for x,y in zip(a,b))
def transpose(A):return list(map(list,zip(*A)))
def mm(A,B):return [[dot(row,col) for col in transpose(B)] for row in A]


def inverse(A):
    n=len(A);M=[[Fraction(x) for x in row]+[Fraction(i==j) for j in range(n)] for i,row in enumerate(A)]
    det=Fraction(1)
    for j in range(n):
        k=next((i for i in range(j,n) if M[i][j]),None);need(k is not None,'exact nonsingular pivot')
        if k!=j:M[j],M[k]=M[k],M[j];det=-det
        a=M[j][j];det*=a;M[j]=[x/a for x in M[j]]
        for i in range(n):
            if i!=j and M[i][j]:
                b=M[i][j];M[i]=[x-b*y for x,y in zip(M[i],M[j])]
    return [row[n:] for row in M],det


def pivots(rows,p):
    basis=rp.Basis(p);chosen=[]
    for i,row in enumerate(rows):
        if basis.add(tuple(int(x)%p for x in row),(i,)):chosen.append(i)
    return chosen


def row_apply(row,mat,n):
    out=[0]*n
    for i,x in enumerate(row):
        if x:
            for j,m in mat[i]:out[j]+=x*m
    return tuple(out)


def col_apply(col,mat):return tuple(sum(m*col[j] for j,m in terms) for terms in mat)


def exact_prefix(seed,u,trans,blocks,word):
    c=seed;row=u
    for a in word:
        c,mat=trans[c,a];row=row_apply(row,mat,len(blocks[c]))
    return c,row


def exact_suffix(c,word,trans,final):
    if not word:return final[c]
    dest,mat=trans[c,word[0]]
    return col_apply(exact_suffix(dest,word[1:],trans,final),mat)


def encode(x):
    if isinstance(x,Fraction):return x.numerator if x.denominator==1 else str(x)
    if isinstance(x,(list,tuple)):return [encode(y) for y in x]
    if isinstance(x,dict):return {str(k):encode(v) for k,v in x.items()}
    return x


def poly_trim(a):
    while len(a)>1 and a[-1]==0:a.pop()
    return a


def poly_mod(a,b,p):
    a=poly_trim([x%p for x in a]);b=poly_trim([x%p for x in b])
    while len(a)>=len(b) and a!=[0]:
        j=len(a)-len(b);m=a[-1]*pow(b[-1],-1,p)%p
        for i,x in enumerate(b):a[i+j]=(a[i+j]-m*x)%p
        poly_trim(a)
    return a


def poly_mul(a,b,p):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]=(out[i+j]+x*y)%p
    return poly_trim(out)


def poly_power(a,n,mod,p):
    out=[1]
    while n:
        if n&1:out=poly_mod(poly_mul(out,a,p),mod,p)
        a=poly_mod(poly_mul(a,a,p),mod,p);n//=2
    return out


def poly_gcd(a,b,p):
    while b!=[0]:a,b=b,poly_mod(a,b,p)
    c=pow(a[-1],-1,p)
    return [(x*c)%p for x in a]


def output_annihilator(M,coeff):
    n=len(M);v=[Fraction(i==0) for i in range(n)];total=[Fraction(0)]*n
    for a in coeff:
        total=[x+a*y for x,y in zip(total,v)]
        v=[dot(row,v) for row in M]
    return total


def main():
    d=2;p=1000003;proposal=rp.probe(d,p)
    seed,blocks,trans,reverse,u,final=rp.build(d)
    data={tuple(b['control']):b for b in proposal['blocks']};R={};O={};Rinv={};Rcols={};queries={};Q={};Hinv={};Hrows={}
    entries=[];closure_checks=0;query_checks=0
    for c,b in data.items():
        prefixes=[tuple(w) for w in b['reachable_prefixes']]
        R[c]=[exact_prefix(seed,u,trans,blocks,w)[1] for w in prefixes]
        need(all(exact_prefix(seed,u,trans,blocks,w)[0]==c for w in prefixes),'actual prefix control')
        cols=pivots(transpose(R[c]),p);need(len(cols)==len(R[c]),'reachable pivot number')
        Rcols[c]=cols;Rinv[c],rdet=inverse([[row[j] for j in cols] for row in R[c]])
        suffixes=[tuple(w) for w in b['observable_suffixes']]
        O[c]=[exact_suffix(c,w,trans,final) for w in suffixes]
        all_H=mm(R[c],transpose(O[c]));chosen_rows=pivots(all_H,p)
        chosen_cols=pivots(transpose([all_H[i] for i in chosen_rows]),p)
        h=len(chosen_rows);need(h==len(chosen_cols)==b['local_hankel_mod_rank'],'Hankel pivots')
        queries[c]=[suffixes[j] for j in chosen_cols]
        O[c]=[O[c][j] for j in chosen_cols]
        Q[c]=mm(R[c],transpose(O[c]));Hrows[c]=chosen_rows
        square=[Q[c][i] for i in chosen_rows];Hinv[c],hdet=inverse(square)
        need(queries[c][0]==() and all(row[0]==dot(r,final[c]) for row,r in zip(Q[c],R[c])),'H is first query')
        entries.append({'control':c,'carrier_coordinates':len(blocks[c]),'exact_reachable_dimension':len(R[c]),
                        'exact_query_dimension':h,'reachable_prefixes':prefixes,'reachable_coordinate_pivots':cols,
                        'reachable_pivot_determinant':rdet,'query_suffixes':queries[c],
                        'hankel_prefixes':[prefixes[i] for i in chosen_rows],
                        'hankel_matrix':square,'hankel_determinant':hdet})
    transfer=[]
    for c in blocks:
        for a in range(len(c)):
            dest,mat=trans[c,a];images=[row_apply(row,mat,len(blocks[dest])) for row in R[c]]
            for image in images:
                coefficients=mm([[image[j] for j in Rcols[dest]]],Rinv[dest])[0]
                need(mm([coefficients],R[dest])[0]==list(image),'exact all-coordinate reachable closure')
                closure_checks+=1
            values=mm(images,transpose(O[dest]))
            small=mm(Hinv[c],[values[i] for i in Hrows[c]])
            need(mm(Q[c],small)==values,'exact query closure on whole reachable space')
            query_checks+=len(R[c])*len(O[dest])
            transfer.append({'source':c,'root_rank':a,'target':dest,'matrix':small,
                             'denominator_lcm':lcm(*(x.denominator for row in small for x in row)),
                             'negative_entries':sum(x<0 for row in small for x in row)})
    initial=[dot(u,col) for col in O[seed]]
    stable=(2,2,1)
    returns={t['root_rank']:t['matrix'] for t in transfer if t['source']==stable and t['target']==stable}
    poly0=[2,-3,1];poly1=[1,-4,2,-6,3,-4,1,-2,1]
    need(all(x==0 for x in output_annihilator(returns[0],poly0)),'root0 exact scalar recurrence')
    need(all(x==0 for x in output_annihilator(returns[1],poly1)),'root1 exact scalar recurrence')
    frobenius=[0,1];rabin4=None
    for k in range(1,9):
        frobenius=poly_power(frobenius,31,poly1,31)
        if k==4:
            difference=frobenius[:]+[0]*max(0,2-len(frobenius));difference[1]=(difference[1]-1)%31
            rabin4=poly_gcd(poly1,poly_trim(difference),31)
    need(frobenius==[0,1] and rabin4==[1],'degree8 irreducible modulo31 certificate')
    c,row=exact_prefix(seed,u,trans,blocks,(1,0));need(c==stable,'stable actual seed')
    y=[dot(row,col) for col in O[stable]];tail=[]
    for _ in range(12):tail.append(y[0]);y=mm([y],returns[1])[0]
    # Sharp global dimension for the word series extended by zero on illegal words.
    offsets={};dimension=0
    for c in blocks:offsets[c]=dimension;dimension+=len(O[c])
    global_matrices=[[[Fraction(0)]*dimension for _ in range(dimension)] for _ in range(4)]
    global_final=[Fraction(0)]*dimension
    for c in blocks:global_final[offsets[c]]=1
    for t in transfer:
        i=offsets[t['source']];j=offsets[t['target']];a=t['root_rank']
        for r,values in enumerate(t['matrix']):
            for s,value in enumerate(values):global_matrices[a][i+r][j+s]=value
    observable=rp.Basis(p);observable.add(tuple(int(v) for v in global_final),())
    global_columns=[global_final];global_suffixes=[()];todo=[(global_final,())]
    while todo:
        col,w=todo.pop(0)
        for a,M in enumerate(global_matrices):
            target=[dot(row,col) for row in M]
            residues=tuple(v.numerator*pow(v.denominator,-1,p)%p for v in target)
            if observable.add(residues,(a,)+w):
                global_columns.append(target);global_suffixes.append((a,)+w);todo.append((target,(a,)+w))
    need(len(global_columns)==dimension,'global observable candidate dimension')
    global_rows=[];global_prefixes=[]
    for c in blocks:
        for i in Hrows[c]:
            row=[Fraction(0)]*dimension;row[offsets[c]:offsets[c]+len(O[c])]=Q[c][i]
            global_rows.append(row);global_prefixes.append(tuple(data[c]['reachable_prefixes'][i]))
    global_H=mm(global_rows,transpose(global_columns))
    need(all(Fraction(z).denominator==1 for row in global_H for z in row),'global actual Hankel integrality')
    _,global_det=inverse(global_H)
    return encode({'status':'PROVED BY EXACT FINITE LINEAR CERTIFICATE; control-dependent query dimensions sharp',
                   'source_sha256':proposal['source_sha256'],'alphabet':'active root rank, first C3 root fixed0',
                   'scope':'actual C3-seed cap2 language; rational query charts; global minimum28 only with illegal words assigned0; no positive-matrix minimum claim',
                   'total_reachable_dimensions':sum(len(v) for v in R.values()),
                   'total_query_chart_dimensions':sum(len(v) for v in O.values()),
                   'all_coordinate_closure_checks':closure_checks,'scalar_query_closure_checks':query_checks,
                   'initial_control':seed,'initial_query_vector':initial,'output':'first coordinate',
                   'charts':entries,'transfers':transfer,
                   'global_guarded_hankel':{'dimension':dimension,'prefixes':global_prefixes,
                      'suffixes':global_suffixes,'matrix':global_H,'determinant':global_det,
                      'guard':'illegal active-rank words have value0'},
                   'stable_return_control':{'control':stable,'root0_annihilator_ascending':poly0,
                     'root1_annihilator_ascending':poly1,'irreducible_modulus':31,
                     'rabin_frobenius8':frobenius,'rabin_gcd4':rabin4,
                     'actual_seed_active_rank_prefix':[1,0],'root1_tail_first12':tail}})


if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
