from audit import *
from subprocess import run,check_output
Q=(238,4,240,118,195,146,34,74)
need(endpoints(Q)==[0]*8,'Q8 zero endpoint parity')
need(delta(Q,0)==0,'Q8 complete coarse zero polynomial')
records=[]; requests=[];actual_meta=[]
for j in range(8):
    Z,B=inflate(Q,[tt(2) if i==j else tt(1) for i in range(8)]);u,b=B[j];n=9
    unchanged,_=inflate(op(Q),[tt(2) if i==j else tt(1) for i in range(8)])
    for k in [2,3,4,5]:
        w=1<<u|((k==2)<<b)
        for target,Zp in [('unchanged',unchanged),('reversed',op(Z))]:
            d=weight_poly(cycle_basis(Z),w)^weight_poly(cycle_basis(Zp),w)
            need(d!=0,'Q8 nonzero')
            if target=='reversed' or k==2:value=L(Z,d,u)
            else:value=L(Z,d,u)^L(Z,d,b)^((d>>((1<<u)|(1<<b)|(1<<(n+b))))&1)
            need(value==1,'Q8 exact functional')
            firstbit=d&-d;coefficient_mask=firstbit.bit_length()-1
            sub=coefficient_mask
            while not evaluate(d,sub):
                need(sub!=0,'nonzero coefficient must have query witness');sub=(sub-1)&coefficient_mask
            row={'j':j,'k':k,'target':target,'functional':value,'first_nonzero_monomial':coefficient_mask,'kernel_query':sub}
            records.append(row)
            if j==0 and k in [2,3]:
                A=[tt(k) if i==j else tt(2) for i in range(8)];C=list(A)
                if target=='reversed':C[j]=op(C[j])
                W,blocks=inflate(Q,A);U,_=inflate(op(Q),C)
                groups=[]
                for i,block in enumerate(blocks):groups.extend([[block[0]],block[1:]] if i==j else [block])
                a=sum(1<<v for i,g in enumerate(groups) if sub>>i&1 for v in g)
                c=sum(1<<v for i,g in enumerate(groups) if sub>>(n+i)&1 for v in g)
                actual_meta.append({'j':j,'k':k,'target':target,'kind':'local','source':W,'target_rows':U,'a':a,'c':c})
                for rows in [W,U]:requests.append(ports(rows,a,c))
                if k==2:
                    AQ=(2,4,1);AW,_=inflate(AQ,[W,tt(1),tt(1)]);AU,_=inflate(AQ,[U,tt(1),tt(1)]);fixed=1<<(len(W)+1)
                    actual_meta.append({'j':j,'k':k,'target':target,'kind':'ambient','source':AW,'target_rows':AU,'a':a|fixed,'c':c|fixed})
                    for rows in [AW,AU]:requests.append(ports(rows,a|fixed,c|fixed))
    print('Q8 vertex',j,'passed',flush=True)
# Only orders <=20: each exact Hamilton count <=20!<2^64.
need(max(map(len,requests))<=20,'uint64 safe order')
run(['g++','-O3','-std=c++17',str(P/'physical_count.cpp'),'-o',str(P/'physical_count')],check=True)
request_text='\n'.join(' '.join(map(str,[len(r),*r])) for r in requests)+'\n'
(P/'q8_physical_input.txt').write_text(request_text)
values=list(map(int,check_output([str(P/'physical_count')],input=request_text,text=True).split()))
need(len(values)==len(requests),'physical result count')
for i,row in enumerate(actual_meta):
    row['source_count'],row['target_count']=values[2*i:2*i+2];row['difference']=row['source_count']-row['target_count'];need(row['difference']%4==2,'Q8 literal separator')
    print('actual',row['kind'],row['k'],row['target'],row['source_count'],row['target_count'],row['difference'],flush=True)
out={'status':'PASS','Q8':Q,'endpoint_parity':[0]*8,'complete_coarse_polynomial_zero':True,'kernel_records':records,'physical_records':actual_meta,'maximum_augmented_order':max(map(len,requests))}
(P/'q8_controls.json').write_text(json.dumps(out,indent=2)+'\n')
