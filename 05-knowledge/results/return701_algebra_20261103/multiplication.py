"""Complete sparse Q-algebra table, radical powers and a nilpotent clock."""
import json
import sympy as s
from common import HERE,need,save,data,rref

def main():
    d=json.loads((HERE/'algebra_exact.json').read_text());R=s.Matrix(d['basis']);P=d['coordinate_pivots'];n=27
    rowmaps=[];colmaps=[]
    for i in range(95):
        M=s.Matrix(n,n,list(R[i,:]));rowmaps.append([{j:M[a,j] for j in range(n) if M[a,j]} for a in range(n)])
        colmaps.append([{a:M[a,j] for a in range(n) if M[a,j]} for j in range(n)])
    table={};entries=[]
    for i in range(95):
        for j in range(95):
            result={}
            for k,p in enumerate(P):
                a,b=divmod(p,n);v=sum(c*colmaps[j][b].get(t,0) for t,c in rowmaps[i][a].items())
                if v:result[k]=v
            if result:
                table[i,j]=result;entries.append([i,j,[[k,int(c) if c.q==1 else str(c)] for k,c in result.items()]])
    def sparse(row):return {i:c for i,c in enumerate(row) if c}
    def product(a,b):
        out={}
        for i,x in a.items():
            for j,y in b.items():
                for k,z in table.get((i,j),{}).items():out[k]=out.get(k,0)+x*y*z
        return [out.get(i,0) for i in range(95)]
    K=s.Matrix(d['trace_kernel_coordinates']);kb=[sparse(K[i,:]) for i in range(K.rows)];W=K;dims=[69];powers=[]
    for exponent in (2,3,4):
        rows=[product(sparse(W[i,:]),b) for i in range(W.rows) for b in kb]
        W,_=rref(s.Matrix(rows));dims.append(W.rows);powers.append(data(W));print('ideal powers',dims,flush=True)
    need(dims[-1]==0 and dims[-2]>0,'exact radical nilpotence indexfour')
    # Pick a reproducible mixed nilpotent with maximal length if one exists.
    clock=None
    for seed in range(1,30):
        v=sum((((i+seed)*(i+seed+1)%7-3)*K[i,:] for i in range(69)),s.zeros(1,95));v2=product(sparse(v),sparse(v));v3=product(sparse(v2),sparse(v));v4=product(sparse(v3),sparse(v))
        if any(v3):
            need(not any(v4),'literal nilpotent clock orderfour');clock={'seed':seed,'coordinates':data(v)[0],
              'square':data(s.Matrix([v2]))[0],'cube':data(s.Matrix([v3]))[0]};break
    # The generators are actual algebra elements; the table must reproduce their known right actions.
    for name,raw in d['generators'].items():
        M=s.Matrix(raw);v=[M[p//n,p%n] for p in P];actual=s.Matrix([product({i:s.Integer(1)},sparse(v)) for i in range(95)])
        need(actual==s.Matrix(d['right_generators'][name]),'sparse presentation agrees with exact generator action')
    save('multiplication.json',{'status':'EXACT_COMPLETE_RATIONAL_MULTIPLICATION_TABLE','basis_dimension':95,
      'identity_coordinates':[int(p//n==p%n) for p in P],
      'nonzero_products':entries,'radical_power_dimensions':dims,'radical_power_bases':powers,'nilpotent_clock':clock})
    print('PASS sparse multiplication products',len(entries),'radical dimensions',dims,'clock seed',clock['seed'] if clock else None,flush=True)
if __name__=='__main__':main()
