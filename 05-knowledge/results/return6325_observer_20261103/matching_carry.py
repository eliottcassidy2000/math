"""Exact crossing matchings and clone-response debt on the actual return."""
import json
import sympy as s
from itertools import combinations
from support import HERE,source,need,save,junction

def matching(T,r):
    O=[i for i in range(len(T)) if T[r]>>i&1];I=[i for i in range(len(T)) if T[i]>>r&1];mate={}
    def aug(a,seen):
        for b in I:
            if not T[a]>>b&1 or b in seen:continue
            seen.add(b)
            if b not in mate or aug(mate[b],seen):mate[b]=a;return True
        return False
    for a in O:aug(a,set())
    edges=[(a,b) for a in O for b in I if T[a]>>b&1]
    cover=next(C for C in combinations(O+I,len(mate)) if all(a in C or b in C for a,b in edges))
    return sorted((a,b) for b,a in mate.items()),cover
def clone(T,r,m):
    names=[(i,j) for i in range(len(T)) for j in range(m if i==r else 1)]
    U=tuple(sum(((a==b and j<k) or (a!=b and (T[a]>>b&1)))<<v for v,(b,k) in enumerate(names)) for a,j in names)
    return U,names
def main():
    A,r,b,names=source();B=tuple(json.loads((HERE/'return6325.json').read_text())['final_rows']);out=[]
    for label,T in [('source',A),('returned',B)]:
        marks={}
        for v in (r,b):
            M,C=matching(T,v);marks[names[v]]={'matching':[[names[a],names[c]] for a,c in M],'cover':[names[a] for a in C]}
        rows=[]
        for m in range(1,6):
            U,ids=clone(T,r,m);z=junction.H(U)
            all_profiles=[junction.profile(U,ids.index((r,j))) for j in range(m)]
            profiles=[all_profiles[0],all_profiles[-1]]
            need(all(a[3]<c[3] for a,c in zip(all_profiles,all_profiles[1:])),'strict endpoint F along twin chain')
            rows.append({'m':m,'H':z,'first_last_profiles':profiles,'all_clone_ready_defects':[v[3]-z for v in all_profiles]})
            junction.endpoints.cache_clear()
        h=rows[0]['H'];c2=(rows[1]['H']-h)//2
        need(len(marks[names[r]]['matching'])==1 and all(d['H']==h+c2*(2**d['m']-2) for d in rows),'matching-one clone law')
        # m>=2: keep the old root as an exterior first/last vertex and clone
        # its adjacent twin. Its matching-one two-cover budget is three.
        coefficients=[]
        for end in (0,1):
            mat=s.Matrix([[j**t for j in (1,2,3)] for t in (1,2,3)])
            vals=s.Matrix([d['first_last_profiles'][end] for d in rows[1:4]])
            C=mat.inv()*vals
            U,ids=clone(T,r,2);v=ids.index((r,1 if end==0 else 0));M,_=matching(U,v)
            need(len(M)==1,'three-mode bound at the unmarked twin')
            for row in rows[1:]:
                need(s.Matrix([[j**(row['m']-1) for j in (1,2,3)]])*C==s.Matrix([row['first_last_profiles'][end]]),'independent extra clone sample')
            coefficients.append([[int(x) if x.q==1 else str(x) for x in C.row(i)] for i in range(3)])
        out.append({'state':label,'marked_boundaries':marks,'scalar_skeleton':[h,c2],'samples':rows,'first_last_profile_modes_m_minus_1':coefficients})
    diff=out[1]['scalar_skeleton'][1]-out[0]['scalar_skeleton'][1]
    save('matching_carry.json',{'status':'EXACT_COEFFICIENT_CERTIFICATE_WITH_INHERITED_MATCHING_BOUND','sources':out,
      'neutral_return_clone_c2_difference':diff,'scalar_return_difference':'delta_c2*(2^m-2)',
      'scope':'cloning the actual old root; matching bound proves the function-space upper bound; no uniform neutral return inferred from m=1'})
    print('PASS matching/carry root clone debt',diff,'*(2^m-2)',flush=True)
if __name__=='__main__':main()
