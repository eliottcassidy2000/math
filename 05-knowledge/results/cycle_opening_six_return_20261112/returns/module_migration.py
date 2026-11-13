"""Labeled, all-height source return with migration of a transitive module."""
import json
from itertools import combinations
import sparse_repair as s
Q=(14,12,8,64,175,79,535,367,127,447)
def target_on_source_names(k):
    expanded=s.inflate(Q,5,k)
    mapping=list(range(7))+[k+6,k+7,k+8]+list(range(7,k+6))
    out=[0]*len(expanded)
    for i,row in enumerate(expanded):
        for j in range(len(expanded)):
            if row>>j&1:out[mapping[i]]|=1<<mapping[j]
    return tuple(out)
def main():
    literal=[]
    for k in range(1,6):
        parent=s.inflate(s.P,1,k);source=s.c.birth(parent,4)
        raw,root=s.c.construct(source,len(source)-1,'LL');target=target_on_source_names(k)
        sv=s.c.state(source,len(source)-1);rv=s.c.state(raw,root);tv=s.c.state(target,root)
        D=2**k
        s.c.need(sv==(24*D+4,1,21*D+3),'ready source state')
        s.c.need(tv==(108*D-4,1,27*D+27),'ready migrated target state')
        s.c.need(rv[0]+rv[2]==tv[0]+tv[2]==3*(sv[0]+sv[2])+2,'literal full target preservation')
        diff=[(i,j) for i,j in combinations(range(len(raw)),2) if ((raw[i]^target[i])>>j)&1]
        fixed=[(0,6),(2,4),(3,4),(3,6),(6,k+8),(k+7,k+8)]
        expected=sorted(fixed+[(i,j) for j in range(7,k+6) for i in (0,1,6)])
        s.c.need(diff==expected,'explicit uniform migration arc recipe')
        s.c.need(len(diff)==3*k+3,'exact arc cost')
        s.c.need(s.c.strong(source) and s.c.strong(raw) and s.c.strong(target),'strong source and endpoints')
        literal.append(dict(k=k,source_state=sv,raw_state=rv,target_state=tv,retained_final_root=root,old_predecessor=k+7,new_predecessor=6,changed_arcs=diff))
    parentQ=tuple(row&511 for row in Q[:9])
    coeff=dict(H=s.exp_coeff(s.module_coeff(Q,5)),d=s.exp_coeff(s.module_coeff(parentQ,5)),b=s.exp_coeff(s.module_coeff(parentQ,5,start=6)))
    s.c.need(coeff==dict(H={1:23,2:135},d={1:27,2:27},b={1:1}),'exact all-height independent coefficient certificate')
    (s.HERE/'module_migration.json').write_text(json.dumps(dict(status='PROVED construction, VERIFIED independent exact compiler and literal k1..5',Q=Q,module=5,coefficients=coeff,literal=literal),indent=2)+'\n')
    print('PASS same-root full-target module migration for k1..5; all-height coefficient identities')
if __name__=='__main__':main()
