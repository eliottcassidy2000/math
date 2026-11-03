"""Exact finite-operator proof of the golden phase formulas after every sealing length."""
from fractions import Fraction
import json,sympy as s
from observer import HERE,old,observer,need
from golden_transport import graded,phase

STATES=tuple((sk,m) for sk in graded.skeletons(2) for m in range(3) if m+sum(b==-1 for a,b in sk)<=2)
INDEX={state:i for i,state in enumerate(STATES)}

def matrix(command):
    return tuple(tuple((INDEX[target],c) for target,c in graded.kernel(sk,m,2,*command).items() if target in INDEX) for sk,m in STATES)

def back(q,M):return tuple(sum(c*q[j] for j,c in row) for row in M)

def output(slot,kind):
    result=[]
    for sk,m in STATES:
        k=len(sk)+m
        h=int(k==1);e=int(k==1 and any(b==slot for a,b in sk));S=int(k in (1,2) and any(a==slot for a,b in sk))
        result.append(h if kind==0 else Fraction(S+(-1 if kind==1 else 1)*e,2))
    return tuple(result)

def suffix(root,cycle):
    seq=[(0,1<<root,1-root),(0,1<<(1-root),root)];root=1-root # N R
    if cycle:
        seq.append((0,1<<root,1-root));root=1-root
        seq.append((0,1<<root,1-root))
        seq.extend([(0,1<<root,2)]*4)
        root=1-root;seq.append((0,1<<root,1-root));root=1-root
        seq.append((0,1<<root,1-root))
    return seq,root

def main():
    w,u=old.chart(3,4);V=observer(old.source(w,u,False),(7,8));initial=tuple(V.get(key,0) for key in STATES)
    formulas=[];certificates=[];rows=json.loads((HERE/'golden_transport.json').read_text())['examples']
    for root in (0,1):
        R=matrix((0,1<<root,1-root))
        for cycle in (False,True):
            commands,slot=suffix(root,cycle)
            key='C5_golden' if cycle else 'sealed_golden';coeffs=[]
            for kind in range(3):
                q=output(slot,kind)
                for cmd in reversed(commands):q=back(q,matrix(cmd))
                powers=[q]
                for _ in range(4):powers.append(back(powers[-1],R))
                # R(R-I)(R-2I)(R-3I) annihilates the actual observation row.
                residual=tuple(a-6*b+11*c-6*d for a,b,c,d in zip(powers[4],powers[3],powers[2],powers[1]))
                need(not any(residual),'universal golden observation annihilator for k>=1')
                values=[sum(a*b for a,b in zip(initial,powers[k])) for k in (1,2,3)]
                C=s.Matrix([[3**k,2**k,1] for k in (1,2,3)]).inv()*s.Matrix(values)
                coeffs.append([int(c) if c.q==1 else str(c) for c in C])
                for row in rows:
                    if row['k']%2==(0 if root==0 else 1):
                        k=row['k'];need(sum(C[i]*base**k for i,base in enumerate((3,2,1)))==row[key][kind],'formula versus full polynomial transport')
            formulas.append({'root':7+root,'stage':key,'bases':[3,2,1],'coefficients':coeffs})
            certificates.append({'root_slot':root,'stage':key,'annihilator':'x*(x-1)*(x-2)*(x-3)','row_identities':3})
    # All-order field clock follows from the formula, not a repeated prefix alone.
    phases=[]
    for k in range(1,7):
        F=next(f for f in formulas if f['stage']=='C5_golden' and f['root']==(7 if k%2==0 else 8))
        v=tuple(sum(int(c)*base**k for c,base in zip(row,(3,2,1))) for row in F['coefficients'])
        phases.append(phase(v))
    need(len(set(tuple(p) for p in phases))==6,'minimal full golden phase period six')
    # C5 multiples of determinant have no claimed relation to a Q4 clock.
    need(phases[1]!=[1,0,4],'matching six-clock positions do not identify different sourced charts')
    out={'status':'PROVED_ALL_K_GE_1','finite_cap_two_end_states':len(STATES),'state_keys':STATES,
         'annihilator_certificates':certificates,'formulas':formulas,'C5_golden_period':6,'C5_golden_phases':phases,
         'selected_scalar_recurrence':'u[k+6]=14*u[k+4]-49*u[k+2]+36*u[k]',
         'ready_family_phase_at_k_2_mod6':[1,0,4]}
    (HERE/'phase_recurrence.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS: universal capped-operator annihilators; all-order golden formulas and exact six-clock')
if __name__=='__main__':main()
