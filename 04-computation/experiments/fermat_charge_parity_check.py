#!/usr/bin/env python3
"""Exact positional-charge fibre parity. Standard library; checks survive -O."""
import json
CHECKS=0

def require(value,label):
    global CHECKS
    CHECKS+=1
    if not value:raise ArithmeticError(label)

def parity_dp(p,S,ell,z):
    n=p-1;size=1<<n;limit=(1<<ell)-1
    require(pow(z,p,ell)==1 and z!=1 and (1+sum(pow(z,s,ell)for s in S))%ell==0,'quotient character')
    incoming=[sum(1<<u for u in range(n)if(v-u)%p in S)for v in range(n)]
    powers=[pow(z,v,ell)for v in range(p)]
    parity=[0]*(size*n);counts=[0]*(size*n)
    for mask in range(1,size):
        position=mask.bit_count();bits=mask
        while bits:
            b=bits&-bits;bits-=b;v=b.bit_length()-1;prior=mask^b;index=mask*n+v
            if not prior:
                if v+1 in S:parity[index]=1<<powers[v+1];counts[index]=1
            else:
                poly=0;number=0;ubs=prior&incoming[v]
                while ubs:
                    ub=ubs&-ubs;ubs-=ub;u=ub.bit_length()-1
                    poly^=parity[prior*n+u];number+=counts[prior*n+u]
                shift=position*powers[v+1]%ell
                parity[index]=((poly<<shift)&limit)|(poly>>(ell-shift))if shift else poly
                counts[index]=number
            require(parity[index].bit_count()%2==counts[index]%2,'each state agrees with ordinary path-count parity')
    out=0
    for poly in parity[(size-1)*n:]:out^=poly
    return out,sum(counts[(size-1)*n:])

def direct_counts(p,S,ell,z):
    powers=[pow(z,v,ell)for v in range(p)];counts=[0]*ell;nodes=0
    def visit(last,used,position,charge):
        nonlocal nodes
        nodes+=1
        if position==p:counts[charge]+=1;return
        for v in range(p):
            if not(used>>v&1)and(v-last)%p in S:
                visit(v,used|(1<<v),position+1,(charge+position*powers[v])%ell)
    visit(0,1,1,0)
    return counts,nodes

def main():
    rows=[]
    for p,S,ell,z in [(17,set(range(2,9))|{16},103,13),(11,{1,6,7,8,9},23,3)]:
        polynomial,total=parity_dp(p,S,ell,z)
        inverse_polynomial,inverse_total=parity_dp(p,{p-s for s in S},ell,pow(z,-1,ell))
        require(polynomial==inverse_polynomial and total==inverse_total,'independent relabelling control v -> -v')
        odd=[i for i in range(ell)if polynomial>>i&1];even=[i for i in range(ell)if not(polynomial>>i&1)]
        require(total%2==1 and len(odd)%2==1,'global Redei parity')
        support=set(odd);squares={a*a%ell for a in range(1,ell)}
        nonsquare=next(a for a in range(1,ell)if a not in squares)
        qr_matches=[(a,b)for a in [1,nonsquare]for b in range(ell)if{(a*x+b)%ell for x in squares}==support]
        affine_symmetries=[(a,b)for a in range(1,ell)for b in range(ell)if(a,b)!=(1,0)and{(a*x+b)%ell for x in support}==support]
        require(not qr_matches,'not an affine quadratic-residue support')
        require(not affine_symmetries,'trivial affine stabilizer of odd support')
        row={'p':p,'S':sorted(S),'quotient_prime':ell,'zeta_image':z,'normalized_path_count':total,'odd_fibre_count':len(odd),'even_fibre_count':len(even),'odd_residues':odd,'even_residues':even,'all_fibres_odd':len(odd)==ell,'affine_QR_representatives_tested':2*ell,'affine_QR_matches':qr_matches,'nonidentity_affine_maps_tested':ell*(ell-1)-1,'nonidentity_affine_symmetries':affine_symmetries}
        if p==17:
            require(total==793904749 and len(odd)==51 and len(even)==52 and 1 in even,'p17 exact parity counterexample')
            witness=[0,16,2,4,3,7,6,11,1,5,8,10,9,13,12,15,14]
            require(set(witness)==set(range(p))and all((v-u)%p in S for u,v in zip(witness,witness[1:]))and sum(i*pow(z,v,ell)for i,v in enumerate(witness))%ell==1,'explicit nonempty even fibre')
            row['charge_one_witness']=witness
        else:
            exact,nodes=direct_counts(p,S,ell,z)
            require(sum(exact)==total and sum((v%2)<<i for i,v in enumerate(exact))==polynomial,'independent exhaustive p11 histogram')
            require(all(exact)and len(odd)==11 and len(even)==12,'p11 positive coverage with even fibres')
            row.update({'direct_fibre_counts':exact,'direct_prefix_nodes':nodes})
        rows.append(row)
    # Exact arithmetic for the elementary all-ones obstruction to the Golay identity.
    require(11*12%3==0 and all(17*18%ell!=0 for ell in [103,239]),'all-ones quadratic obstruction controls')
    print(json.dumps({'status':'PASS','cases':rows,'active_checks':CHECKS},indent=2,sort_keys=True))
if __name__=='__main__':main()
