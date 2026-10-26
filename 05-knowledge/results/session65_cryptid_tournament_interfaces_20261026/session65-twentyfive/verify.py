#!/usr/bin/env python3
"""Exact source parsing, control relabeling, finite trace and arithmetic checks.
No Lean build, infinite-run conclusion, or conjecture proof is claimed.
"""
from pathlib import Path
from collections import deque
import hashlib,json,re
ROOT=Path(__file__).resolve().parent
LEAN=ROOT/"inputs/goldbach_tm"
WIKI="1RB1RD_1LC1RB_0RA1LC_0LQ1RE_0LF1RG_0LC1LF_0LF0LH_1LQ1LI_0RJ0LI_1RK0LJ_0RL0RS_1RL0RM_1RN1RM_0LO0LU_0LP1LO_1RH1LX_1LR1LQ_0RK0LT_1LR1RS_---1RC_1LV1LU_0LW0LJ_0RK0LW_1RY1LX_1RE1RY"
def require(x,msg):
    if not x: raise RuntimeError(msg)
def parse_lean():
    s=(LEAN/"GoldbachTm/Tm25/TuringMachine25.lean").read_text()
    pat=r"\| ⟨\s*(\d+), _⟩, Γ\.(zero|one)\s*=> some ⟨⟨\s*(\d+), by omega⟩, ⟨Turing.Dir\.(right|left), Γ\.(zero|one)⟩⟩"
    ans={}
    for q,r,t,d,w in re.findall(pat,s):
        ans[int(q),int(r=="one")]=(int(w=="one"),1 if d=="right" else -1,int(t))
    for q,r in re.findall(r"\| ⟨\s*(\d+), _⟩, Γ\.(zero|one)\s*=> none",s):
        ans[int(q),int(r=="one")]=None
    require(len(ans)==50,"all 50 branches must parse")
    return ans
def parse_wiki():
    rows=WIKI.split("_")
    require(len(rows)==25 and all(len(r)==6 for r in rows),"25 binary states")
    return {(q,b):None if row[3*b:3*b+3]=="---" else
            (int(row[3*b]),1 if row[3*b+1]=="R" else -1,ord(row[3*b+2])-65)
            for q,row in enumerate(rows) for b in (0,1)}
def conjugacy(a,b):
    rho={0:0}; todo=deque([0]); checks=0
    while todo:
        q=todo.popleft()
        for r in (0,1):
            aa,bb=a[q,r],b[rho[q],r];checks+=1
            if aa is None or bb is None:
                require(aa==bb,"halting preserved");continue
            require(aa[:2]==bb[:2],"write/move preserved")
            x,y=aa[2],bb[2]
            if x in rho: require(rho[x]==y,"consistent state map")
            else: rho[x]=y;todo.append(x)
    require(len(rho)==25 and len(set(rho.values()))==25,"bijection")
    return rho,checks
def trace(a,rho,maxsteps=2000000):
    ones=set();q=0;head=0;macro=[];seen=set();cfg45=False
    for step in range(maxsteps+1):
        seen.add(q)
        if step==45:
            cfg45=(q==21 and head not in ones and ones==set(range(head-4,head)))
        if q==21 and head not in ones and ones and ones==set(range(head-len(ones),head)):
            n=len(ones)
            if not macro or macro[-1]["n"]!=n:macro.append({"n":n,"step":step})
        if step==maxsteps:break
        t=a[q,int(head in ones)]
        require(t is not None,"finite trace unexpectedly halts")
        w,d,q=t
        if w:ones.add(head)
        else:ones.discard(head)
        head+=d
    require(cfg45,"formal cfg45 matches independent execution")
    require([x["n"] for x in macro]==list(range(4,macro[-1]["n"]+1,2)),"consecutive even macro blocks")
    return {"steps":maxsteps,"halted":False,"visited_control_states":len(seen),
            "cfg45":cfg45,"macro_blocks":macro}
def prime(n):return n>=2 and all(n%d for d in range(2,__import__("math").isqrt(n)+1))
def arithmetic(nmax=500):
    ps=[p for p in range(2,nmax+1) if prime(p)];coeff=[0]*(2*nmax+1)
    for p in ps:
        for q in ps:coeff[p+q]+=1
    checks=0; samples={}
    for n in range(4,nmax+1,2):
        pairs=[(p,n-p) for p in ps if p<=n and prime(n-p)]
        require(coeff[n]==len(pairs),"convolution witness identity")
        require(coeff[n]>0,"finite Goldbach control")
        if n<=20:samples[n]=pairs
        checks+=1
    # Pair marginals lose compatibility: they cannot certify a shared witness.
    A={(0,0),(1,1)}; B={(0,1),(1,0)}
    require({x for x,y in A}=={x for x,y in B},"first marginals same")
    require({y for x,y in A}=={y for x,y in B},"second marginals same")
    require(not(A&B),"shared-pair witness absent")
    # Most elementary state-only hostile: q24 halts on 0 but not 1.
    return {"even_numbers_tested":checks,"maximum_even":nmax,
            "convolution_counts_match":True,"small_witnesses":samples,
            "same_marginals_disjoint_pair_fibres":{"A":sorted(A),"B":sorted(B)}}
def padding(a):
    # Replace the initially visited (state0, symbol0) transition.
    w,d,target=a[0,0]
    require(d==1,"chosen transition moves right")
    b=dict(a);b[0,0]=(w,d,25)
    for bit in (0,1):
        b[25,bit]=(bit,1,26)
        b[26,bit]=(bit,-1,target)
    def run(M,ones,n):
        ones=set(ones);q=0;h=0
        for _ in range(n):
            v,move,q=M[q,int(h in ones)]
            if v:ones.add(h)
            else:ones.discard(h)
            h+=move
        return q,h,ones
    for bits in range(4):
        ones={i+1 for i in (0,1) if bits>>i&1}
        require(run(a,ones,1)==run(b,ones,3),"padding preserves full macroconfiguration")
    return {"original_states":25,"padded_states":27,"new_states_reached_from_blank":True,
            "local_tape_patterns_checked":4,"extra_steps_per_replaced_transition":2}
def main():
    a,b=parse_lean(),parse_wiki();rho,count=conjugacy(a,b)
    static={}
    for p in sorted(LEAN.rglob("*.lean")):
        s=p.read_text()
        static[str(p.relative_to(LEAN))]=re.findall(r"\b(?:sorry|admit|axiom|native_decide|sorryAx)\b",s)
    require(not any(static.values()),"no flagged proof placeholders in saved source")
    ps=[p for p in range(3,201,2) if prime(p)]
    for p in ps:
        for q in ps:
            E,O=p+q,p+2*q
            require((2*E-O,O-E)==(p,q),"ordered-pair inverse")
            require(q+2*p==3*E-O,"swap reflection")
            require((2*O==3*E)==(p==q),"diagonal fixed locus")
    V=[26*i for i in range(1,13)]+[339]
    from math import gcd
    periods={q:[p for p in range(q) if gcd(p,q)==1 and
                 all(14*min((p*v)%q,q-(p*v)%q)>=q for v in V)]
             for q in range(1,28)}
    require(all(not periods[q] for q in range(1,27)),"no period through 26")
    require(periods[27]==[2,25],"first period 27")
    require(a[24,0] is None and a[24,1] is not None,"state alone loses read-symbol semantics")
    out={"status":"FINITE-EXACT; static source audit, not a Lean rebuild",
         "external_commit":"6cd292062516a0a14ea1b34f2ab75154cae7ab1e",
         "control_isomorphism":{"lean_state_to_wiki_state":{str(k):chr(65+v) for k,v in sorted(rho.items())},
                                "branches_checked":count,"defined_instructions":49,"halting_branch":[24,0]},
         "trace":trace(a,rho),"arithmetic":arithmetic(),
         "reachable_padding":padding(a),
         "pair_shadow_checks":{"ordered_prime_pairs_through_200":len(ps)**2,
                               "inverse_swap_and_diagonal":True},
         "lrc_cutoff_hostile":{"speeds":V,"first_good_reduced_denominator":27,
                               "good_numerators_at_27":[2,25]},
         "static_flag_counts":{p:len(v) for p,v in static.items()}}
    print(json.dumps(out,indent=2,sort_keys=True))
if __name__=="__main__":main()
