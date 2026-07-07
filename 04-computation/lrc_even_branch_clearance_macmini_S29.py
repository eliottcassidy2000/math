#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S29 -- toward a PROOF: the even-branch clearance lemma.

HYP-4572 trichotomy: for F(N) = {1..N}\{N-1} u {3(N-1)}, N even =>
M(F(N)) = 3/(3N-1) (binder speed 2 @ Q=3N-1).  Since 3/(3N-1) > 2/(2N+1) always
(3(2N+1)=6N+3 > 6N-2=2(3N-1)), proving the LOWER bound M(F(N)) >= 3/(3N-1)
suffices to show F(N) is NOT a gap member for even N.

M >= 3/(3N-1) needs ONE witness t=b/(3N-1) with ||v t|| >= 3/(3N-1) for ALL v in
F(N) -- i.e. all residues v*b mod (3N-1) avoid {0,+-1,+-2}.  This script:
 (1) find the closed form of b (binder 2*b == +-3 mod 3N-1);
 (2) print the FULL residue pattern for even N=6..20 -- look for a clean rule;
 (3) identify which v (if any) sit at exactly dist 3 (the binders) and confirm
     NO v sits in {0,+-1,+-2} -> the clearance holds -> proof template.
"""
from math import gcd
from fractions import Fraction as F


def family(N):
    base = [v for v in range(1, N + 1) if v != N - 1]
    return tuple(sorted(set(base) | {3 * (N - 1)}))


def even_branch_witness(N):
    """Q=3N-1 (odd for N even). binder speed 2: 2b == +-3 mod Q. Try both signs;
    return the b whose residue picture clears {0,+-1,+-2} for ALL of F(N)."""
    Q = 3 * N - 1
    W = family(N)
    inv2 = pow(2, -1, Q)
    for sign in (3, Q - 3):  # 2b == 3 or 2b == -3
        b = (sign * inv2) % Q
        res = [((v * b) % Q) for v in W]
        dists = [min(r, Q - r) for r in res]
        if all(d >= 3 for d in dists):
            return Q, b, W, res, dists
    return Q, None, W, None, None


def main():
    print("=" * 88)
    print("EVEN-BRANCH CLEARANCE: witness t=b/(3N-1), all residues avoid {0,+-1,+-2}?")
    print("=" * 88)
    print(f"  {'N':>3} {'Q=3N-1':>7} {'b':>4} {'b closed form?':>22} {'min dist':>9} {'binders (v@dist3)':>20}")
    ok_all = True
    binder_data = []
    for N in range(6, 22, 2):  # even N
        Q, b, W, res, dists = even_branch_witness(N)
        if b is None:
            print(f"  {N:>3} {Q:>7}  -- NO clearing witness at Q=3N-1 (unexpected)")
            ok_all = False
            continue
        mind = min(dists)
        binders = [W[i] for i in range(len(W)) if dists[i] == 3]
        # closed-form guesses for b:  3*inv2 mod Q with inv2=(Q+1)/2=(3N)/2
        inv2 = (Q + 1) // 2
        b_plus = (3 * inv2) % Q
        b_minus = ((Q - 3) * inv2) % Q
        cf = ("b=3*inv2=(9N/2)%Q" if b == b_plus else
              ("b=-3*inv2%Q" if b == b_minus else "??"))
        print(f"  {N:>3} {Q:>7} {b:>4} {cf:>22} {str(F(mind,Q)):>9} {str(binders):>20}")
        binder_data.append((N, binders, mind))
        if mind != 3:
            ok_all = False
    print()
    print("  => for every even N, a witness at Q=3N-1 clears {0,+-1,+-2}, min dist = 3")
    print(f"     => M(F(N)) >= 3/(3N-1) > 2/(2N+1)  => F(N) NOT a gap member.  holds: {ok_all}")
    print(f"     binders (the speeds at dist exactly 3): {[(N,b) for N,b,_ in binder_data]}")

    print()
    print("=" * 88)
    print("FULL RESIDUE PATTERN at even N=12 (our LRC-14 case) -- the proof witness")
    print("=" * 88)
    Q, b, W, res, dists = even_branch_witness(12)
    print(f"  F(12) = {list(W)}")
    print(f"  witness t = {b}/{Q}   (Q=35=5*7, b={b})")
    print(f"  {'v':>4} {'v*b mod 35':>11} {'dist to 0':>10} {'>=3?':>5}")
    for i, v in enumerate(W):
        print(f"  {v:>4} {res[i]:>11} {dists[i]:>10} {'OK' if dists[i]>=3 else 'BAD':>5}")
    print(f"  min dist = {min(dists)} => min ||v t|| = {F(min(dists),Q)} = 3/35")
    print(f"  3/35 = {float(F(3,35)):.5f} > 2/25 = {float(F(2,25)):.5f}  => M(F(12)) >= 3/35 > 2/25")
    print(f"  => F(12) is NOT a gap member (its M is ABOVE the gap).  [Lean-able finite check]")

    # Verify the closed form b = 3*inv2 mod Q holds uniformly, and residue formula
    print()
    print("=" * 88)
    print("CLOSED-FORM RESIDUE RULE (for the proof): v*b mod Q with b=3*inv2")
    print("=" * 88)
    print("  b = 3*inv2 = 3*(3N/2) = 9N/2 mod (3N-1).  Since 9N/2 = (3(3N-1)+3)/2,")
    print("  and residues v*b = v*3*inv2 = (3v)*inv2 = 3v/2 mod Q.")
    print("  So residue(v) = 3v * inv2 mod Q; forbidden <=> 3v == {0,+-1,+-2}*2 = {0,+-2,+-4} mod Q.")
    print("  i.e. 3v in {0,2,4,Q-2,Q-4} mod Q.  Check which v in F(N) hit these:")
    for N in [12]:
        Q = 3 * N - 1
        W = family(N)
        forbidden_3v = {0, 2, 4, Q - 2, Q - 4}
        for v in W:
            hit = (3 * v) % Q in forbidden_3v
            tag = "  <-- would collide" if hit else ""
            if hit or v in (2, 3*(N-1)):
                print(f"    v={v}: 3v mod {Q} = {(3*v)%Q}  {'FORBIDDEN' if hit else 'ok'}{tag}")
        print(f"    (binders: 3v==+-4 => v with 3v=4 or Q-4; e.g. speed 2: 3*2=6, hmm)")


if __name__ == "__main__":
    main()
