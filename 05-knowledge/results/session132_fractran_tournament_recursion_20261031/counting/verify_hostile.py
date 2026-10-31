#!/usr/bin/env python3
"""Negative controls for plausible but false enumeration shortcuts."""
from fractions import Fraction
from math import factorial
import json
from check_counts import (Census, fixed_count_from_signed_pair_cycles, pair_action,
                          partition_count, require, verify_extension)


def main():
    cases = []
    # A vertex transposition reverses its unique edge: no ordinary invariant.
    a = pair_action(2,(1,0))
    require(fixed_count_from_signed_pair_cycles(a) == 0,
            'ordinary even-cycle exclusion')
    require(fixed_count_from_signed_pair_cycles(a,True) == 2,
            'twist makes the transposition consistent')
    cases.append('ordinary cycle exclusion differs from converse exclusion')
    # An odd 3-cycle fixes two cyclic tournaments, but no twisted tournament.
    a = pair_action(3,(1,2,0))
    require(fixed_count_from_signed_pair_cycles(a) == 2 and
            fixed_count_from_signed_pair_cycles(a,True) == 0,
            'ordinary odd cycles are not converse-admissible')
    cases.append('odd nontrivial cycles rejected in the twisted sum')
    # Even is not enough: cycles divisible by four fail the antipodal gate.
    require(fixed_count_from_signed_pair_cycles(pair_action(4,(1,2,3,0)),True) == 0,
            'divisible-by-four anti cycle rejected')
    cases.append('anti cycle lengths need 2 modulo 4, not merely even')
    # Restricting to anti-involutions in a Burnside average is incorrect.
    involution_only = Fraction(15*2**9,factorial(6))
    require(involution_only == Fraction(32,3) and involution_only != 12,
            'all anti-isomorphisms, not just anti-involutions, must be summed')
    cases.append('anti-involution-only average at n=6 is 32/3, not 12')
    # Formal F(q,x) need not have integer coefficients for arbitrary q.
    try:
        partition_count(2,1,1)
    except RuntimeError as e:
        require('integral' in str(e),'expected integrality failure')
    else:
        raise RuntimeError('q=1 must not silently round 1/2')
    cases.append('non-combinatorial q=1 rejected rather than rounded')
    c1,c2,c3 = Census(1),Census(2),Census(3)
    r = verify_extension(c2,c3,False)
    require(r['extension_orbits'] == 4 and r['naive_divide_by_nplus1'] == '4/3'
            and r['weighted_unpointed'] == 2,'root-count division hostile')
    cases.append('four root classes divided by three do not give two child types')
    plain = verify_extension(c1,c2,False)
    signed = verify_extension(c1,c2,True)
    require(plain['extension_orbits'] == 2 and signed['extension_orbits'] == 1,
            'attachment complements must be included in signed action')
    cases.append('dropping the converse sign doubles the first merged root count')
    print(json.dumps({'status':'PASS','rejected_shortcuts':cases},indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
