# Independent review: five-parameter positive-clock family

**PASS. No correction required.** Review date2026-11-05. Frozen source directory: `session197-hybrid-constructions/positive-clock-family/`.

- `PROOF.md`: `be3c701188a4a354cee62f21b6ee76bf172f67351028d88a872f4c95c9b865c3`
- `family.py`: `ab43dbe27e78bc11a176c5a2ba5f93c63b982c488612fd7a16570ca8b886bc45`
- `graphs.py`: `3204eeb1e9d2b05cac9bcc21c98f25697d10f925d6ff743393e945ce4e73615c`
- `verify.py`: `7de585230b22e2b030abfc509bad71fa17c1f29e5c945c0ed5d7c3cad78abf17`
- `VALIDATION.json`: `5e87d9fc75347afaa817dcaa3740afe96419f054d50b8d492190ed6d9de98ad4`

## Proof review

The seven-vertex quotient and its displayed Hamilton cycle are valid. Every crossing-cover vertex is one of the fixed singleton blocks a,b. The five covers therefore remain bounded when all five variable transitive blocks grow. The resulting separate run bounds3,2,2,2,2 imply the48-dimensional joint tensor span by successive invertible positive-point interpolation. The complete grid consequently proves the stated formula at every positive multiplicity.

Writing x=2^p,y=2^q,z=2^u,w=2^v gives exactly

    A=2zw(x-1)(y-1),
    B=w(y-1)(4-x)+z(x-1)(4-y),
    C=(x-3)(y-3).

The root deletion count is not the constant coefficient of this positive-clone expansion. In particular the source deletion count73 and constant1 are correctly distinguished.

The cancellation classification is complete. The same-side cases have definite sign. In the only remaining mixed case x=2,y>=8, cancellation would require z/w=2(y-1)/(y-4), strictly between2 and4, which is never an integral power of two, including negative powers. Thus B=0 exactly when p=q=2. Linear independence of1,2^m,3^m gives the asserted eventual-pure-ray conclusion. The restriction to this quotient family is essential and is stated correctly.

For p=q=2,u=1, the original root's outneighbors are precisely the four fixed ports1,2,3,4. The specified right birth appends a later transitive root clone and forgets only its temporary exposure. The displayed relabeling preserves every old port and every incidence. The returned scalar operation is3H-2 on the strong family72*2^(v-1)*3^m+1. Neither zero-length blocks nor arbitrary nontransitive substitutions are claimed.

## Independent controls

`positive-family-audit/independent.py` imports no producer code. It builds the quotient blow-ups directly and counts all48 grid graphs using signed walk inclusion–exclusion, independently of the source endpoint DP. Every physical count agrees with the formula. Literal outside-permutation/root-slot enumeration independently gives the full source vector `(73,217,216,72)` and all higher zeros.

The checker also verifies the five symbolic quotient covers,72 physical returning graphs and strongness through m=9,v=8, and20,736 cancellation cases through multiplicity12. The all-order cancellation proof above does not rely on that finite grid.

Copied-source normal and optimized runs agree with each other and the frozen receipt byte-for-byte. Original files were not modified.
