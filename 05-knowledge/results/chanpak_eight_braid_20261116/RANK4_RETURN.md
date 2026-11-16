# 1024 dyadic targets while returning the actual four-plane

The [five-letter all-strong bank](RANK4_STEERING.md) preserves rank but does not generally return to the same subspace. Here is a different, finite exact result that does return.

Fix the strong seed (0,0,3,5,2,4), H=25. In the literal after-first model let P be the action of its last five letters S=(0,3,5,2,4), and E=image(P), of rank four. Use a common refresh suffix R=(3,0,3,5,2,4) that ends in S.

**Exact return theorem.** For each of the **1024 odd residues modulo 2048**, some word u of length at most five, with only nonzero letters, followed by R has both properties:

1. Its count when appended to the seed is the specified residue.
2. Its induced operator maps E isomorphically **onto E**, not just to another four-dimensional image.

Consequently every finite or infinite concatenation of these return blocks preserves rank four throughout and returns to exactly the same image at every block boundary. The zero in R may temporarily break strongness, so this is separate from the all-strong modulo-1024 result.

**Proof and replay.** The image of the full refresh is inside E because its last five edges are S; the verifier directly checks the same rational inclusion in semantic coordinates. Give E a basis by four selected columns of the seed action on the reachable space. For each residue the sidecar lists u in increasing residue order. Propagate the four columns through u and R: they stay independent after projection modulo 10007. Thus their rational rank is at least four, at most four by inclusion, and the map is an isomorphism on E. Multiplication cannot increase rank, so all intermediate ranks are four too. In parallel the exact Hamiltonian path transfer modulo 2048 gives the specified physical readout. `check_rank4_return.py` verifies `RETURN_MOD2048.json` with the same rational quotient and integer modular operations used in the preceding certificate.

**What repeats and what does not.** The return is independent of the current vector in E. The stated menu of count residues was certified for the **seed vector**, not uniformly on all future vectors. In fact even the smaller checked mod-32 fixed menu loses options after a single turn (see [Three moving strands](THREE_LEASHES.md)). This is the precise distinction between bringing back the subspace, bringing back a state up to a unit, and achieving arbitrary online choices of the scalar while time keeps moving. The result is not a density theorem or an absolute target selector. It supplies a concrete family of invertible four-dimensional return actions on which that stronger state-control question can be asked.
