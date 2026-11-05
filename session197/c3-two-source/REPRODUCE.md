# Reproduction and inherited dependencies

The mathematical and executable snapshot is pinned by RELEASE_FILES.sha256. Run from a copy of this directory:

    python3 verify.py
    python3 -O verify.py

The verifier uses only Python's standard library. Both runs reproduce VALIDATION-normal.json and VALIDATION-optimized.json. It rewrites VALIDATION.json and SUPPORT.json locally. It requires no repository checkout or network access. The complete source bank is generated, not supplied as an oracle.

The proof copies in dependencies/ are inherited theorems, not newly claimed results. Exact immutable upstream copies include:

- [Complete C3 bridge theorem](https://github.com/eliottcassidy2000/math/blob/6ed2d7a0f1e1bfb26b9cc865ddc9f35dd0ad4567/session197/lowcycle/C3_COMPLETE.md).
- [Canonical start-three grammar](https://github.com/eliottcassidy2000/math/blob/6ed2d7a0f1e1bfb26b9cc865ddc9f35dd0ad4567/session197/lowcycle/root-amalgamation/START_THREE.md).
- [Ordered root amalgamation](https://github.com/eliottcassidy2000/math/blob/6ed2d7a0f1e1bfb26b9cc865ddc9f35dd0ad4567/session197/lowcycle/root-amalgamation/ROOT_AMALGAMATION.md).
- [Complete C2/source inverse](https://github.com/eliottcassidy2000/math/blob/f01df1766885c807aa4acae17685c2f4777f3ead/session197/bilinear-c2-inverse/PROOF.md).

The new result is a complete, possibly redundant, eight-template C3 source union and a constructive exact scalar-membership algorithm. The 501-bound numerical certificate is finite; it is not a proof of universal odd coverage. No output-root return or unique C3 presentation is claimed.

The independent review and its standalone C++ checker are in audit/. The checker uses only its three included fixture files; exact optimized/debug commands and fixture provenance are in audit/REPORT.md. The review binds the unchanged mathematical and Python snapshot.
