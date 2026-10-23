# Session38: incidence, ports and electrical constraints

Start with the [audited synthesis](SYNTHESIS.md). This packet combines complete order-eight complement-line wiring, finite incidence-signature determination, and all-order proofs with explicit limits.

- [Complete wiring](wiring/README.md) and [independent audit](wiring-audit/AUDIT.md): every one of 2,097,152 fixed-path presentations; loops, black twins and parallel multiplicities retained.
- [Incidence signatures](incidence/README.md), [converse repair](incidence/CONVERSE_SIGNATURES.md), [ternary compression](incidence/TERNARY_SIGNATURES.md), and [aligned decoding](incidence/THREE_COLOR_DECODING.md).
- [Electrical laws](electrical/PROOF.md) and [bridge localization](electrical/BRIDGE_LOCALIZATION.md): the exact two-primary obstruction and its bridge-free core.
- [Binary counterfeits](binary-counterfeits/README.md): every two-colour quotient in the declared nine-symbol family has a certified alternative literal matching.
- [Transport rules](rule-audit/AUDIT.md) and [reflected counterfeit audit](rule-audit/REFLECTION_AUDIT.md).
- [Transitive-profile obstruction](theory/TRANSITIVE_PROFILE_OBSTRUCTION.md) with [independent audit](transitive-audit/AUDIT.md), and [growing-alphabet boundary](theory/GROWING_ALPHABET_BOUNDARY.md).

## Reproduction

Use each packet's instructions. Producer wiring and independent C++ enumeration generate their own complete ledgers. The raw generated class map and audit intermediates are not publication inputs. Deterministic gzip files are verified byte-for-byte.

From this packet root, incidence/check.py takes --classes wiring/classes_input_n8.tsv --ports wiring/ports.tsv.gz --capacities wiring/endpoint_capacities.tsv.gz. Run its structural checker separately. Electrical scripts accept --base wiring (check_electrical.py instead uses --n8 wiring); run check_electrical.py without --n8 for the small controls. For transitive-audit/check.py use --producer theory. Run theory/check_residual_counterfeits.py with --incidence incidence --wiring wiring and binary-counterfeits/check.py with --wiring wiring. Independent wiring-audit scripts consume generated audit intermediates and sibling packet data as described in wiring-audit/README.md.

Some frozen electrical receipts retain historical absolute source-path keys as provenance. Compare those keys by basename after relocation while retaining every hash value and every mathematical field; do not rewrite the stored receipt chain. Producer wiring stdout contains elapsed seconds, which are normalized only for transcript comparison. Mathematical receipts contain no timing tolerance. Run normal Python for assertions where present; optimized-mode agreement alone does not establish that removed assertions executed.

All results retain their stated finite universes, signature menus and graph conventions. No general black-core connectivity theorem, unrestricted minimal encoding theorem, complete odd critical-group computation or LRC(14) proof is claimed. Private attached-paper critiques and previously withheld material remain excluded.
