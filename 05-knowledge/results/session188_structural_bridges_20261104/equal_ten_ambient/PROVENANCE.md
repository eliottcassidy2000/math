# Provenance of the finite Hamilton inputs

Checked on4 November2026. The official [OEIS entryA038375](https://oeis.org/A038375), maintained as the sequence of maximal spanning-path counts in tournaments, records5 at order4 and45 at order6. It cites J.W.Moon, Topics on Tournaments(1968), page28, and Reid--Beineke's Tournaments chapter(1978). We inspected the sequence entry, not those books, and make no unsupported claim about the first proof or exact location of each value in them.

The same record gives the sharper value15,745 at order10. That larger-order enumeration is not needed or imported into our proof. We instead independently replay both complete small labelled universes and use the elementary6+4 split bound47,250, whose weaker value already suffices for the response certificate.

No new extremal Hamilton count is claimed. The exact same-parity profile cone is classical Grinberg--Stanley; the certificate and ambient lifting dependencies are separately cited in the proof. The historical context and stronger asymptotic theory are recorded in `../equal_nine_ambient/PROVENANCE.md`.
