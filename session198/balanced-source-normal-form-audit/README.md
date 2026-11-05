# Independent balanced-source audit

The report is `balanced-replay/AUDIT.md`. Its frozen primary proof and verifier
hashes agree with the sibling `balanced-source-normal-form` packet.

Run both independent controls with normal Python (their audit assertions are
intended to remain enabled):

    python3 balanced-independent.py
    python3 balanced-repairs-independent.py

The second checker reads the sibling balanced packet and the already published
sibling `c1-monotone` packet. Only those two relative input paths were adapted
from the original scratch audit; all mathematics and checks are unchanged.
The adapted checker was replayed and reproduced the original audit output.
The main packet's own normal and optimized replays are separately certified.
