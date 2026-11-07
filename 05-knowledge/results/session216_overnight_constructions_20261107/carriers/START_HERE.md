# Reconstruct this size-bounded release

Every published transport file is at most131072 bytes. Oversized original banks, receipts and compressed files are omitted and represented by ordered base64 parts of at most120001 bytes. No mathematical content was removed.

Before running any research or audit command, run:

    python3 reconstruct.py

The helper works from any current directory, checks every part and reconstructed original against SHA256, refuses to overwrite a different existing file, and verifies RECONSTRUCTED_SHA256SUMS. It recreates originals at the exact relative paths expected by the scripts and nested audit manifests. Python3.9 or later is required.

For transport-only validation without creating oversized files:

    python3 reconstruct.py --check
    sha256sum -c SHA256SUMS

SHA256SUMS and PUBLICATION_FILES.txt describe the published size-bounded transport. ORIGINAL_RELEASE_SHA256SUMS and ORIGINAL_PUBLICATION_FILES.txt preserve prior release metadata when present; their historical names may differ from the new transport names. Nested frozen-audit manifests are unchanged and become directly checkable after reconstruction. RECONSTRUCTED_SHA256SUMS governs the full pre-split payload paths.

Do not publish the reconstructed oversized originals. Keep only the files listed in PUBLICATION_FILES.txt, plus SHA256SUMS and PUBLICATION_FILES.txt themselves. Reconstruction creates local working copies, not additional publication payload.
