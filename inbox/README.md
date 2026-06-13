# Inbox

Drop files here to contribute to the research system. This is the intake point for:
- Chat exports from Claude sessions on the web interface
- LaTeX, Markdown, HTML, or plain text documents from any Claude instance
- Code scripts used for verification
- Any other relevant mathematical work

---

## Workflow

1. **Drop file(s) into this directory** (inbox/)
2. **Run the processor:**
   ```bash
   python3 /path/to/math-research/inbox/processor.py
   ```
   This generates `inbox/PROCESSING-REPORT.md`
3. **Tell Claude:** "Process the inbox" or "There are new files in the inbox"
4. **Claude reads** `PROCESSING-REPORT.md` and integrates the content
5. Processed files are moved to `inbox/processed/YYYY-MM-DD/`

---

## What the Processor Does

- **Exact duplicate:** Detected by SHA-256 hash. No action needed; file is archived.
- **Near-duplicate:** Detected by text similarity > 85%. A diff is shown for Claude to review. Discrepancies between similar files are flagged for potential court cases.
- **New file:** Content preview is shown. Claude integrates it into the appropriate parts of the system (canon theorems, MISTAKES, TANGENTS, etc.).

---

## File Format Notes

- **Preferred input:** Markdown (.md)
- **Also accepted:** Plain text (.txt), LaTeX (.tex — archived to 03-artifacts/drafts/, content extracted to Markdown), Python scripts (.py — archived to 03-artifacts/code/, indexed)
- **HTML/DOCX:** Accepted but will be converted to Markdown during processing

---

## For Other Claude Instances Submitting Work

Until a more polished inter-Claude submission system is set up, the recommended approach is:
1. Export your session's key findings as a Markdown file
2. Name it descriptively: `[INSTANCE-ID]-[DATE]-[topic].md`
3. The human will drop it in this inbox
4. The next Claude session will integrate it

Use the instance naming convention: `[ACCOUNT]-[DATE]-S[N]`
