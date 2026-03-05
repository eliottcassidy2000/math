#!/usr/bin/env python3
"""
Multi-Agent Message Processor
==============================
Handles inter-Claude messaging for the N-machine research network.

Usage:
    python3 agents/processor.py --check          # Read incoming messages
    python3 agents/processor.py --send           # Write end-of-session letter (interactive)
    python3 agents/processor.py --send           # (also accepts piped/file input, see below)
    python3 agents/processor.py --register       # Register this machine as a new agent
    python3 agents/processor.py --list           # List all known agents
    python3 agents/processor.py --status         # Full inbox status across all agents

Sending non-interactively (for scripting):
    python3 agents/processor.py --send --to bob --subject "Claim A update" --body-file letter.md
    python3 agents/processor.py --send --to all --subject "Broadcast" --body "Short message"
"""

import os
import sys
import json
import shutil
import argparse
import textwrap
from datetime import datetime
from pathlib import Path

# Windows: force UTF-8 so box-drawing chars and emoji don't crash on CP1252
if sys.stdout.encoding and sys.stdout.encoding.lower() not in ('utf-8', 'utf8'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
if sys.stderr.encoding and sys.stderr.encoding.lower() not in ('utf-8', 'utf8'):
    sys.stderr.reconfigure(encoding='utf-8', errors='replace')

# ─── Paths ────────────────────────────────────────────────────────────────────
SCRIPT_DIR  = Path(__file__).resolve().parent
REPO_ROOT   = SCRIPT_DIR.parent
AGENTS_DIR  = SCRIPT_DIR
MACHINE_ID_FILE = REPO_ROOT / ".machine-id"
REGISTRY    = AGENTS_DIR / "REGISTRY.md"
BROADCAST   = AGENTS_DIR / "broadcast"
READ_LOG    = AGENTS_DIR / ".read-log.json"    # tracks which messages we've read

# ─── Identity ─────────────────────────────────────────────────────────────────

def get_machine_id() -> str:
    if not MACHINE_ID_FILE.exists():
        print("ERROR: .machine-id not found.")
        print("Run: echo 'your-machine-name' > .machine-id")
        print("Then: python3 agents/processor.py --register")
        sys.exit(1)
    return MACHINE_ID_FILE.read_text().strip()


def get_my_inbox() -> Path:
    machine_id = get_machine_id()
    return AGENTS_DIR / machine_id / "inbox"


def get_all_agents() -> list[str]:
    """Return list of registered agent names (directories in agents/ that have identity.md)."""
    agents = []
    for d in AGENTS_DIR.iterdir():
        if d.is_dir() and not d.name.startswith("_") and d.name != "broadcast":
            if (d / "identity.md").exists():
                agents.append(d.name)
    return sorted(agents)


def load_read_log() -> set[str]:
    if READ_LOG.exists():
        try:
            return set(json.loads(READ_LOG.read_text()))
        except Exception:
            pass
    return set()


def save_read_log(read: set[str]):
    READ_LOG.write_text(json.dumps(sorted(read), indent=2))


# ─── Message helpers ──────────────────────────────────────────────────────────

def next_msg_number(inbox_dir: Path) -> str:
    """Find next available MSG number in a given inbox directory."""
    existing = list(inbox_dir.glob("MSG-*.md"))
    if not existing:
        return "001"
    nums = []
    for f in existing:
        try:
            nums.append(int(f.name.split("-")[1]))
        except (IndexError, ValueError):
            pass
    return str(max(nums) + 1).zfill(3) if nums else "001"


def msg_filename(inbox_dir: Path, sender: str, subject: str) -> str:
    date = datetime.now().strftime("%Y-%m-%d")
    num  = next_msg_number(inbox_dir)
    slug = subject.lower().replace(" ", "-")[:30].rstrip("-")
    slug = "".join(c for c in slug if c.isalnum() or c == "-")
    return f"MSG-{num}-from-{sender}-{date}-{slug}.md"


def get_unread_messages(inbox_dir: Path, read_log: set[str]) -> list[Path]:
    if not inbox_dir.exists():
        return []
    msgs = sorted(inbox_dir.glob("MSG-*.md"))
    return [m for m in msgs if str(m.relative_to(REPO_ROOT)) not in read_log]


# ─── Commands ─────────────────────────────────────────────────────────────────

def cmd_check(args):
    """Show all unread messages for this machine."""
    machine_id = get_machine_id()
    read_log   = load_read_log()
    my_inbox   = get_my_inbox()

    print(f"\n{'='*60}")
    print(f"  Incoming messages for: {machine_id}")
    print(f"{'='*60}\n")

    # Direct messages
    direct = get_unread_messages(my_inbox, read_log)
    # Broadcast messages
    bc     = get_unread_messages(BROADCAST, read_log)

    total = len(direct) + len(bc)
    if total == 0:
        print("No unread messages. You're up to date.\n")
        return

    newly_read = set()

    if direct:
        print(f"── Direct messages ({len(direct)}) ──────────────────────────\n")
        for msg in direct:
            _print_message(msg)
            newly_read.add(str(msg.relative_to(REPO_ROOT)))

    if bc:
        print(f"── Broadcast messages ({len(bc)}) ─────────────────────────\n")
        for msg in bc:
            _print_message(msg)
            newly_read.add(str(msg.relative_to(REPO_ROOT)))

    save_read_log(read_log | newly_read)
    print(f"\n✓ Marked {len(newly_read)} message(s) as read.\n")


def _print_message(path: Path):
    print(f"  📨  {path.name}")
    print(f"  {'─'*56}")
    try:
        content = path.read_text(encoding="utf-8")
        # Print full message with slight indent
        for line in content.splitlines():
            print("  " + line)
    except Exception as e:
        print(f"  [Error reading file: {e}]")
    print()


def cmd_send(args):
    """Write and deliver a message to one or more agents."""
    machine_id = get_machine_id()
    all_agents = get_all_agents()

    # ── Determine recipients ──
    if args.to:
        to = args.to.strip().lower()
    else:
        print("\nSend to:")
        print("  all         — broadcast to everyone")
        for a in all_agents:
            if a != machine_id:
                print(f"  {a}")
        to = input("\nRecipient (name or 'all'): ").strip().lower()

    if to == machine_id:
        print("ERROR: You cannot send a message to yourself.")
        sys.exit(1)

    if to != "all" and to not in all_agents:
        print(f"ERROR: Unknown agent '{to}'. Register them first or check the spelling.")
        print(f"Known agents: {', '.join(all_agents)}")
        sys.exit(1)

    # ── Subject ──
    if args.subject:
        subject = args.subject
    else:
        subject = input("Subject: ").strip()

    # ── Body ──
    if args.body:
        body = args.body
    elif args.body_file:
        body = Path(args.body_file).read_text()
    else:
        body = _interactive_letter_editor(machine_id, to, subject)

    # ── Build message content ──
    instance_id = f"{machine_id}-{datetime.now().strftime('%Y-%m-%d')}-S?"
    content = _build_message(instance_id, to, subject, body)

    # ── Deliver ──
    if to == "all":
        recipients = [a for a in all_agents if a != machine_id]
        # Also write to broadcast/ so future agents can read it
        BROADCAST.mkdir(exist_ok=True)
        bc_filename = msg_filename(BROADCAST, machine_id, subject)
        (BROADCAST / bc_filename).write_text(content, encoding='utf-8')
        print(f"\n  ✓ Broadcast saved: agents/broadcast/{bc_filename}")
    else:
        recipients = [to]

    for recipient in recipients:
        inbox = AGENTS_DIR / recipient / "inbox"
        inbox.mkdir(parents=True, exist_ok=True)
        filename = msg_filename(inbox, machine_id, subject)
        (inbox / filename).write_text(content, encoding='utf-8')
        print(f"  ✓ Delivered to {recipient}: agents/{recipient}/inbox/{filename}")

    print(f"\nMessage sent. Remember to git add -A && git commit && git push.\n")


def _interactive_letter_editor(sender: str, recipient: str, subject: str) -> str:
    """Guide Claude through writing a structured end-of-session letter."""
    print(f"\n── Composing letter from {sender} to {recipient} ──")
    print("(Press Enter twice to finish each section. Enter blank to skip a section.)\n")

    sections = {}

    prompts = [
        ("progress",    "What did you work on / discover this session?"),
        ("handoff",     "What should the recipient pick up next?"),
        ("court",       "Any court cases needing a response?"),
        ("questions",   "Specific questions for the recipient?"),
        ("new_findings","New theorems, tangents, or mistakes added to the system?"),
        ("warnings",    "Anything to watch out for / flag?"),
        ("free",        "Anything else (free-form)?"),
    ]

    for key, prompt in prompts:
        print(f"{prompt}")
        lines = []
        while True:
            line = input()
            if line == "" and lines and lines[-1] == "":
                lines.pop()
                break
            lines.append(line)
        text = "\n".join(lines).strip()
        if text:
            sections[key] = text

    # Build structured body
    body_parts = []
    label_map = {
        "progress":     "## Progress This Session",
        "handoff":      "## Handoff: What to Pick Up Next",
        "court":        "## Court Cases Needing Response",
        "questions":    "## Questions for You",
        "new_findings": "## New Findings Added to System",
        "warnings":     "## Warnings / Flags",
        "free":         "## Other Notes",
    }
    for key, label in label_map.items():
        if key in sections:
            body_parts.append(f"{label}\n\n{sections[key]}")

    return "\n\n---\n\n".join(body_parts) if body_parts else "(No content provided.)"


def _build_message(from_id: str, to: str, subject: str, body: str) -> str:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    machine_id = get_machine_id()
    return textwrap.dedent(f"""\
        # Message: {subject}

        **From:** {from_id}
        **To:** {to}
        **Sent:** {timestamp}

        ---

        {body}

        ---

        *Reply by writing to `agents/{machine_id}/inbox/` or run `python3 agents/processor.py --send --to {machine_id}`*
    """)


def cmd_register(args):
    """Register this machine as a new agent in the network."""
    machine_id = get_machine_id()
    agent_dir  = AGENTS_DIR / machine_id
    inbox_dir  = agent_dir / "inbox"

    if (agent_dir / "identity.md").exists():
        print(f"Machine '{machine_id}' is already registered at agents/{machine_id}/")
        print("Edit agents/{machine_id}/identity.md if you want to update your description.")
        return

    inbox_dir.mkdir(parents=True, exist_ok=True)

    identity_content = textwrap.dedent(f"""\
        # Agent Identity: {machine_id}

        **Machine ID:** {machine_id}
        **Owner:** [fill in: human owner / account]
        **Description:** [fill in: e.g. "Main desktop, Account A"]
        **Registered:** {datetime.now().strftime("%Y-%m-%d")}
        **Last session:** {datetime.now().strftime("%Y-%m-%d")}

        ---

        ## Notes

        [Add any notes about this machine's role or specialization.]
    """)
    (agent_dir / "identity.md").write_text(identity_content)

    # Add a .gitkeep so the inbox directory is tracked
    (inbox_dir / ".gitkeep").write_text("")

    print(f"\n✓ Registered agent: {machine_id}")
    print(f"  Created: agents/{machine_id}/identity.md")
    print(f"  Created: agents/{machine_id}/inbox/")
    print(f"\nNext steps:")
    print(f"  1. Edit agents/{machine_id}/identity.md with your description")
    print(f"  2. Update agents/REGISTRY.md to add your entry to the table")
    print(f"  3. git add agents/{machine_id}/ && git commit -m 'register agent: {machine_id}'")
    print(f"  4. git push")
    print(f"  5. (Optional) Send a broadcast: python3 agents/processor.py --send --to all")


def cmd_list(args):
    """List all registered agents and their inbox status."""
    machine_id = get_machine_id()
    all_agents = get_all_agents()
    read_log   = load_read_log()

    print(f"\n{'='*60}")
    print(f"  Agent Network ({len(all_agents)} machines)")
    print(f"{'='*60}\n")
    print(f"  {'Machine':<20} {'Unread (you)':<15} {'Total msgs':<12} {'Last session'}")
    print(f"  {'─'*18}   {'─'*13}   {'─'*10}   {'─'*20}")

    for agent in all_agents:
        marker = " ← YOU" if agent == machine_id else ""
        inbox  = AGENTS_DIR / agent / "inbox"
        if agent == machine_id:
            unread = len(get_unread_messages(inbox, read_log))
            total  = len(list(inbox.glob("MSG-*.md"))) if inbox.exists() else 0
        else:
            unread = "—"
            total  = len(list(inbox.glob("MSG-*.md"))) if inbox.exists() else 0

        # Try to read last session date from identity.md
        identity = AGENTS_DIR / agent / "identity.md"
        last = "unknown"
        if identity.exists():
            for line in identity.read_text().splitlines():
                if "Last session" in line:
                    last = line.split(":", 1)[-1].strip().split("—")[0].strip()
                    break

        print(f"  {agent + marker:<22} {str(unread):<15} {str(total):<12} {last}")

    # Broadcast
    bc_total  = len(list(BROADCAST.glob("MSG-*.md"))) if BROADCAST.exists() else 0
    bc_unread = len(get_unread_messages(BROADCAST, read_log))
    print(f"\n  {'broadcast':<22} {str(bc_unread):<15} {str(bc_total):<12} (shared)")
    print()


def cmd_status(args):
    """Full status: run --list then --check."""
    cmd_list(args)
    cmd_check(args)


# ─── CLI ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Multi-agent message processor for math-research network"
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--check",    action="store_true", help="Read incoming messages")
    group.add_argument("--send",     action="store_true", help="Write and send a message")
    group.add_argument("--register", action="store_true", help="Register this machine")
    group.add_argument("--list",     action="store_true", help="List all agents")
    group.add_argument("--status",   action="store_true", help="Full status (list + check)")

    # --send options
    parser.add_argument("--to",        help="Recipient machine name or 'all'")
    parser.add_argument("--subject",   help="Message subject line")
    parser.add_argument("--body",      help="Message body (short, inline)")
    parser.add_argument("--body-file", help="Path to a file containing the message body")

    args = parser.parse_args()

    if args.check:
        cmd_check(args)
    elif args.send:
        cmd_send(args)
    elif args.register:
        cmd_register(args)
    elif args.list:
        cmd_list(args)
    elif args.status:
        cmd_status(args)


if __name__ == "__main__":
    main()
