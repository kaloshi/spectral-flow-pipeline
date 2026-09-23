#!/usr/bin/env python3
"""Pre-publication check for this repository.

Scans every file below a directory for content that must not be published:
  * terms from a private deny list (one term per line; patient or sample identifiers, personal names, project names).
    The deny list is kept outside the repository and passed with --deny.
  * absolute paths of local machines, e-mail addresses and sample-name patterns,
  * study groups or patient information, also as anonymised codes (sex, diagnosis, genotype, factor codes such as A0/B1),
  * data files (FCS, workspaces, spreadsheets, tables, R data) and files larger than 1 MB,
  * German text in code comments or documentation (the repository is written in English),
  * numbers that look like study results (p values, effect sizes) in comments or documentation, for manual review.

Errors make the script exit with status 1. Warnings are listed for manual review but do not fail the check.

Usage: python tools/check_public.py <directory> [--deny private_terms.txt] [--allow-name "Name"]
"""
import argparse
import re
import sys
from pathlib import Path

DATA_SUFFIXES = {".fcs", ".wsp", ".acs", ".xlsx", ".xls", ".xlsb", ".csv", ".tsv", ".gz", ".rds", ".rdata", ".parquet",
                 ".h5", ".h5ad", ".pdf", ".svg", ".png", ".jpg", ".zip"}
TEXT_SUFFIXES = {".r", ".py", ".sh", ".md", ".txt", ".cff", ".yml", ".yaml", ".toml", ".cfg", ""}
PATH_PATTERN = re.compile(r"(/Volumes/|/Users/|/home/|/run/media/|/private/tmp/|[A-Z]:\\\\)")
MAIL_PATTERN = re.compile(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}")
SAMPLE_PATTERN = re.compile(r"\b[A-Z]_ID\d{3,4}\b|\bID(?!7000)\d{4}\b")   # ID7000 = instrument name
GERMAN_PATTERN = re.compile(r"[äöüÄÖÜß]|\b(und|nicht|oder|wird|werden|sind|ist|der|die|das|mit|fuer|für|auf|bei|nach)\b")
# Study groups and patient information must not appear, not even as anonymised codes.
GROUP_PATTERN = re.compile(r"\b(sex|sexes|female|male|women|men|woman|man|gender|patient\w*|cohort\w*|donor\w*|"
                           r"genotype\w*|diagnos\w*|disease\w*|tumou?r\w*|cancer\w*|age_decade|factor_a|factor_b|"
                           r"[AB][01])\b", re.IGNORECASE)   # study-specific terms belong in the private deny list
RESULT_PATTERN = re.compile(r"\b(p|q)\s*[=<]\s*0[.,]\d{2,}|\bdelta\s*[=:]\s*[-+−]?0[.,]\d{2}|\brho\s*[=:]\s*[-+−]?0[.,]\d{2}")


def comment_lines(path, text):
    """Lines that are comments or documentation (German and result checks apply only there)."""
    if path.suffix.lower() in {".md", ".txt", ".cff", ""}:
        return list(enumerate(text.splitlines(), 1))
    out, in_doc = [], False
    for i, line in enumerate(text.splitlines(), 1):
        s = line.strip()
        if path.suffix.lower() == ".py" and s.startswith(('"""', "'''")):
            in_doc = not in_doc if s.count('"""') + s.count("'''") == 1 else in_doc
            out.append((i, line)); continue
        if in_doc or s.startswith("#"):
            out.append((i, line))
        elif "#" in line and path.suffix.lower() in {".r", ".py", ".sh"}:
            out.append((i, line[line.index("#"):]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("directory")
    ap.add_argument("--deny", help="private file with forbidden terms, one per line")
    ap.add_argument("--allow-name", action="append", default=[], help="names allowed in LICENSE, README and CITATION")
    ap.add_argument("--allow-term", action="append", default=[],
                    help="deny-list term allowed everywhere because it has a technical meaning here (e.g. a random seed)")
    a = ap.parse_args()
    root = Path(a.directory)
    deny = []
    if a.deny:
        deny = [t.strip() for t in Path(a.deny).read_text(encoding="utf-8").splitlines() if t.strip() and not t.startswith("#")]
    errors, warnings = [], []
    for f in sorted(p for p in root.rglob("*") if p.is_file() and ".git" not in p.parts):
        rel = f.relative_to(root)
        suf = f.suffix.lower()
        if suf in DATA_SUFFIXES:
            errors.append(f"{rel}: data or figure file ({suf})")
            continue
        if f.stat().st_size > 1_000_000:
            errors.append(f"{rel}: larger than 1 MB")
        if suf not in TEXT_SUFFIXES:
            warnings.append(f"{rel}: unknown file type, check manually")
            continue
        text = f.read_text(encoding="utf-8", errors="replace")
        allow_names = f.name in {"LICENSE", "README.md", "CITATION.cff"}
        for term in deny:
            if (allow_names and term in a.allow_name) or term in a.allow_term:
                continue
            for m in re.finditer(r"(?<![A-Za-z0-9])" + re.escape(term) + r"(?![A-Za-z0-9])", text, flags=re.IGNORECASE):
                line = text.count("\n", 0, m.start()) + 1
                errors.append(f"{rel}:{line}: private term '{term}'")
        if f.resolve() == Path(__file__).resolve():
            continue   # this file defines the patterns; the deny list above still applies to it
        for pat, what in ((PATH_PATTERN, "absolute path"), (MAIL_PATTERN, "e-mail address"), (SAMPLE_PATTERN, "sample identifier"),
                          (GROUP_PATTERN, "group or patient term")):
            for m in pat.finditer(text):
                line = text.count("\n", 0, m.start()) + 1
                errors.append(f"{rel}:{line}: {what} '{m.group(0)}'")
        for i, line in comment_lines(f, text):
            if GERMAN_PATTERN.search(line):
                warnings.append(f"{rel}:{i}: German text? {line.strip()[:90]}")
            if RESULT_PATTERN.search(line):
                warnings.append(f"{rel}:{i}: looks like a result, check: {line.strip()[:90]}")
    for e in errors:
        print("ERROR   " + e)
    for w in warnings:
        print("WARNING " + w)
    print(f"\n{len(errors)} errors, {len(warnings)} warnings in {root}")
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()
