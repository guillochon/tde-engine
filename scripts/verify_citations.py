#!/usr/bin/env python3
"""
Pull titles and abstracts from NASA ADS for every reference cited in the paper,
so the citations can be checked for existence and relevance.

    export ADS_TOKEN=...              # https://ui.adsabs.harvard.edu/user/settings/token
    python scripts/verify_citations.py

    python scripts/verify_citations.py --dry-run     # no network; show what it would query
    python scripts/verify_citations.py --all         # include uncited bib entries
    python scripts/verify_citations.py -o out.md     # choose output path

Resolution order for each entry: bibcode (from adsurl) -> DOI -> arXiv eprint ->
first-author + year + title search. Anything that fails all four is reported at
the end under UNRESOLVED so it can be checked by hand.

Output is markdown: one block per reference with the key, what was queried, the
ADS title, the author list, year, and the abstract.
"""

import argparse
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

API = "https://api.adsabs.harvard.edu/v1/search/query"
FIELDS = "bibcode,title,author,year,pub,doi,abstract,citation_count"
HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.join(HERE, "..", "paper")


# ----------------------------------------------------------------- parsing --
def cited_keys(texpath):
    """Every key appearing in any \\cite... command, in order of first use."""
    src = open(texpath, encoding="utf-8").read()
    src = re.sub(r"(?<!\\)%.*", "", src)          # strip comments, keep \%
    keys, seen = [], set()
    pat = r"\\[Cc]ite[a-zA-Z]*\s*(?:\[[^\]]*\]\s*)*\{([^}]+)\}"
    for m in re.finditer(pat, src):
        for k in m.group(1).split(","):
            k = k.strip()
            if k and k not in seen:
                seen.add(k)
                keys.append(k)
    return keys


def parse_bib(bibpath):
    """key -> {field: value}. Tolerant of nested braces in values."""
    src = open(bibpath, encoding="utf-8").read()
    out = {}
    for m in re.finditer(r"@(\w+)\s*\{\s*([^,]+),", src):
        key = m.group(2).strip()
        i = src.index("{", m.start())
        depth, j = 0, i
        while j < len(src):
            if src[j] == "{":
                depth += 1
            elif src[j] == "}":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        body = src[i + 1:j]
        fields = {}
        for fm in re.finditer(r"(\w+)\s*=\s*(\{)", body):
            name = fm.group(1).lower()
            k0 = fm.end() - 1
            d, k = 0, k0
            while k < len(body):
                if body[k] == "{":
                    d += 1
                elif body[k] == "}":
                    d -= 1
                    if d == 0:
                        break
                k += 1
            fields[name] = body[k0 + 1:k].strip()
        for fm in re.finditer(r'(\w+)\s*=\s*"([^"]*)"', body):
            fields.setdefault(fm.group(1).lower(), fm.group(2).strip())
        # bare values, e.g. `year = 1946,` -- common in hand-written entries
        for fm in re.finditer(r"(\w+)\s*=\s*([^\s{}\",][^,\n]*)", body):
            fields.setdefault(fm.group(1).lower(), fm.group(2).strip().rstrip(","))
        out[key] = fields
    return out


def clean(s):
    """Strip LaTeX noise from a bib field so it can go into a query."""
    s = re.sub(r"\\[a-zA-Z]+\s*", " ", s)
    s = s.replace("{", "").replace("}", "").replace("\\", "")
    s = re.sub(r"\s+", " ", s)
    return s.strip()


def query_for(key, fields):
    """(description, ADS query string) using the most reliable identifier."""
    adsurl = fields.get("adsurl", "")
    m = re.search(r"/abs/([^/]+)", adsurl)
    if m:
        bc = urllib.parse.unquote(m.group(1))
        return f"bibcode {bc}", f'bibcode:"{bc}"'
    if fields.get("bibcode"):
        return f"bibcode {fields['bibcode']}", f'bibcode:"{fields["bibcode"]}"'
    doi = fields.get("doi", "")
    if doi and not doi.startswith("10.48550"):
        return f"doi {doi}", f'doi:"{doi}"'
    ep = fields.get("eprint", "")
    if ep:
        return f"arXiv {ep}", f'arxiv:"{ep}"'
    if doi:
        return f"doi {doi}", f'doi:"{doi}"'
    au = clean(fields.get("author", "")).split(" and ")[0]
    au = au.split(",")[0].strip()
    yr = clean(fields.get("year", ""))
    ti = clean(fields.get("title", ""))
    if au and ti:
        q = f'author:"{au}"'
        if yr:
            q += f" year:{yr}"
        q += f' title:"{ti[:80]}"'
        return f"author/title: {au} {yr}", q
    return None, None


# -------------------------------------------------------------------- ADS --
def ads_get(query, token, rows=5, retries=3):
    url = API + "?" + urllib.parse.urlencode(
        {"q": query, "fl": FIELDS, "rows": rows})
    req = urllib.request.Request(url, headers={
        "Authorization": f"Bearer {token}",
        "User-Agent": "verify_citations/1.0",
    })
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=45) as r:
                remaining = r.headers.get("X-RateLimit-Remaining")
                return json.loads(r.read().decode())["response"]["docs"], remaining
        except urllib.error.HTTPError as e:
            if e.code == 429:
                sys.stderr.write("rate limited; sleeping 60s\n")
                time.sleep(60)
                continue
            if e.code == 401:
                sys.exit("ADS rejected the token (401). Check $ADS_TOKEN.")
            if attempt == retries - 1:
                return None, f"HTTP {e.code}"
            time.sleep(2 ** attempt)
        except Exception as e:                     # noqa: BLE001
            if attempt == retries - 1:
                return None, str(e)
            time.sleep(2 ** attempt)
    return None, "failed"


# ------------------------------------------------------------------- main --
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tex", default=os.path.join(PAPER, "engine.tex"))
    ap.add_argument("--bib", default=os.path.join(PAPER, "engineNotes.bib"))
    ap.add_argument("-o", "--out", default="citations_ads.md")
    ap.add_argument("--all", action="store_true",
                    help="include bib entries that are never cited")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the queries without contacting ADS")
    ap.add_argument("--delay", type=float, default=0.15)
    args = ap.parse_args()

    bib = parse_bib(args.bib)
    keys = cited_keys(args.tex)
    missing = [k for k in keys if k not in bib]
    if args.all:
        keys += [k for k in bib if k not in keys]

    print(f"{len(bib)} bib entries; {len(keys)} to resolve", file=sys.stderr)
    if missing:
        print(f"WARNING: cited but absent from the .bib: {missing}", file=sys.stderr)

    if args.dry_run:
        for k in keys:
            desc, q = query_for(k, bib.get(k, {}))
            print(f"{k:28s} {desc or 'NO IDENTIFIER'}")
        return

    token = os.environ.get("ADS_TOKEN") or os.environ.get("ADS_DEV_KEY")
    if not token:
        sys.exit("Set ADS_TOKEN (https://ui.adsabs.harvard.edu/user/settings/token)")

    blocks, unresolved = [], []
    for n, k in enumerate(keys, 1):
        f = bib.get(k, {})
        desc, q = query_for(k, f)
        if not q:
            unresolved.append((k, "no usable identifier in the .bib"))
            continue
        docs, note = ads_get(q, token)
        if not docs:
            unresolved.append((k, f"{desc} -> no ADS match ({note})"))
            continue
        d = docs[0]
        au = d.get("author", [])
        aus = "; ".join(au[:4]) + (" et al." if len(au) > 4 else "")
        blocks.append(
            f"## {k}\n\n"
            f"- **queried by:** {desc}\n"
            f"- **bibcode:** {d.get('bibcode','?')}\n"
            f"- **bib title:** {clean(f.get('title','(none)'))}\n"
            f"- **ADS title:** {(d.get('title') or ['(none)'])[0]}\n"
            f"- **authors:** {aus}\n"
            f"- **year / pub:** {d.get('year','?')} / {d.get('pub','?')}\n"
            f"- **citations:** {d.get('citation_count','?')}\n\n"
            f"**Abstract:** {d.get('abstract','(no abstract in ADS)')}\n"
        )
        print(f"  [{n}/{len(keys)}] {k}", file=sys.stderr)
        time.sleep(args.delay)

    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write("# ADS records for citations in engine.tex\n\n")
        fh.write(f"Resolved {len(blocks)} of {len(keys)}.\n\n")
        if missing:
            fh.write("## Cited but absent from the .bib\n\n"
                     + "\n".join(f"- {k}" for k in missing) + "\n\n")
        fh.write("\n".join(blocks))
        if unresolved:
            fh.write("\n\n# UNRESOLVED - check these by hand\n\n")
            for k, why in unresolved:
                fh.write(f"- **{k}**: {why}\n")
    print(f"\nwrote {args.out} ({len(blocks)} resolved, "
          f"{len(unresolved)} unresolved)", file=sys.stderr)


if __name__ == "__main__":
    main()
