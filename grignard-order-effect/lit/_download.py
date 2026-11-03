#!/usr/bin/env python3
"""Try to download fulltext PDFs from legitimate channels; record what worked."""
import json, os, urllib.parse, urllib.request

UA = (
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/120.0 Safari/537.36"
)
MAILTO = "research@example.org"
OUT = os.path.dirname(os.path.abspath(__file__))

CANDIDATES = {
    "Chen2012_OBC_GrignardFormation": [
        "https://pubs.rsc.org/en/content/articlepdf/2012/ob/c2ob26658j",
        "https://pubs.rsc.org/en/content/articlehtml/2012/ob/c2ob26658j",
    ],
    "Yamazaki2002_JOC_Addition": [
        "https://pubs.acs.org/doi/pdf/10.1021/jo026017c",
    ],
    "Peltzer2017_JPCB_Schlenk": [
        "https://pubs.acs.org/doi/pdf/10.1021/acs.jpcb.7b02716",
        "http://hdl.handle.net/10852/64642",
    ],
    "ChemSci2025_FellowshipReview": [
        "https://europepmc.org/articles/PMC12039541?pdf=render",
        "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12039541/fullTextPDF",
        "https://hal.umontpellier.fr/hal-05069103v1/document",
        "https://pubs.rsc.org/en/content/articlepdf/2025/sc/d5sc01078k",
    ],
    "JACS2020_GrignardPuzzle": [
        "https://pubs.acs.org/doi/pdf/10.1021/jacs.9b11829",
        "http://urn.nb.no/URN:NBN:no-86650",
    ],
    "JCAT2025_MLP_Schlenk": [
        "https://hal.umontpellier.fr/hal-05441570v1/document",
    ],
    "ClusterModel2013_CTC": [
        "https://www.sciencedirect.com/science/article/pii/S2210271X13000820",
    ],
    "Blomberg1977_Synthesis_Barbier": [
        "https://www.thieme-connect.com/products/ejournals/abstract/10.1055/s-1977-24261",
    ],
}


def fetch(url):
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": UA,
            "Accept": "application/pdf,text/html;q=0.9,*/*;q=0.8",
            "Accept-Language": "en-US,en;q=0.9",
        },
    )
    with urllib.request.urlopen(req, timeout=60) as r:
        ctype = r.headers.get("Content-Type", "")
        data = r.read()
    return ctype, data


def semantic_scholar(doi):
    u = (
        "https://api.semanticscholar.org/graph/v1/paper/DOI:"
        f"{urllib.parse.quote(doi)}?fields=title,openAccessPdf,externalIds"
    )
    try:
        req = urllib.request.Request(u, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=30) as r:
            return json.load(r)
    except Exception as e:  # noqa: BLE001
        return {"error": str(e)}


if __name__ == "__main__":
    report = {}
    for name, urls in CANDIDATES.items():
        entry = {"tried": [], "ok": False, "file": None, "note": None}
        for url in urls:
            try:
                ctype, data = fetch(url)
            except Exception as e:  # noqa: BLE001
                entry["tried"].append({"url": url, "error": str(e)[:120]})
                continue
            entry["tried"].append(
                {"url": url, "ctype": ctype.split(";")[0], "bytes": len(data)}
            )
            if "pdf" in ctype.lower() and data[:4] == b"%PDF":
                path = os.path.join(OUT, name + ".pdf")
                with open(path, "wb") as f:
                    f.write(data)
                entry.update(ok=True, file=os.path.basename(path))
                break
        report[name] = entry
    print(json.dumps(report, indent=2))
