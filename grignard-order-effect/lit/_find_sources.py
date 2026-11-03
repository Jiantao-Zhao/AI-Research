#!/usr/bin/env python3
"""Find fulltext sources for the target DOIs across legitimate OA channels."""
import json, urllib.parse, urllib.request

MAILTO = "research@example.org"
DOIS = {
    "peltzer2017": "10.1021/acs.jpcb.7b02716",
    "chemsci2025": "10.1039/d5sc01078k",
    "jacs2020": "10.1021/jacs.9b11829",
    "jcat2025": "10.1016/j.jcat.2025.116619",
    "chen2012": "10.1039/c2ob26658j",
    "yamazaki2002": "10.1021/jo026017c",
    "blomberg1977": "10.1055/s-1977-24261",
    "cluster2013": "10.1016/j.comptc.2013.02.011",
}


def get_json(url):
    req = urllib.request.Request(url, headers={"User-Agent": f"lit ({MAILTO})"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)


def unpaywall(doi):
    u = f"https://api.unpaywall.org/v2/{doi}?email={MAILTO}"
    try:
        d = get_json(u)
    except Exception as e:  # noqa: BLE001
        return {"error": str(e)}
    locs = []
    for loc in d.get("oa_locations") or []:
        locs.append(
            {
                "url": loc.get("url_for_pdf") or loc.get("url"),
                "host": loc.get("host_type"),
                "version": loc.get("version"),
                "repo": loc.get("repository_institution"),
            }
        )
    return {"is_oa": d.get("is_oa"), "status": d.get("oa_status"), "locs": locs}


def europepmc(doi):
    q = urllib.parse.quote(f'DOI:"{doi}"')
    u = (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/search?"
        f"query={q}&format=json&resultType=core"
    )
    try:
        d = get_json(u)
    except Exception as e:  # noqa: BLE001
        return {"error": str(e)}
    out = []
    for hit in d.get("resultList", {}).get("result", []):
        out.append(
            {
                "pmcid": hit.get("pmcid"),
                "isOpen": hit.get("isOpenAccess"),
                "hasPDF": hit.get("hasPDF"),
                "title": hit.get("title"),
            }
        )
    return out


res = {}
for k, doi in DOIS.items():
    res[k] = {"doi": doi, "unpaywall": unpaywall(doi), "europepmc": europepmc(doi)}
print(json.dumps(res, indent=2, ensure_ascii=False))
