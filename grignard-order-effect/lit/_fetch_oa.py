#!/usr/bin/env python3
"""Query OpenAlex for OA locations of target DOIs; emit a manifest."""
import json, sys, urllib.parse, urllib.request

MAILTO = "research@example.org"

TARGETS = {
    "c2ob26658j": "10.1039/c2ob26658j",
    "jo026017c": "10.1021/jo026017c",
    "jpcb7b02716": "10.1021/acs.jpcb.7b02716",
    "chemsci2025_grignard": "10.1039/D5SC01078K",
    "blomberg1977": "10.1055/s-1977-24261",
}


def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": f"lit-sweep ({MAILTO})"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)


def by_doi(doi):
    url = f"https://api.openalex.org/works/doi:{urllib.parse.quote(doi)}?mailto={MAILTO}"
    return get(url)


def search(q):
    url = (
        "https://api.openalex.org/works?"
        f"search={urllib.parse.quote(q)}&per_page=5&mailto={MAILTO}"
    )
    return get(url)


def summarize(w):
    best = w.get("best_oa_location") or {}
    locs = []
    for loc in (w.get("locations") or []):
        if loc.get("pdf_url") or loc.get("is_oa"):
            locs.append(
                {
                    "pdf_url": loc.get("pdf_url"),
                    "landing": loc.get("landing_page_url"),
                    "is_oa": loc.get("is_oa"),
                    "version": loc.get("version"),
                    "host": (loc.get("source") or {}).get("display_name"),
                }
            )
    return {
        "id": w.get("doi"),
        "title": w.get("title"),
        "year": w.get("publication_year"),
        "venue": ((w.get("primary_location") or {}).get("source") or {}).get("display_name"),
        "is_oa": w.get("open_access", {}).get("is_oa"),
        "oa_status": w.get("open_access", {}).get("oa_status"),
        "best_pdf": best.get("pdf_url"),
        "best_landing": best.get("landing_page_url"),
        "locations": locs,
        "cited_by": w.get("cited_by_count"),
    }


if __name__ == "__main__":
    out = {}
    for key, doi in TARGETS.items():
        try:
            out[key] = summarize(by_doi(doi))
        except Exception as e:  # noqa: BLE001
            out[key] = {"error": str(e), "doi": doi}
    searches = {
        "cluster_model": "Quantum-chemical study of the Grignard reaction mechanism within the cluster model of reaction center",
        "schlenk_ml_2025": "influence of the organic residue and the solvent in the Schlenk equilibrium for Grignard reagents in THF",
    }
    for key, q in searches.items():
        try:
            res = search(q)
            out[key] = [summarize(w) for w in res.get("results", [])]
        except Exception as e:  # noqa: BLE001
            out[key] = {"error": str(e)}
    print(json.dumps(out, indent=2, ensure_ascii=False))
