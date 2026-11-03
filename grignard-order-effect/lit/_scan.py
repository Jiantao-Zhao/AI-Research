#!/usr/bin/env python3
"""Systematic OpenAlex scan for prior work on order-of-addition / speciation."""
import json, urllib.parse, urllib.request

MAILTO = "research@example.org"
QUERIES = [
    ("order_of_addition", 'grignard "order of addition"'),
    ("barbier_mechanism", "Barbier reaction mechanism Grignard one-pot"),
    ("barbier_kinetics", "Barbier Grignard kinetics selectivity comparison"),
    ("schlenk_kinetics", "Schlenk equilibrium kinetics Grignard speciation"),
    ("grignard_aggregation", "Grignard reagent aggregation structure solution THF"),
    ("insitu_grignard", "in situ Grignard formation selectivity"),
    ("grignard_surface", "magnesium surface Grignard formation electron transfer mechanism"),
    ("microkinetic_om", "microkinetic modeling organometallic catalysis speciation"),
    ("grignard_aging", "aging Grignard solution reactivity changes"),
    ("mixing_order", "mixing order addition sequence selectivity reaction outcome"),
]


def q(term, per_page=8):
    url = (
        "https://api.openalex.org/works?search="
        + urllib.parse.quote(term)
        + f"&per_page={per_page}&sort=relevance_score:desc&mailto={MAILTO}"
    )
    req = urllib.request.Request(url, headers={"User-Agent": f"lit ({MAILTO})"})
    with urllib.request.urlopen(req, timeout=40) as r:
        return json.load(r)


out = {}
for key, term in QUERIES:
    try:
        d = q(term)
        rows = []
        for w in d.get("results", []):
            rows.append(
                {
                    "year": w.get("publication_year"),
                    "title": w.get("title"),
                    "doi": (w.get("doi") or "").replace("https://doi.org/", ""),
                    "venue": ((w.get("primary_location") or {}).get("source") or {}).get(
                        "display_name"
                    ),
                    "oa": w.get("open_access", {}).get("is_oa"),
                    "cited": w.get("cited_by_count"),
                }
            )
        out[key] = rows
    except Exception as e:  # noqa: BLE001
        out[key] = [{"error": str(e)}]

print(json.dumps(out, indent=1, ensure_ascii=False))
