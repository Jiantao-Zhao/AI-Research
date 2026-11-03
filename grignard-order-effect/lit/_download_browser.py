#!/usr/bin/env python3
"""Download fulltexts using a real Chrome (Playwright) to pass bot checks.

Only legitimate channels are attempted: publisher OA, PMC, HAL, institutional
repositories. Paywalled articles are reported, not circumvented.
"""
import json, os, sys, time
from playwright.sync_api import sync_playwright

OUT = os.path.dirname(os.path.abspath(__file__))

TARGETS = [
    dict(
        name="Chen2012_OBC_GrignardFormation",
        landing="https://pubs.rsc.org/en/content/articlelanding/2012/ob/c2ob26658j",
        pdfs=["https://pubs.rsc.org/en/content/articlepdf/2012/ob/c2ob26658j"],
    ),
    dict(
        name="Yamazaki2002_JOC_Addition",
        landing="https://pubs.acs.org/doi/10.1021/jo026017c",
        pdfs=["https://pubs.acs.org/doi/pdf/10.1021/jo026017c"],
    ),
    dict(
        name="Peltzer2017_JPCB_Schlenk",
        landing="https://pubs.acs.org/doi/10.1021/acs.jpcb.7b02716",
        pdfs=[
            "https://pubs.acs.org/doi/pdf/10.1021/acs.jpcb.7b02716",
            "https://pubs.acs.org/doi/epdf/10.1021/acs.jpcb.7b02716",
        ],
    ),
    dict(
        name="ChemSci2025_FellowshipReview",
        landing="https://pubs.rsc.org/en/content/articlelanding/2025/sc/d5sc01078k",
        pdfs=["https://pubs.rsc.org/en/content/articlepdf/2025/sc/d5sc01078k"],
    ),
    dict(
        name="JACS2020_GrignardPuzzle",
        landing="https://pubs.acs.org/doi/10.1021/jacs.9b11829",
        pdfs=[
            "https://pubs.acs.org/doi/pdf/10.1021/jacs.9b11829",
            "https://pubs.acs.org/doi/epdf/10.1021/jacs.9b11829",
        ],
    ),
    dict(
        name="JCAT2025_MLP_Schlenk",
        landing="https://hal.science/hal-05441570",
        pdfs=[
            "https://hal.science/hal-05441570/document",
            "https://hal.science/hal-05441570v1/document",
        ],
    ),
    dict(
        name="ClusterModel2013_CTC",
        landing="https://www.sciencedirect.com/science/article/pii/S2210271X13000820",
        pdfs=[
            "https://www.sciencedirect.com/science/article/pii/S2210271X13000820/pdfft"
        ],
    ),
    dict(
        name="Blomberg1977_Synthesis_Barbier",
        landing="https://www.thieme-connect.com/products/ejournals/abstract/10.1055/s-1977-24261",
        pdfs=[],
    ),
]


def settle(page, seconds=6.0):
    """Wait for Cloudflare / Anubis interstitials to clear."""
    deadline = time.time() + 45
    while time.time() < deadline:
        try:
            title = page.title()
        except Exception:  # noqa: BLE001
            title = ""
        low = title.lower()
        if "just a moment" not in low and "bot" not in low and "attention" not in low:
            break
        time.sleep(1.5)
    page.wait_for_timeout(int(seconds * 1000))


def main(headless=True):
    report = {}
    with sync_playwright() as p:
        browser = p.chromium.launch(
            channel="chrome",
            headless=headless,
            args=["--disable-blink-features=AutomationControlled"],
        )
        ctx = browser.new_context(
            viewport={"width": 1440, "height": 900},
            user_agent=(
                "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
                "(KHTML, like Gecko) Chrome/156.0.0.0 Safari/537.36"
            ),
            accept_downloads=True,
        )
        page = ctx.new_page()
        for t in TARGETS:
            entry = {"landing_status": None, "attempts": [], "ok": False, "file": None}
            try:
                resp = page.goto(t["landing"], wait_until="domcontentloaded", timeout=60000)
                entry["landing_status"] = resp.status if resp else None
                settle(page)
                entry["final_url"] = page.url
                entry["title"] = page.title()[:120]
                # Does the page advertise open access?
                body = page.content()[:400000]
                entry["shows_open_access"] = any(
                    k in body for k in ("Open Access", "open access", "OPEN ACCESS")
                )
            except Exception as e:  # noqa: BLE001
                entry["landing_error"] = str(e)[:200]

            for url in t["pdfs"]:
                try:
                    r = ctx.request.get(url, timeout=90000, headers={"Referer": page.url})
                    ct = r.headers.get("content-type", "")
                    body = r.body()
                    entry["attempts"].append(
                        {"url": url, "status": r.status, "ctype": ct.split(";")[0],
                         "bytes": len(body), "magic": body[:4].decode("latin-1")}
                    )
                    if body[:4] == b"%PDF":
                        path = os.path.join(OUT, t["name"] + ".pdf")
                        with open(path, "wb") as f:
                            f.write(body)
                        entry.update(ok=True, file=os.path.basename(path))
                        break
                except Exception as e:  # noqa: BLE001
                    entry["attempts"].append({"url": url, "error": str(e)[:200]})
            report[t["name"]] = entry
            print(json.dumps({t["name"]: entry}, ensure_ascii=False), flush=True)
        browser.close()
    with open(os.path.join(OUT, "_download_report.json"), "w") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    return report


if __name__ == "__main__":
    main(headless="--headed" not in sys.argv)
