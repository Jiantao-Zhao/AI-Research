#!/usr/bin/env python3
"""Second download wave: additional OA sources relevant to the design."""
import json, os, time
from playwright.sync_api import sync_playwright

OUT = os.path.dirname(os.path.abspath(__file__))

TARGETS = [
    dict(
        name="AccChemRes2026_DynamicSpeciation_Zn",
        landing="https://pubs.acs.org/doi/10.1021/acs.accounts.6c00598",
        pdfs=[
            "https://pubs.acs.org/doi/pdf/10.1021/acs.accounts.6c00598",
            "https://pubs.acs.org/achre4/article-pdf/doi/10.1021/acs.accounts.6c00598/68729164/acs.accounts.6c00598.pdf",
        ],
    ),
    dict(
        name="PAC1980_MechanismGrignardKetones",
        landing="https://www.degruyter.com/document/doi/10.1351/pac198052030545/html",
        pdfs=["https://www.degruyter.com/document/doi/10.1351/pac198052030545/pdf"],
    ),
    dict(
        name="RCE2023_PackedBed_ResidenceTime",
        landing="https://pubs.rsc.org/en/content/articlelanding/2023/re/d3re00191a",
        pdfs=["https://pubs.rsc.org/en/content/articlepdf/2023/re/d3re00191a"],
    ),
]


def settle(page):
    deadline = time.time() + 45
    while time.time() < deadline:
        try:
            t = page.title().lower()
        except Exception:
            t = ""
        if "just a moment" not in t and "bot" not in t:
            break
        time.sleep(1.5)
    page.wait_for_timeout(6000)


def run(headless):
    rep = {}
    with sync_playwright() as p:
        b = p.chromium.launch(channel="chrome", headless=headless,
                              args=["--disable-blink-features=AutomationControlled"])
        ctx = b.new_context(viewport={"width": 1440, "height": 900}, accept_downloads=True,
                            user_agent=("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                                        "AppleWebKit/537.36 (KHTML, like Gecko) "
                                        "Chrome/156.0.0.0 Safari/537.36"))
        page = ctx.new_page()
        for t in TARGETS:
            e = {"ok": False, "file": None, "attempts": []}
            try:
                page.goto(t["landing"], wait_until="domcontentloaded", timeout=60000)
                settle(page)
                e["final_url"] = page.url[:120]
            except Exception as ex:  # noqa: BLE001
                e["landing_error"] = str(ex)[:150]
            for url in t["pdfs"]:
                try:
                    r = ctx.request.get(url, timeout=90000, headers={"Referer": page.url})
                    body = r.body()
                    e["attempts"].append({"url": url[:80], "status": r.status,
                                          "ctype": r.headers.get("content-type", "").split(";")[0],
                                          "bytes": len(body)})
                    if body[:4] == b"%PDF":
                        open(os.path.join(OUT, t["name"] + ".pdf"), "wb").write(body)
                        e.update(ok=True, file=t["name"] + ".pdf")
                        break
                except Exception as ex:  # noqa: BLE001
                    e["attempts"].append({"url": url[:80], "error": str(ex)[:120]})
            rep[t["name"]] = e
            print(json.dumps({t["name"]: e}, ensure_ascii=False), flush=True)
        b.close()
    return rep


if __name__ == "__main__":
    import sys
    r = run(headless="--headed" not in sys.argv)
    failed = [k for k, v in r.items() if not v["ok"]]
    if failed:
        print("retrying headed:", failed)
        run(headless=False)
