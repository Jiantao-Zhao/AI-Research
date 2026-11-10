"""从 ChemRxiv / Wiley 抓 SI。"""
import re, sys, pathlib
from playwright.sync_api import sync_playwright

URLS = [
    "https://chemrxiv.org/engage/chemrxiv/article-details/64467bfa83fa35f8f62d8a81",
    "https://onlinelibrary.wiley.com/doi/10.1002/anie.202305775",
]
out = pathlib.Path("si_download")

with sync_playwright() as p:
    b = p.chromium.launch(channel="chrome", headless=True)
    ctx = b.new_context(viewport={"width": 1400, "height": 900})
    pg = ctx.new_page()
    for u in URLS:
        try:
            pg.goto(u, wait_until="domcontentloaded", timeout=60000)
            pg.wait_for_timeout(6000)
            html = pg.content()
            print(f"[{u[:50]}] 标题: {pg.title()[:70]}")
            links = set(re.findall(r'href="([^"]+\.pdf[^"]*)"', html))
            links |= set(re.findall(r'href="([^"]*(?:supp|si|misc)[^"]*)"', html, re.I))
            for l in sorted(links)[:12]:
                print("   ", l[:160])
            # 直接找 pdf 资源
            for l in sorted(links):
                if any(k in l.lower() for k in ("supp", "si", "misc", "supporting")):
                    full = l if l.startswith("http") else "https://chemrxiv.org" + l
                    try:
                        r = ctx.request.get(full, timeout=60000)
                        if r.ok and "pdf" in r.headers.get("content-type", ""):
                            data = r.body()
                            if data[:4] == b"%PDF" and len(data) > 100000:
                                f = out.with_suffix(".pdf")
                                f.write_bytes(data)
                                print(f"     -> 保存 {f} ({len(data)/1e6:.2f} MB)")
                    except Exception as e:
                        print("     dl fail", str(e)[:60])
        except Exception as e:
            print("goto fail", u[:40], str(e)[:80])
    b.close()
