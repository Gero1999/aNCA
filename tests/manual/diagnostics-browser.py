"""Manual browser integration check using an installed Python Playwright and Chromium.

From the checkout, activate R and run:
    python3 tests/manual/diagnostics-browser.py --rscript /path/to/Rscript

The script starts only the synthetic app in diagnostics-fixture.R, uses two browser
contexts, and terminates its own server. It never loads clinical datasets. Evidence
and downloaded reports are written outside the checkout.
"""

import argparse
import json
import os
import re
from pathlib import Path
import subprocess
import time
from urllib.request import urlopen

from playwright.sync_api import expect, sync_playwright


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--chromium", default="/usr/bin/chromium")
    parser.add_argument("--port", type=int, default=3849)
    parser.add_argument("--evidence", default="/tmp/anca-diagnostics-browser")
    args = parser.parse_args()
    evidence = Path(args.evidence).resolve()
    evidence.mkdir(parents=True, exist_ok=True)
    root = Path(__file__).resolve().parents[2]
    env = dict(os.environ, DIAGNOSTICS_TEST_PORT=str(args.port),
               DIAGNOSTICS_TEST_LOG=str(evidence / "logger.log"))
    env.setdefault("XDG_CACHE_HOME", str(evidence / "cache"))
    outcomes = []
    errors = []

    def passed(name):
        outcomes.append(name)
        print("PASS:", name, flush=True)

    with (evidence / "server.log").open("w") as log:
        server = subprocess.Popen(
            [args.rscript, "tests/manual/diagnostics-fixture.R"], cwd=root,
            env=env, stdout=log, stderr=subprocess.STDOUT,
        )
        try:
            url = f"http://127.0.0.1:{args.port}/"
            for _ in range(225):
                if server.poll() is not None:
                    raise RuntimeError("Synthetic server exited; see server.log")
                try:
                    with urlopen(url, timeout=1) as response:
                        if response.status == 200:
                            break
                except OSError:
                    time.sleep(0.2)
            else:
                raise RuntimeError("Synthetic server did not become ready")

            with sync_playwright() as playwright:
                browser = playwright.chromium.launch(
                    executable_path=args.chromium, headless=True,
                    args=["--no-sandbox", "--disable-dev-shm-usage"],
                )
                first_context = browser.new_context(accept_downloads=True)
                second_context = browser.new_context(accept_downloads=True)
                first = first_context.new_page()
                second = second_context.new_page()
                for page in (first, second):
                    page.on("pageerror", lambda error: errors.append(str(error)))
                    page.goto(url)
                    page.wait_for_function(
                        "window.Shiny && Shiny.shinyapp && "
                        "Shiny.shinyapp.$socket && Shiny.shinyapp.$socket.readyState === 1"
                    )
                    expect(page.locator("#fixture_loading")).to_have_text("Idle")

                def open_panel(page):
                    if not page.locator("#diagnostics-panel").is_visible():
                        page.locator("#diagnostics-open").click()
                    expect(page.locator("#diagnostics-panel")).to_be_visible()

                def report(page, name):
                    open_panel(page)
                    with page.expect_download() as downloaded:
                        page.locator("#diagnostics-download").click()
                    download = downloaded.value
                    assert download.url.startswith("blob:"), download.url
                    target = evidence / name
                    download.save_as(target)
                    result = json.loads(target.read_text())
                    assert isinstance(result.get("entries"), list), result
                    return result

                first.locator("#fixture_run").click()
                expect(first.locator("#fixture_result")).to_have_text("Synthetic result: 42")
                expect(first.locator("#fixture_loading")).to_have_text("Idle")
                passed("normal synthetic operation")

                first.locator("#fixture_invalid").check()
                first.locator("#fixture_run").click()
                expect(first.locator("#fixture_result")).to_have_text("No current result")
                expect(first.locator("#fixture_loading")).to_have_text("Idle")
                expect(first.locator("#fixture_run")).to_be_enabled()
                notification = first.locator(".shiny-notification").filter(
                    has_text=re.compile("diagnostics", re.IGNORECASE)
                )
                expect(notification.first).to_be_visible()
                notification.first.locator("button,a").first.click()
                expect(first.locator("#diagnostics-panel")).to_be_visible()
                failed = report(first, "failed-operation.json")
                failure_entries = [e for e in failed["entries"] if e["level"] == "ERROR"]
                assert failure_entries, failed
                reference = failure_entries[-1]["id"]
                assert reference
                serialized = json.dumps(failed)
                assert all(s not in serialized for s in
                           ("do-not-export", "/srv/private", "ALPHA-42"))
                passed("recoverable failure, stale-result invalidation, cleanup and diagnostics route")

                other = report(second, "second-session.json")
                assert reference not in json.dumps(other)
                assert not any(e["level"] == "ERROR" for e in other["entries"])
                passed("concurrent-session isolation and safe report fields")

                first.locator("#diagnostics-close").click()
                first.locator("#fixture_invalid").uncheck()
                first.locator("#fixture_run").click()
                expect(first.locator("#fixture_result")).to_have_text("Synthetic result: 42")
                expect(first.locator("#fixture_loading")).to_have_text("Idle")
                expect(first.locator("#fixture_run")).to_be_enabled()
                passed("retry succeeds in the same live session")

                first.locator("#fixture_literal").click()
                open_panel(first)
                expect(first.locator("#diagnostics-events")).to_contain_text("<img")
                assert first.locator("#diagnostics-events img").count() == 0
                assert not first.evaluate("Boolean(window.diagnosticsInjected)")
                passed("HTML-like diagnostic messages render literally")
                first.locator("#diagnostics-close").click()

                first.locator("#fixture_flood").click()
                first.wait_for_timeout(200)
                bounded = report(first, "bounded-buffer.json")
                assert 1 <= len(bounded["entries"]) <= 100
                assert all(len(str(value)) <= 500 for entry in bounded["entries"]
                           for value in entry.values())
                assert len(json.dumps(bounded)) <= 150_000
                passed("bounded browser buffer and bounded entry fields")
                first.locator("#diagnostics-close").click()

                first.locator("#fixture_disconnect").click()
                first.wait_for_function(
                    "!Shiny.shinyapp.$socket || Shiny.shinyapp.$socket.readyState !== 1"
                )
                expect(first.locator("#diagnostics-status")).to_contain_text("disconnect")
                first.locator("#diagnostics-close").click()
                open_panel(first)
                offline_requests = []
                first.on("request", lambda request: offline_requests.append(request.url))
                offline = report(first, "after-disconnect.json")
                assert offline["disconnected"] is True
                assert offline["entries"] == bounded["entries"]
                assert not any("download" in url and not url.startswith("blob:")
                               for url in offline_requests)
                assert "before disconnection" in json.dumps(offline).lower()
                passed("client download after actual Shiny websocket disconnection")
                assert not errors, errors
                passed("no uncaught browser JavaScript errors")
                browser.close()
        finally:
            server.terminate()
            try:
                server.wait(timeout=10)
            except subprocess.TimeoutExpired:
                server.kill()
                server.wait(timeout=10)
            (evidence / "results.json").write_text(json.dumps(
                {"passed": outcomes, "browser_errors": errors}, indent=2
            ))


if __name__ == "__main__":
    main()
