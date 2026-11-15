/* Session-only memory. No raw server log, persistent storage, or live download handler. */
(function () {
  "use strict";
  let snapshot = { package_version: "", entries: [] };
  let disconnected = false;
  const text = (value, max) => typeof value === "string" ? value.slice(0, max) : "";
  function report() {
    return {
      notice: "Only information received before disconnection is available. " +
        "This report cannot diagnose process termination, out-of-memory or native-code crashes.",
      disconnected: disconnected,
      package_version: snapshot.package_version,
      entries: snapshot.entries
    };
  }
  function render() {
    const panel = document.getElementById("diagnostics-panel");
    if (!panel) return;
    const errors = snapshot.entries.filter(entry => entry.level === "ERROR");
    const last = errors[errors.length - 1];
    document.getElementById("diagnostics-summary").textContent = last ?
      last.summary + "\nReference: " + last.id : "No error has been received in this session.";
    document.getElementById("diagnostics-events").textContent = snapshot.entries.map(entry =>
      entry.time + " " + entry.level + " " + entry.operation + "\n" + entry.summary +
      "\nReference: " + entry.id).join("\n\n");
    document.getElementById("diagnostics-status").textContent = disconnected ?
      "Disconnected. Only information received before disconnection is available. " +
      "Download this report before reloading the page." :
      "Recent session events only (up to 100). Downloads are generated in this browser.";
  }
  function open() {
    render();
    const panel = document.getElementById("diagnostics-panel");
    if (panel && !panel.open) panel.showModal();
  }
  document.addEventListener("click", function (event) {
    const button = event.target.closest("button");
    if (!button) return;
    if (button.id === "diagnostics-open" || button.dataset.ancaDiagnostics === "open") open();
    if (button.id === "diagnostics-close") document.getElementById("diagnostics-panel").close();
    if (button.id === "diagnostics-download") {
      const url = URL.createObjectURL(new Blob([JSON.stringify(report(), null, 2)],
        { type: "application/json" }));
      const link = document.createElement("a");
      link.href = url;
      link.download = "aNCA-diagnostics.json";
      document.body.appendChild(link);
      link.click();
      link.remove();
      setTimeout(() => URL.revokeObjectURL(url), 1000);
    }
  });
  $(function () {
    Shiny.addCustomMessageHandler("anca-diagnostics", function (message) {
      // Whitelist fields and bound both entry count and field lengths, even if
      // a malformed message is sent. Never interpret report text as markup.
      if (!message || typeof message !== "object") return;
      snapshot = {
        package_version: text(message.package_version, 40),
        entries: (Array.isArray(message.entries) ? message.entries : []).slice(-100)
          .filter(entry => entry && typeof entry === "object").map(entry => ({
          id: text(entry.id, 80), time: text(entry.time, 30), level: text(entry.level, 10),
          operation: text(entry.operation, 80), summary: text(entry.summary, 400)
        }))
      };
      render();
    });
    $(document).on("shiny:connected", function () { disconnected = false; render(); });
    $(document).on("shiny:disconnected", function () {
      disconnected = true;
      // Native dialog top layer stays above Shiny's disconnected overlay.
      const button = document.getElementById("diagnostics-open");
      if (button) {
        document.body.appendChild(button);
        button.classList.add("diagnostics-offline");
      }
      open();
    });
  });
})();
