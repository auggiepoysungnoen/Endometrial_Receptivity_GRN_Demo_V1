#!/usr/bin/env bash
# Serve the Nandini_Website folder on localhost for development.
#
# Serves the whole Nandini_Website/ folder (not just a site subfolder) so the
# page can fetch the exported network data under data/ and the reference PNGs
# under figures/ with plain relative URLs, the same way it will on GitHub Pages.
# Browsers block fetch() from file:// URLs, so opening the HTML directly will
# not load the data -- always go through this server.
#
# Usage (from anywhere):
#   ./serve_local.sh            # http://127.0.0.1:8000/
#   ./serve_local.sh 8080       # custom port
#
# Binds to 127.0.0.1 only, so nothing is exposed to the network.
set -euo pipefail

PORT="${1:-8000}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WEBSITE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

if ! command -v python3 >/dev/null 2>&1; then
    echo "ERROR: python3 not found on PATH." >&2
    exit 1
fi

if [[ ! -d "${WEBSITE_ROOT}/data" ]]; then
    echo "ERROR: ${WEBSITE_ROOT}/data not found -- is this script inside Nandini_Website/local_host_development/script/?" >&2
    exit 1
fi

echo "Serving ${WEBSITE_ROOT}"
echo "Site:  http://127.0.0.1:${PORT}/local_host_development/site/   (Ctrl+C to stop)"
exec python3 -m http.server "${PORT}" --bind 127.0.0.1 --directory "${WEBSITE_ROOT}"
