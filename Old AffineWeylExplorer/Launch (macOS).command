#!/bin/bash
# Affine Weyl Group Explorer — macOS launcher
# Double-click to start. First run: right-click → Open (Gatekeeper bypass)
cd "$(dirname "$0")"

# Deactivate conda if active (prevents interference with Java)
conda deactivate 2>/dev/null || true

# Check Java
if ! command -v java &>/dev/null; then
    osascript -e 'display alert "Java not found" message "Please install Java 11+ from https://adoptium.net" as critical'
    exit 1
fi

# Fix emoji issue in Java file (safe to run multiple times)
python3 -c "
src = open('AffineWeylExplorer.java').read()
src = src.replace('\\\\U0001F537', '').replace('\\\\U0001F310', '')
open('AffineWeylExplorer.java', 'w').write(src)
" 2>/dev/null || true

echo "Starting Affine Weyl Group Explorer..."
java --source 11 AffineWeylExplorer.java
