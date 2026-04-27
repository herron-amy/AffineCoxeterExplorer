#!/bin/bash
# Affine Weyl Group Explorer — Linux launcher
# Make executable first: chmod +x "Launch (Linux).sh"
cd "$(dirname "$0")"

if ! command -v java &>/dev/null; then
    echo "Java not found. Install with:"
    echo "  sudo apt install default-jre   # Debian/Ubuntu"
    echo "  sudo dnf install java-latest-openjdk  # Fedora"
    exit 1
fi

echo "Starting Affine Weyl Group Explorer..."
java --source 11 AffineWeylExplorer.java
