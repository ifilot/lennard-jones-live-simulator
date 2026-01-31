#!/usr/bin/env bash
set -euo pipefail

# ============================
# Configuration
# ============================
APP_NAME="ljsim"
APP_EXE="ljsim.exe"
BUILD_DIR="build"
DIST_ROOT="dist"
DIST_DIR="${DIST_ROOT}/${APP_NAME}"
BUILD_TYPE="Release"

# ============================
# Sanity checks
# ============================
command -v cmake >/dev/null || { echo "cmake not found"; exit 1; }
command -v makensis >/dev/null || { echo "makensis not found"; exit 1; }

if [[ -z "${MINGW_PREFIX:-}" ]]; then
  echo "MINGW_PREFIX not set (are you in a MinGW shell?)"
  exit 1
fi

echo "[INFO] Using MINGW_PREFIX=${MINGW_PREFIX}"

# ============================
# Clean
# ============================
echo "[INFO] Cleaning previous output"
rm -rf "${DIST_ROOT}"

# ============================
# Configure
# ============================
echo "[INFO] Configuring (${BUILD_TYPE})"
cmake -S . -B "${BUILD_DIR}" \
  -G Ninja \
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}"

# ============================
# Build
# ============================
echo "[INFO] Building"
cmake --build "${BUILD_DIR}"

# ============================
# Stage files
# ============================
echo "[INFO] Staging files"
mkdir -p "${DIST_DIR}"
cp "${BUILD_DIR}/${APP_EXE}" "${DIST_DIR}/"

# ============================
# Locate windeployqt
# ============================
: "${MINGW_PREFIX:=/mingw64}"

if command -v windeployqt >/dev/null 2>&1; then
  WINDEPLOYQT="$(command -v windeployqt)"
elif [[ -x "$MINGW_PREFIX/bin/windeployqt-qt5.exe" ]]; then
  WINDEPLOYQT="$MINGW_PREFIX/bin/windeployqt-qt5.exe"
elif [[ -x "$MINGW_PREFIX/bin/windeployqt-qt6.exe" ]]; then
  WINDEPLOYQT="$MINGW_PREFIX/bin/windeployqt-qt6.exe"
else
  echo "ERROR: windeployqt not found"
  exit 1
fi

: "${MINGW_PREFIX:=/mingw64}"

QMAKE_QT5="$MINGW_PREFIX/bin/qmake-qt5.exe"
QMAKE_SHIM="$MINGW_PREFIX/bin/qmake.exe"

if [[ ! -x "$QMAKE_QT5" ]]; then
  echo "qmake-qt5.exe not found at $QMAKE_QT5"
  exit 1
fi

# Provide qmake.exe for windeployqt (MSYS2 requirement)
if [[ ! -x "$QMAKE_SHIM" ]]; then
  echo "[INFO] Creating qmake.exe shim"
  ln -s qmake-qt5.exe "$QMAKE_SHIM" 2>/dev/null || cp "$QMAKE_QT5" "$QMAKE_SHIM"
fi

# Canonical qmake path from here on
QMAKE="$QMAKE_SHIM"

"$WINDEPLOYQT" \
  --no-angle \
  --no-opengl-sw \
  "${DIST_DIR}/${APP_EXE}"

echo "[INFO] Bundling MinGW runtime DLLs"

ldd "${DIST_DIR}/${APP_EXE}" | awk '{print $3}' | while read -r path; do
  case "$path" in
    /mingw64/bin/*.dll)
      dll="$(basename "$path")"
      echo "  + $dll"
      cp -n "$path" "${DIST_DIR}/"
      ;;
  esac
done

# ============================
# Build installer
# ============================
echo "[INFO] Building NSIS installer"
makensis installer.nsi

echo
echo "[SUCCESS] Windows installer created"