#!/bin/sh
cd ${TRAVIS_BUILD_DIR}/build
git clone -b pkg-latest https://${GITHUB_USER}:${GITHUB_TOKEN}@github.com/${GITHUB_USER}/dace.git pkg-latest >/dev/null 2>&1
cd pkg-latest
cp ${TRAVIS_BUILD_DIR}/build/packages/* .
git add *
git commit -a -m "Automatic update of build artifacts for commit ${TRAVIS_COMMIT}"
git push
cd ..
rm -rf pkg-latest
