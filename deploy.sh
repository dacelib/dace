#!/bin/sh
cd ${TRAVIS_BUILD_DIR}/build
  git clone -b pkg-latest https://abgandar:${GITHUB_TOKEN}@github.com/abgandar/dace.git pkg-latest
cd pkg-latest
cp ${TRAVIS_BUILD_DIR}/build/packages/* .
git add *
git commit -a -m "Automatic update of build artifacts for commit ${TRAVIS_COMMIT}"
git push
cd ..
rm -rf pkg-latest
