#!/bin/bash
mkdir -p build-test
cd build-test
#export PKG_CONFIG_PATH="/usr/lib/x86_64-linux-gnu/pkgconfig:/usr/share/pkgconfig"
cmake ../CMakeLists.txt
RET="$?"
if [ "$RET" -eq 0 ]
then
    make
    RET="$?"
else
    echo "cmake failed"
fi
cd ..
#return "$RET"
