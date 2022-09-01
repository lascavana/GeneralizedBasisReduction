rm -r build
mkdir build
cmake -H. -Bbuild
cmake --build build

./build/ls_reduction example.lp
