if [ $(uname) = "Linux" ]; then
  export LD_LIBRARY_PATH=./arma_libs/lib/:$LD_LIBRARY_PATH
fi
./bin/BinarySearch.exe
