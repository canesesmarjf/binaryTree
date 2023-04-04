if [ $(uname) = "Linux" ]; then
  export LD_LIBRARY_PATH=./arma_libs/lib/:$LD_LIBRARY_PATH
fi

FORGE_PATH=/home/jfcm/arm/forge/22.0.1/bin
export PATH=${FORGE_PATH}:$PATH

ddt ./bin/BinarySearch.exe
