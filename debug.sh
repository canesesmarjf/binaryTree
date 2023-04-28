if [ $(uname) = "Linux" ]; then
  # export LD_LIBRARY_PATH=./arma_libs/lib/:$LD_LIBRARY_PATH
  export LD_LIBRARY_PATH=/home/jfcm/Repos/binaryTree/arma_libs/lib:$LD_LIBRARY_PATH
fi

FORGE_PATH=/home/jfcm/arm/forge/22.0.1/bin
export PATH=${FORGE_PATH}:$PATH

if [ "$1" = "1" ]; then
  ddt ./bin/main_1.exe
elif [ "$1" = "2" ]; then
  ddt ./bin/main_2.exe
elif [ "$1" = "3" ]; then
  ddt ./bin/main_3.exe
elif [ "$1" = "4" ]; then
  ddt ./bin/main_4.exe
else
  echo "Invalid argument. Usage: $0 [1|2|3|4]"
  exit 1
fi
