if [ $(uname) = "Linux" ]; then
  export LD_LIBRARY_PATH=./arma_libs/lib/:$LD_LIBRARY_PATH
fi

if [ "$1" = "1" ]; then
  ./bin/main_1.exe
elif [ "$1" = "2" ]; then
  ./bin/main_2.exe
elif [ "$1" = "3" ]; then
  ./bin/main_3.exe
else
  echo "Invalid argument. Usage: $0 [1|2|3]"
  exit 1
fi
