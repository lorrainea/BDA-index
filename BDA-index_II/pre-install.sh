#! /bin/sh

tar -xvf sdsl-lite.tar.gz
cd sdsl-lite
./install.sh "$(pwd)"/libsdsl
mv libsdsl/ ..

cd psascan
make

cd ../sparsePhi/src
make

