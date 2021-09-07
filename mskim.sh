export CXX="mpicxx.mpich"
export CC="mpicc.mpich"
export CPP="cpp"

export CXXFLAGS="-O3 -Wall"
export SU2_INSTALLPATH="/home/ubuntu/mskim/SU2"
#export LAPACK_LIB="/home/mskim/lapack-3.9.1"
#export LAPACK_INCLUDE="/home/mskim/lapack-3.9.1"
export LAPACK_LIB="/usr/lib"
export LAPACK_INCLUDE="/usr/include"
#export LAPACK_INCLUDE="/home/mskim/lapack-3.9.1/LAPACKE/include"

./preconfigure.py --enable-autodiff --enable-mpi --enable-PY_WRAPPER --prefix=/home/ubuntu/mskim/SU2

#./preconfigure.py --enable-autodiff --enable-mpi --enable-PY_WRAPPER --with-LAPACK-lib=$LAPACK_LIB --with-LAPACK-include=$LAPACK_INCLUDE --prefix=/home/mskim/SU2
#./preconfigure.py --enable-autodiff --enable-mpi --enable-PY_WRAPPER --prefix=/home/mskim/SU2

make -j 6 > log_make
make install -j 6 > log_make_install

