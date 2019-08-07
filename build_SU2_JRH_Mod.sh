#!bash
echo "JRH - Cleaning..."

echo "Moving ~/SU2/bin directory to bin_backup"
#rm -rf bin_backup
#mv bin bin_backup

rm -rf bin

make clean
make distclean

echo "********JRH - cleanining complete...Configuring...**********"

#./configure --prefix=$SU2_HOME --enable-PY_WRAPPER --enable-metis --enable-codi-reverse --enable-codi-forward --enable-mpi

#Do this for full build
#CXX CC AND CPP ADDED 2232019 IN ATTEMPT TO FIX ISSUE WITH PARMETIS MPI IN SU2_CFD_AD
export CXX="mpicxx.mpich"
export CC="mpicc.mpich"
export CPP="cpp"

export CXXFLAGS="-O3 -Wall"
export SU2_INSTALLPATH="/home/jon/SU2"
export LAPACK_LIB="/usr/lib"
export LAPACK_INCLUDE="/usr/include"
./preconfigure.py --enable-autodiff --enable-mpi --enable-PY_WRAPPER --with-LAPACK-lib=$LAPACK_LIB --with-LAPACK-include=$LAPACK_INCLUDE --prefix=$SU2_INSTALLPATH

#./preconfigure.py --enable-autodiff --enable-mpi --enable-PY_WRAPPER --enable-cgns --prefix=$SU2_INSTALLPATH

#Use #--enable-PY_WRAPPER for python wrap
#make -j 4

echo "********JRH - Configuration Complete...Installing*************"

#Install - sudo permissions required
make install -j 4

echo "******** JRH - Complete***********"

