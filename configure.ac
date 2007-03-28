AC_INIT(myconfig, version-0.1)
cd ./inst/Sabre/src
echo "                             Testing for FORTRAN compiler"
AC_PROG_FC
echo "                             Modifying Makefile"
AC_OUTPUT(Makefile)
echo "                             Building Sabre 5.0"
make
echo "                             Installing Sabre 5.0 binary"
cd ../../..
cp ./inst/Sabre/bin/sabre ./exec
echo "                             Tidying up"
cd ./inst/Sabre/src
make clean
rm -rf Makefile
cd ../../..
echo "                             Done"
