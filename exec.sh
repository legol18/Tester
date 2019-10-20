touch 1.dat
echo 0.0 > 1.dat

touch 2.dat
echo 0.8 > 2.dat

touch 3.dat
echo -0.8 > 3.dat

#### for testing
# ifort random.f90 test_up.f90 -O2  -stand f03 -assume realloc_lhs  -check all  -traceback  -warn all  -fstack-protector  -assume protect_parens  -implicitnone -o tt

gfortran random.f90 ref_prog.f90 -o tt -r8 
./tt < 1.dat
cp fort.11 zero_w.dat

./tt < 2.dat
cp fort.11 p8_w.dat

./tt < 3.dat
cp fort.11 mpzero6_w.dat

rm 1.dat
rm 2.dat
rm 3.dat
rm fort.*

gnuplot plotting.gpl

rm *.mod
rm tt
