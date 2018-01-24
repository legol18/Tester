touch 1.dat
echo 0.0 > 1.dat

touch 2.dat
echo 0.8 > 2.dat

touch 3.dat
echo -0.8 > 3.dat

ifort random.f90 test_up.f90 -o tt -r8
./tt < 1.dat
cp fort.11 zero.dat
cp fort.21 zero_sp_rk_p01.dat
rm fort.11
rm fort.21

./tt < 2.dat
cp fort.11 p8.dat
cp fort.21 p8_sp_rk_p01.dat
rm fort.11
rm fort.21

./tt < 3.dat
cp fort.11 mp8.dat
cp fort.21 mp8_sp_rk_p01.dat
rm fort.11
rm 1.dat
rm 2.dat
rm 3.dat

gnuplot plotting.gpl
rm *.mod
