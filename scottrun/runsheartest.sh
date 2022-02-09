#! /bin/bash
make sheartest;
for i in 0.02 0.04;
	do echo "-"${i} | ./sheartest >> logfiles/vt${i}.txt &
done
