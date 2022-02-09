#! /bin/bash
make bulktest;
for i in 0.005 0.01 0.02 0.03 0.04 0.05;
	do echo "-"${i} | ./bulktest >> logfiles/vt${i}.txt &
done
