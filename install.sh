#!/bin/bash

#git clone
git clone https://github.com/ysycloud/paraGSEA.git
cd paraGSEA
#make
make all
#install
make install
#clean
make clean
#Test
quick_search_serial -i data/data_for_test.txt -n 8 -s data/data_for_test_cidnum.txt -r data/Reference