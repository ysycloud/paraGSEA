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
quick_search_serial