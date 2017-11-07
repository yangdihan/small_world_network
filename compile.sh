#!/bin/bash
make clean
mpic++ -std=c++11 main.cpp network.cpp sac_network.cpp mpi_network.cpp crack.cpp helper_funcs.cpp -o $1
