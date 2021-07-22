clear all
clc

Dimension = 3;          % dimensi diganti sesuai dengan jumlah variabel yang dioptimasi
UB = [0.625 1300 40];           % Upper Bounds diganti sesuai dengan constraint fungsi objektif
LB = [0.2 1071 31];            % Lower Bounds diganti sesuai dengan constraint fungsi objektif
save ('propco2egr.mat')

clear all;
close all;
clc;