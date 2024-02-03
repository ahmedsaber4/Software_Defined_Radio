clear all ;
clc;
filters = ["Blackman", "Chebychev", "Equiripple", "Least Squares"];

orders = [199, 200, 273, 200];
passbandRipples = [0.02, 0.02, 1, 0.201];
stopbandAttenuations = [-79.4, -112.3, -39.39, -29.4973];
transitionBW = [87, 87, 24.8, 30.23];

optimum_filter = getOptimumFilter(filters, orders, passbandRipples, stopbandAttenuations, transitionBW);