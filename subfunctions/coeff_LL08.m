% Coefficients of LL08 = [period  c1  ...  c7  sig]
function C_LL08 = coeff_LL08()
C_LL08 = [
0.00  -2.500  1.205 -1.905 0.51552 0.63255 0.0075 0.275 0.5268
0.01  -2.500  1.205 -1.895 0.51552 0.63255 0.0075 0.275 0.5218
0.02  -2.490  1.200 -1.880 0.51552 0.63255 0.0075 0.275 0.5189
0.03  -2.280  1.155 -1.875 0.51552 0.63255 0.0075 0.275 0.5235
0.04  -2.000  1.100 -1.860 0.51552 0.63255 0.0075 0.275 0.5352
0.05  -1.900  1.090 -1.855 0.51552 0.63255 0.0075 0.275 0.5370
0.06  -1.725  1.065 -1.840 0.51552 0.63255 0.0075 0.275 0.5544
0.09  -1.265  1.020 -1.815 0.51552 0.63255 0.0075 0.275 0.5818
0.10  -1.220  1.000 -1.795 0.51552 0.63255 0.0075 0.275 0.5806
0.12  -1.470  1.040 -1.770 0.51552 0.63255 0.0075 0.275 0.5748
0.15  -1.675  1.045 -1.730 0.51552 0.63255 0.0075 0.275 0.5817
0.17  -1.846  1.065 -1.710 0.51552 0.63255 0.0075 0.275 0.5906
0.20  -2.170  1.085 -1.675 0.51552 0.63255 0.0075 0.275 0.6059
0.24  -2.585  1.105 -1.630 0.51552 0.63255 0.0075 0.275 0.6315
0.30  -3.615  1.215 -1.570 0.51552 0.63255 0.0075 0.275 0.6656
0.36  -4.160  1.255 -1.535 0.51552 0.63255 0.0075 0.275 0.7010
0.40  -4.595  1.285 -1.500 0.51552 0.63255 0.0075 0.275 0.7105
0.46  -5.020  1.325 -1.495 0.51552 0.63255 0.0075 0.275 0.7148
0.50  -5.470  1.365 -1.465 0.51552 0.63255 0.0075 0.275 0.7145
0.60  -6.095  1.420 -1.455 0.51552 0.63255 0.0075 0.275 0.7177
0.75  -6.675  1.465 -1.450 0.51552 0.63255 0.0075 0.275 0.7689
0.85  -7.320  1.545 -1.450 0.51552 0.63255 0.0075 0.275 0.7787
1.00  -8.000  1.620 -1.450 0.51552 0.63255 0.0075 0.275 0.7983
1.50  -9.240  1.705 -1.440 0.51552 0.63255 0.0075 0.275 0.8411
2.00 -10.200  1.770 -1.430 0.51552 0.63255 0.0075 0.275 0.8766
3.00 -11.470  1.830 -1.370 0.51552 0.63255 0.0075 0.275 0.8590
4.00 -12.550  1.845 -1.260 0.51552 0.63255 0.0075 0.275 0.8055
5.00 -13.390  1.805 -1.135 0.51552 0.63255 0.0075 0.275 0.7654
];
%{
C_LL08_stiffsoil = [
0.00  -0.900  1.000 -1.900 0.99178 0.52632 0.0040 0.310 0.6277
0.01  -2.200  1.085 -1.750 0.99178 0.52632 0.0040 0.310 0.5800
0.02  -2.290  1.085 -1.730 0.99178 0.52632 0.0040 0.310 0.5730
0.03  -2.340  1.095 -1.720 0.99178 0.52632 0.0040 0.310 0.5774
0.04  -2.215  1.090 -1.730 0.99178 0.52632 0.0040 0.310 0.5808
0.05  -1.895  1.055 -1.755 0.99178 0.52632 0.0040 0.310 0.5937
0.06  -1.110  1.010 -1.835 0.99178 0.52632 0.0040 0.310 0.6123
0.09  -0.210  0.945 -1.890 0.99178 0.52632 0.0040 0.310 0.6481
0.10  -0.055  0.920 -1.880 0.99178 0.52632 0.0040 0.310 0.6535
0.12   0.055  0.935 -1.895 0.99178 0.52632 0.0040 0.310 0.6585
0.15  -0.040  0.955 -1.880 0.99178 0.52632 0.0040 0.310 0.6595
0.17  -0.340  1.020 -1.885 0.99178 0.52632 0.0040 0.310 0.6680
0.20  -0.800  1.045 -1.820 0.99178 0.52632 0.0040 0.310 0.6565
0.24  -1.575  1.120 -1.755 0.99178 0.52632 0.0040 0.310 0.6465
0.30  -3.010  1.315 -1.695 0.99178 0.52632 0.0040 0.310 0.6661
0.36  -3.680  1.380 -1.660 0.99178 0.52632 0.0040 0.310 0.6876
0.40  -4.250  1.415 -1.600 0.99178 0.52632 0.0040 0.310 0.7002
0.46  -4.720  1.430 -1.545 0.99178 0.52632 0.0040 0.310 0.7092
0.50  -5.220  1.455 -1.490 0.99178 0.52632 0.0040 0.310 0.7122
0.60  -5.700  1.470 -1.445 0.99178 0.52632 0.0040 0.310 0.7280
0.75  -6.450  1.500 -1.380 0.99178 0.52632 0.0040 0.310 0.7752
0.85  -7.250  1.565 -1.325 0.99178 0.52632 0.0040 0.310 0.7931
1.00  -8.150  1.605 -1.235 0.99178 0.52632 0.0040 0.310 0.8158
1.50 -10.300  1.800 -1.165 0.99178 0.52632 0.0040 0.310 0.8356
2.00 -11.620  1.860 -1.070 0.99178 0.52632 0.0040 0.310 0.8474
3.00 -12.630  1.890 -1.060 0.99178 0.52632 0.0040 0.310 0.8367
4.00 -13.420  1.870 -0.990 0.99178 0.52632 0.0040 0.310 0.7937
5.00 -13.750  1.835 -0.975 0.99178 0.52632 0.0040 0.310 0.7468
];
%}