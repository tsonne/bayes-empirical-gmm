% Coefficients of Am05 = [period  a1  ...  a10  s1A s1B s2A s2B]
function C_Am05 = coeff_Am05()
C_Am05 = [
0.000  2.522 -0.142 -3.184  0.314 7.6  0.137  0.050 -0.084  0.062 -0.044  0.665  0.065  0.222  0.022
0.050  3.247 -0.225 -3.525  0.359 7.4  0.098  0.005 -0.096  0.078 -0.048  0.708  0.069  0.249  0.024
0.055  3.125 -0.206 -3.418  0.345 7.1  0.085  0.004 -0.096  0.072 -0.050  0.672  0.063  0.235  0.022
0.060  3.202 -0.212 -3.444  0.347 7.4  0.079  0.002 -0.103  0.073 -0.047  0.687  0.065  0.237  0.023
0.065  3.442 -0.242 -3.571  0.365 7.7  0.069  0.001 -0.104  0.076 -0.035  0.693  0.067  0.241  0.023
0.070  3.504 -0.249 -3.576  0.367 7.9  0.064 -0.002 -0.114  0.068 -0.043  0.647  0.059  0.225  0.021
0.075  3.472 -0.240 -3.521  0.358 8.0  0.064 -0.003 -0.121  0.063 -0.046  0.674  0.063  0.227  0.021
0.080  3.526 -0.248 -3.520  0.358 8.1  0.069 -0.002 -0.116  0.074 -0.040  0.756  0.076  0.252  0.025
0.085  3.320 -0.215 -3.381  0.336 8.0  0.067  0.010 -0.116  0.075 -0.039  0.750  0.076  0.258  0.026
0.090  3.309 -0.211 -3.353  0.332 7.9  0.064  0.014 -0.119  0.065 -0.048  0.727  0.072  0.249  0.025
0.095  3.479 -0.240 -3.420  0.345 7.8  0.062  0.014 -0.107  0.073 -0.051  0.772  0.079  0.262  0.027
0.100  3.596 -0.258 -3.511  0.360 7.9  0.065  0.025 -0.095  0.076 -0.047  0.747  0.075  0.249  0.025
0.110  3.453 -0.239 -3.398  0.345 7.9  0.077  0.041 -0.082  0.072 -0.052  0.810  0.084  0.256  0.027
0.120  3.330 -0.214 -3.300  0.329 8.0  0.070  0.045 -0.081  0.065 -0.046  0.753  0.075  0.240  0.024
0.130  3.249 -0.195 -3.254  0.321 8.2  0.069  0.043 -0.084  0.056 -0.059  0.712  0.068  0.236  0.023
0.140  2.993 -0.154 -3.088  0.297 8.2  0.065  0.042 -0.074  0.053 -0.067  0.650  0.059  0.218  0.020
0.150  2.725 -0.111 -2.909  0.270 8.3  0.067  0.044 -0.074  0.067 -0.060  0.634  0.057  0.223  0.020
0.160  2.738 -0.120 -2.912  0.274 8.2  0.085  0.049 -0.069  0.090 -0.061  0.734  0.072  0.251  0.025
0.170  2.692 -0.114 -2.907  0.275 8.2  0.091  0.053 -0.059  0.087 -0.055  0.760  0.077  0.257  0.026
0.180  2.665 -0.110 -2.907  0.276 8.1  0.098  0.049 -0.057  0.087 -0.054  0.736  0.073  0.251  0.025
0.190  2.713 -0.118 -2.989  0.288 8.1  0.112  0.059 -0.050  0.090 -0.054  0.752  0.076  0.250  0.025
0.200  2.632 -0.109 -2.990  0.289 8.1  0.124  0.070 -0.033  0.090 -0.039  0.784  0.080  0.251  0.026
0.220  2.483 -0.088 -2.941  0.281 7.9  0.136  0.078 -0.033  0.086 -0.024  0.778  0.079  0.244  0.025
0.240  2.212 -0.051 -2.823  0.265 7.6  0.156  0.087 -0.037  0.090 -0.020  0.770  0.077  0.235  0.024
0.260  2.058 -0.036 -2.787  0.263 7.3  0.179  0.077 -0.024  0.120  0.010  0.917  0.101  0.278  0.030
0.280  1.896 -0.010 -2.732  0.251 7.5  0.193  0.074 -0.023  0.112  0.027  0.947  0.104  0.285  0.031
0.300  1.739  0.009 -2.667  0.244 7.1  0.192  0.069 -0.034  0.104  0.012  0.890  0.095  0.267  0.028
0.320  1.728  0.001 -2.688  0.251 7.1  0.207  0.073 -0.021  0.118  0.008  0.917  0.098  0.273  0.029
0.340  1.598  0.020 -2.667  0.246 7.2  0.216  0.078 -0.010  0.118  0.005  0.896  0.095  0.261  0.028
0.360  1.477  0.034 -2.641  0.244 6.9  0.230  0.091 -0.013  0.107 -0.011  0.846  0.087  0.254  0.026
0.380  1.236  0.071 -2.534  0.227 6.7  0.247  0.100 -0.010  0.106 -0.018  0.803  0.080  0.250  0.025
0.400  1.070  0.091 -2.474  0.219 6.3  0.256  0.097 -0.013  0.115 -0.020  0.793  0.078  0.244  0.024
0.420  0.998  0.096 -2.469  0.220 5.9  0.259  0.100 -0.021  0.116 -0.024  0.757  0.072  0.233  0.022
0.440  1.045  0.085 -2.540  0.231 6.3  0.269  0.114 -0.016  0.114 -0.028  0.787  0.077  0.241  0.024
0.460  0.980  0.093 -2.564  0.234 6.3  0.278  0.122 -0.011  0.108 -0.029  0.766  0.074  0.238  0.023
0.480  0.874  0.103 -2.530  0.231 6.2  0.286  0.130  0.001  0.118 -0.024  0.778  0.076  0.240  0.023
0.500  0.624  0.139 -2.410  0.212 6.1  0.289  0.133  0.004  0.126 -0.026  0.798  0.079  0.246  0.024
0.550  0.377  0.174 -2.317  0.196 6.1  0.293  0.137 -0.004  0.118 -0.035  0.841  0.085  0.268  0.027
0.600  0.359  0.158 -2.343  0.206 5.4  0.311  0.136  0.008  0.118 -0.028  0.919  0.099  0.308  0.033
0.650  0.130  0.182 -2.294  0.202 5.0  0.318  0.149  0.005  0.107 -0.031  0.867  0.090  0.301  0.031
0.700 -0.014  0.198 -2.305  0.205 4.8  0.327  0.154 -0.011  0.105 -0.032  0.803  0.080  0.298  0.030
0.750 -0.307  0.236 -2.201  0.191 4.7  0.318  0.148 -0.001  0.114 -0.032  0.774  0.076  0.278  0.027
0.800 -0.567  0.279 -2.083  0.170 5.2  0.332  0.178 -0.003  0.083 -0.062  0.661  0.059  0.240  0.021
0.850 -0.519  0.262 -2.177  0.186 4.9  0.341  0.183  0.005  0.085 -0.070  0.694  0.064  0.253  0.023
0.900 -0.485  0.249 -2.246  0.199 4.5  0.354  0.191 -0.003  0.072 -0.082  0.714  0.067  0.263  0.025
0.950 -1.133  0.369 -1.957  0.143 5.5  0.353  0.204 -0.025  0.024 -0.109  0.309  0.000  0.121  0.000
1.000 -1.359  0.403 -1.848  0.124 6.0  0.357  0.211 -0.013  0.024 -0.101  0.305  0.000  0.120  0.000
1.100 -1.675  0.437 -1.711  0.108 5.5  0.373  0.213 -0.029 -0.007 -0.108  0.306  0.000  0.118  0.000
1.200 -1.982  0.477 -1.636  0.095 5.4  0.389  0.226 -0.014 -0.017 -0.095  0.297  0.000  0.120  0.000
1.300 -2.226  0.511 -1.605  0.089 5.5  0.395  0.215 -0.004 -0.025 -0.085  0.296  0.000  0.119  0.000
1.400 -2.419  0.533 -1.541  0.080 6.0  0.408  0.237  0.028 -0.040 -0.091  0.290  0.000  0.115  0.000
1.500 -2.639  0.550 -1.443  0.074 4.9  0.405  0.229  0.020 -0.053 -0.133  0.292  0.000  0.111  0.000
1.600 -2.900  0.587 -1.351  0.060 5.2  0.387  0.216  0.019 -0.056 -0.131  0.296  0.000  0.114  0.000
1.700 -2.695  0.564 -1.564  0.086 6.5  0.380  0.212  0.001 -0.081 -0.141  0.302  0.000  0.117  0.000
1.800 -3.209  0.630 -1.410  0.069 5.4  0.391  0.174  0.012 -0.035 -0.154  0.291  0.000  0.128  0.000
1.900 -3.313  0.647 -1.424  0.067 5.9  0.386  0.175  0.030 -0.033 -0.145  0.290  0.000  0.133  0.000
2.000 -3.063  0.586 -1.372  0.070 4.2  0.421  0.177  0.008 -0.019 -0.174  0.282  0.000  0.134  0.000
2.100 -3.043  0.578 -1.435  0.080 4.3  0.404  0.171  0.002 -0.026 -0.164  0.281  0.000  0.134  0.000
2.200 -3.068  0.575 -1.448  0.083 4.2  0.394  0.160 -0.007 -0.034 -0.169  0.283  0.000  0.136  0.000
2.300 -3.996  0.740 -0.829 -0.025 5.1  0.349  0.135 -0.010 -0.031 -0.125  0.282  0.000  0.137  0.000
2.400 -4.108  0.758 -0.755 -0.038 5.3  0.338  0.119 -0.024 -0.050 -0.147  0.284  0.000  0.137  0.000
2.500 -4.203  0.768 -0.714 -0.044 5.1  0.325  0.103 -0.026 -0.063 -0.155  0.285  0.000  0.137  0.000
];