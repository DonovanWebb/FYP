import matplotlib.pyplot as plt
import numpy as np

###  exp data
R_L_e = 1000*np.array([
    0.115,
    0.515,
    0.765,
    1.015,
    2.015,
    3.015,
    4.015,
    5.015,
    6.015,
    7.015,
    8.015,
    9.015,
    10.015,
    100.015,
    200.015,
    300.015,
    400.015,
    500.015,
    750.015,
    1000.015,
    ])
P_L_e = np.array([
    0.00001476034783,
    0.00001005483806,
    0.00001373366013,
    0.00001753257143,
    0.00002521347891,
    0.00002951339635,
    0.00003161885181,
    0.00003269072981,
    0.00003264133167,
    0.00003232593585,
    0.00003169257642,
    0.00003076068331,
    0.00003049093959,
    0.000006484108284,
    0.000003428490913,
    0.000002275791844,
    0.000001708104846,
    0.000001383740968,
    0.0000009472882676,
    0.0000007097956031,
    ])


R_L = np.linspace(0,10**6, 10**6)
r = 5550
V = 0.85

P_L = (V-V*r/(r+R_L))**2/R_L
P_LI = (V/(r+R_L))**2*R_L
R_T = r+R_L
b = 1
n = (b/(R_T**2 + b*R_T))*R_L

w = 15000
E = 0.165
R_1 = 100
L_21 = 0.0471

P_L__ = (w**2*L_21**2*E**2)/(2*R_L)*1/(R_1 + (w**2*L_21**2/R_L))**2

plt.figure(1)
# plt.plot((R_L), (P_LI))
plt.plot((R_L_e), (P_L_e))
# plt.plot((R_L), (n))
plt.plot((R_L), (P_L__))

plt.figure(2)
plt.plot((R_L)[0:10**5], (P_L)[0:10**5])
plt.plot((R_L)[0:10**5], (P_L__)[0:10**5])
plt.plot((R_L_e)[1:13], (P_L_e)[1:13])

plt.show()
