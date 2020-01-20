# Load-flow-analysis-in-Power-systems
A MATLAB programs for solving the power-flow equations using either of this three methods : Gauss-Seidel (G-S)  , Newton-Raphson (N-R) &amp; Fast Decoupled Load Flow (FLDF).

The key information presented in power flow analysis is the magnitude and phase angle of voltage at each bus and the real and reactive power flowing in each transmission lines. The purpose of any load flow analysis is to compute precise steady-state voltages magnitudes and angles of all buses in the network, the real and reactive power flows into every line and transformer, under the assumption of known generation and load. Power flow analysis is the backbone of power system analysis and design and are necessary for planning, operation, economic scheduling and exchange of power between utilities. 

1) Formation of Y bus matrix is imperative to proceed with this task. One can also get it from here: https://github.com/Beingink
2)Matlab codes for all the three methods were written based on Hadi Saadat’s Power System Analysis and from the flow charts given below.
3)Both the program methods gives Y Bus Matrix ,Line Flow and Line Losses. 
4)However this code is devoid of tap changing transformer effects and charging effects of shunt capacitors. 

Any contribution to the development of this code is very much welcome and appreciated.
