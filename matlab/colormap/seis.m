function h = seis(m)
cmap=[
175.98,0,0;
187.93,0,0;
199.88,0,0;
211.84,0,0;
223.79,0,0;
235.74,0,0;
247.7,0,0;
255,4.6484,0;
255,16.602,0;
255,28.555,0;
255,40.508,0;
255,52.461,0;
255,64.414,0;
255,76.367,0;
255,88.32,0;
255,100.27,0;
255,112.23,0;
255,124.18,0;
255,136.13,0;
255,148.09,0;
255,160.04,0;
255,171.99,0;
255,183.95,0;
255,195.9,0;
255,207.85,0;
255,219.8,0;
255,231.76,0;
255,243.71,0;
255,255,0;
255,255,0;
255,255,0;
255,255,0;
255,255,0;
255,255,0;
255,255,0;
255,255,0;
233.09,255,3.9844;
209.88,255,8.2031;
186.68,255,12.422;
163.48,255,16.641;
140.27,255,20.859;
117.07,255,25.078;
93.867,255,29.297;
79.453,253.24,39.375;
66.797,251.13,50.625;
54.141,249.02,61.875;
41.484,246.91,73.125;
28.828,244.8,84.375;
16.172,242.7,95.625;
3.5156,240.59,106.87;
0,223.75,124.73;
0,201.25,145.12;
0,178.75,165.51;
0,156.25,185.9;
0,133.75,206.29;
0,111.25,226.68;
0,88.75,247.07;
0,73.125,250.7;
0,61.875,243.67;
0,50.625,236.64;
0,39.375,229.61;
0,28.125,222.58;
0,16.875,215.55;
0,5.625,208.52;
]/255;
if nargin < 1, m = 64; end
h = cmap;
