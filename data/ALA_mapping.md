
Symbolic Resolution of Labeled Measurements
===========================================

Contents
========

* [SymNmrMs](#symnmrms)
	* [Problem setup](#problem-setup)
	* [Solution](#solution)
	* [System Status](#system-status)
	* [Redundant measurements](#redundant-measurements)
	* [Measurable elementary combinations](#measurable-elementary-combinations)
	* [Least Squares](#least-squares)
  
<a name="symnmrms"></a>
# SymNmrMs
  
<a name="problem-setup"></a>
## Problem setup
  
<a name="explicit-parameters"></a>
### Explicit Parameters


```
mm='./data/ALA_mapping.tsv' colsel='[1, 2, 3, 4, 5, 6, 7, 8]' tim=None data='data.frame 21x3' w=True s=None rd=None
```  
<a name="available-measurements"></a>
### Available measurements
  

||D1|D2|D3|D4|D5|D6|D7|D8|
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
|000||||e|t||o|g|
|001||k||e|u||p|h|
|100|||r|e|t||o|h|
|101||k|r|e|u||p|i|
|010|a|||f|t|m|p|h|
|011|b|l||f|u|m|q|i|
|110|c||s|f|t|n|p|i|
|111|d|l|s|f|u|n|q|j|
  
<a name="equation-system"></a>
### Equation system
  
i000 + i001 + i100 + i101 + i010 + i011 + i110 + i111 = 1  
a·(i010 + i011 + i110 + i111) - i010 = 0  
b·(i010 + i011 + i110 + i111) - i011 = 0  
c·(i010 + i011 + i110 + i111) - i110 = 0  
d·(i010 + i011 + i110 + i111) - i111 = 0  
-i001 - i101 + k·(i001 + i011 + i101 + i111) = 0  
-i011 - i111 + l·(i001 + i011 + i101 + i111) = 0  
-i100 - i101 + r·(i100 + i101 + i110 + i111) = 0  
-i110 - i111 + s·(i100 + i101 + i110 + i111) = 0  
e - i000 - i001 - i100 - i101 = 0  
f - i010 - i011 - i110 - i111 = 0  
-i000 - i010 - i100 - i110 + t = 0  
-i001 - i011 - i101 - i111 + u = 0  
-i010 - i011 + m·(i010 + i011 + i110 + i111) = 0  
-i110 - i111 + n·(i010 + i011 + i110 + i111) = 0  
-i000 - i100 + o = 0  
-i001 - i010 - i101 - i110 + p = 0  
-i011 - i111 + q = 0  
g - i000 = 0  
h - i001 - i010 - i100 = 0  
i - i011 - i101 - i110 = 0  
-i111 + j = 0  
<a name="solution"></a>
## Solution
  
<a name="isotopomers"></a>
### Isotopomers

|variable|formula|
| :---: | :---: |
|000|g|
|001|-a·f + g + h - o|
|100|-g + o|
|101|a·f + e - g - h|
|010|a·f|
|011|-j + q|
|110|-a·f - o + t|
|111|j|
  
<a name="emu"></a>
### EMU

|variable|Exx|xEx|xxE|EEx|ExE|xEE|EEE|
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
|M+0|2·g + h - j - o + q|e|t|-a·f + 2·g + h - o|a·f + g|o|g|
|M+1|-2·g - h + j + o - q + 1|f|u|2·a·f - 2·g - h - j + 2·o + u|f·(1 - 2·a) + h - j|p|h|
|M+2||||-a·f + j - o + t|a·f + e - g - h + j|q|i|
|M+3|||||||j|
  
<a name="cumomers"></a>
### Cumomers

|variable|formula|
| :---: | :---: |
|1xx|-2·g - h + j + o - q + 1|
|x1x|f|
|xx1|u|
|11x|-a·f + j - o + t|
|1x1|a·f + e - g - h + j|
|x11|q|
|111|j|
  
<a name="system-status"></a>
## System Status


The system is *over-determined*.

The reason is that following measurements were not used: *b*, *d*, *e*, *f*, *k*, *l*, *m*, *n*, *r*, *s*.

Number of redundant measurements: 10.  
<a name="redundant-measurements"></a>
## Redundant measurements

|Measurement|Formula|Methods|
| :---: | :---: | :--- |
|b|-(j - o + t)/(-o + q + t) + 1|D1: *D5, D7, D8*|
|d|j/(-o + q + t)|D1: *D5, D7, D8*|
|e|o - q + u|D4: *D5, D7*|
|f|-o + q + t|D4: *D5, D7*|
|k|-q/u + 1|D2: *D5, D7*|
|l|q/u|D2: *D5, D7*|
|m|a - (j - o + t)/(-o + q + t) + 1|D6: *D1, D5, D7, D8*|
|n|-a + (j - o + t)/(-o + q + t)|D6: *D1, D5, D7, D8*|
|r|-a - (a·(-2·g - h + j + 2·t + u) - j + o - t)/(2·g + h - j - o + q - 1) + 1|D3: *D1, D5, D7, D8*|
|s|a + (a·(-2·g - h + j + 2·t + u) - j + o - t)/(2·g + h - j - o + q - 1)|D3: *D1, D5, D7, D8*|
  
<a name="measurable-elementary-combinations"></a>
## Measurable elementary combinations

|N°|Combination|Accumulated|Formula|Methods|
| :---: | :---: | :---: | :---: | :---: |
|1|000|000|g|D8|
|2|001|001|-a·f + g + h - o|D1, D4, D7, D8|
|3|100|100|-g + o|D7, D8|
|4|101|101|a·f + e - g - h|D1, D4, D8|
|5|010|010|a·f|D1, D4|
|6|011|011|-j + q|D7, D8|
|7|110|110|-a·f - o + t|D1, D4, D5, D7|
|8|111|111|j|D8|
  
<a name="least-squares"></a>
## Least Squares
  
<a name="problem-summary"></a>
### Problem summary

|Item|Value|Comment|
| :---: | :---: | :---: |
|number of measurements|21|m|
|system rank (number of statistically defined parameters)|7|p|
|degree of freedom|14|dof=m - p|
|χ²(14)|106.679|Σ((estimated-measured)/sd)²<br/>95% confidence interval is [0; 24]|
|χ²(14)/14|7.62|Reduced χ², χ²(dof)/dof: if ≫ 1 then poor fitting, ~ 1 means 'good' fitting and ≪ 1 means over-fitting or overestimating sd|
|p-value of one-tail χ² test, i.e.<br/>P(χ²(14) > 106.679)|0.0|Value close to 0 (e.g. under 0.05) means poor fitting. Value close to 1 can be an evidence for over-fitting or that sd are overestimated. It can be NaN (not a number) for dof=0|
  
<a name="measured-values"></a>
### Measured values
  

|name|sd|value|estimated|(estimated - value)/sd|
| :---: | :---: | :---: | :---: | :---: |
|a|0.02004|0.254|0.241|-0.61778|
|b|0.02004|0.053|0.157|5.17031|
|c|0.02004|0.448|0.334|-5.72715|
|d|0.02004|0.244|0.268|1.17461|
|k|0.02|0.497|0.496|-0.08698|
|l|0.02|0.503|0.504|0.08698|
|r|0.03|0.502|0.482|-0.66398|
|s|0.03|0.498|0.518|0.66398|
|e|0.03|0.51|0.482|-0.91611|
|f|0.03|0.49|0.518|0.91611|
|t|0.03|0.496|0.564|2.25558|
|u|0.03|0.504|0.436|-2.25558|
|m|0.03|0.496|0.398|-3.25308|
|n|0.03|0.504|0.602|3.25308|
|o|0.01|0.259|0.266|0.70237|
|p|0.01|0.495|0.514|1.92553|
|q|0.01|0.246|0.22|-2.6279|
|g|0.01|0.117|0.114|-0.30342|
|h|0.01|0.355|0.355|0.07771|
|i|0.01|0.387|0.392|0.45884|
|j|0.01|0.141|0.139|-0.23313|
  
<a name="redundant-measurements"></a>
### Redundant measurements

|Measurement|Measured Value|Formula Applied|Δ=FA-MV|
| :---: | :---: | :---: | :---: |
|b|0.053|0.218|0.164|
|d|0.244|0.292|0.047|
|e|0.51|0.516|0.007|
|f|0.49|0.484|-0.007|
|k|0.497|0.511|0.014|
|l|0.503|0.489|-0.014|
|m|0.496|0.471|-0.025|
|n|0.504|0.529|0.025|
|r|0.502|0.548|0.046|
|s|0.498|0.452|-0.046|
  
<a name="isotopomers"></a>
### Isotopomers

|variable|value|sd|
| :---: | :---: | :---: |
|000|0.114|0.00849|
|001|0.078|0.01289|
|100|0.152|0.01016|
|101|0.138|0.01089|
|010|0.125|0.00862|
|011|0.081|0.00679|
|110|0.173|0.008|
|111|0.139|0.00594|
  
<a name="cumomers"></a>
### Cumomers

|variable|value|sd|
| :---: | :---: | :---: |
|1xx|0.602|0.01814|
|x1x|0.518|0.00997|
|xx1|0.436|0.01233|
|11x|0.312|0.00889|
|1x1|0.277|0.01401|
|x11|0.22|0.00598|
|111|0.139|0.00594|
  
<a name="emu-values"></a>
### EMU values

|variable|Exx|xEx|xxE|EEx|ExE|xEE|EEE|
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
|M+0|0.398|0.482|0.564|0.192|0.239|0.266|0.114|
|M+1|0.602|0.518|0.436|0.496|0.484|0.514|0.355|
|M+2||||0.312|0.277|0.22|0.392|
|M+3|||||||0.139|
  
<a name="emu-standard-deviations"></a>
### EMU standard deviations

|variable|Exx|xEx|xxE|EEx|ExE|xEE|EEE|
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
|M+0|0.01814|0.00997|0.01233|0.01692|0.01215|0.00748|0.00849|
|M+1|0.01814|0.00997|0.01233|0.01693|0.01858|0.00729|0.00904|
|M+2||||0.00889|0.01401|0.00598|0.0085|
|M+3|||||||0.00594|
  
<a name="measurable-combinations"></a>
### Measurable combinations

|Combination|Accumulated|value|sd|
| :---: | :---: | :---: | :---: |
|000|000|0.114|0.00849|
|001|001|0.078|0.01289|
|100|100|0.152|0.01016|
|101|101|0.138|0.01089|
|010|010|0.125|0.00862|
|011|011|0.081|0.00679|
|110|110|0.173|0.008|
|111|111|0.139|0.00594|
