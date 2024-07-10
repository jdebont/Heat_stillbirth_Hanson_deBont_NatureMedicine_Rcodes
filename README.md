# A time-stratified, case-crossover study of heat exposure and perinatal mortality from 16 hospitals in Sub-Saharan Africa

In this repository we have uploaded all the Rcodes we used to analyse the ALERT data focusing on the associations between heat exposure and perinatal mortality and antepartum and intrapartum stillbirths

Data Used:

- This study included singleton births from a prospective observational study in 16 hospitals in Benin, Malawi, Tanzania and Uganda collected as part of the Action Leveraging Evidence to Reduce perinatal morTality and morbidity (ALERT) study

- We obtained daily mean temperatures at 2 meters from the European Centre for Medium-Range Weather Forecasts (ECMWF) at a 9 x 9 km resolution (0.1° x 0.1°).

- We applied a time-stratified case-crossover design (lag 0-6)


Description of our codes:

01_data management : Data management creating the case-crossover dataset.

02_temperature_stillbirth.R: Analysis between mean temperature (lag 0-6) and perinatal mortality (Figure 2). The same codes can be used to reproduce the results for maximum and minimum temperature.

03_temperature_stillbirth_meta-regression: Meta-regression between mean temperature (lag 0-6) and perinatal mortality (Figure 3). The same codes can be used to reproduce the results for maximum and minimum temperature.

04_summarize_Rcode: This codes was used to summarize all the results in tables and figures for the manuscript and supplementary material. We have added a PDF with the results from this code.

On request to jeroen.de.bont@ki.se, all the sensitivity analysis codes can be provided.

