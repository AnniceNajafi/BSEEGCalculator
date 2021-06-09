# BSEEG Calculator (An application for quantifying delirium in rodents)
This is an application for the detection of delirium in mice using 24-hour EEG recordings. Upon submission, the power spectral density of the EEG file is calculated. The transformed data is then parsed into 10-minute windows and filtered for the range of 3-3.5Hz(diffuse slow waves) and 9.5-10 Hz(Alpha wave) and the ratio of the two is found and averaged per hour. Finally, the results of all days is calculated and binded then a plot is normalized based on the first day. 
Please try the application at "https://shinozakilab.shinyapps.io/bseegcalculator/" and use username:'research' and password:'1234' to login.

