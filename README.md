# RoboArms
### Background
Linear Discriminant Analysis or LDA is one of many approaches used in 
supervised machine learning to solve multi-class classification 
problems. It uses many different methods from statistics and data 
processing to make conclusions about given criteria using the input
data.

### Introduction
In this specific project, LDA is being used to classify
various postures of hand movements based on Electromyography (EMG) 
data gathered from the forearms of multiple subjects. The data 
classified was gathered using 8 EMG electrodes measuring 11 postures 
over 6 trials. The raw data was taken and passed through our 
proprietary open source algorithm to classify each posture with an 
intermediate level of accuracy. Future improvements can be made to 
significantly increase the accuracy. Below is a detailed UML of how
the program works.

## UML Diagrams

### Data Pre-Processing
![UML](https://github.com/sayounan/RoboArms/blob/main/Media/LDA%20Landscape/Page-1.png)
### Data Processing
![UML](https://github.com/sayounan/RoboArms/blob/main/Media/LDA%20Landscape/Page-2.png)
### Confusion Matrix
![UML](https://github.com/sayounan/RoboArms/blob/main/Media/Picture1.png)
This matrix shows the accuracy of our model by showing the number of calculated
classifications of a posture relative to the expected posture. For example, the TAB
posture was correctly classified 80.25% of the time on the high end while HC was
classified correctly only 38.07% of the time relative to the expected.
