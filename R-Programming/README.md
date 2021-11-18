# R-Project

We were given several datasets and we had to condect an Exploratory Data Analysis and create an algorithm to predict solar energy in the state of Oklahoma. 

We had 98 solar stations and had to forecast their future sun exposure based on several variables. 

In this EDA we use LeafLet to plot the Stations on a map, and then conducted a clustering analysis in order to divide the different weather stations into clusters based on their geographical position (latitude, longitude, elevation).

Then we used Time Series models (SARIMA model) in order to try to forecast the solar energy in each zone.
The idea is to consider that close stations should receive the same amount of Sun, then predict the solar reception of the geographical mean of each cluster and consider it is the same for the all cluster. 

At the en we created a basic Support Vector Machines algorithm to predict the missing values.

We had to render our report as an R Markdown that you can dowload in order visually see our clusters, map and time series. 
