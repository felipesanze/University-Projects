# R-Project

We were given several datasets and we had to conduct an Exploratory Data Analysis and create a model to predict solar energy reception in the state of Oklahoma. 

We had 98 solar stations (sensors) and had to forecast their future sun exposure based on several variables. 

In this EDA we used LeafLet to plot the Stations on a map, and then conducted a clustering analysis in order to divide the different weather stations into clusters based on their 3D geographical position (latitude, longitude, elevation).

Then we used Time Series models (SARIMA model) to forecast the solar energy in each zone.
The idea is to consider that close stations should receive the same amount of Sun, then predict the solar reception of the geographical mean of each cluster and consider it is the same for the all cluster.

At the en we created a basic Support Vector Machine algorithm to predict the missing values.

We had to render our report as an R Markdown that you can dowload in order visually see our clusters, map and time series. 
