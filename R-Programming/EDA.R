library(data.table)
library(Hmisc)
library(corrplot)
library(leaflet)
library(fBasics)
library(forecast) 
library(ggplot2)
library(dplyr)
library(tibble)
library(gridExtra)
library(tidyverse)  # data manipulation
library(cluster)  # clustering algorithms
library(factoextra) # clustering algorithms & visualization

### path ####
folder_path <- "~/IE/Term 2/R prgramming/Group Assignment"

###Read Data###

dat <- readRDS(file.path(folder_path, "solar_dataset.RData"));
dim(dat)
head(dat)

df <- as.data.frame(dat)
dt <- as.data.table(dat)

station_dat <- read.csv(file.path(folder_path, "station_info.csv"))
dim(station_dat)
head(station_dat)

df_stat <- as.data.frame(station_dat)
dt_stat <- as.data.table(station_dat)

additional_var_dat <- readRDS(file.path(folder_path, "additional_variables.RData"));
dim(additional_var_dat)
head(additional_var_dat)


############################## EDA ##########################
# Function to count NAs

count_nas <- function(x){
  ret <- sum(is.na(x));
  return(ret);
}

# Call function for each column
sapply(dat, count_nas);

# In percentages
sapply(dat, function(x){100*sum(is.na(x))/length(x)}); # All weather stations have about 26% of missing values! Which we will attempt to predict

### Statistics
sapply(dat, summary)

###removing rows with missing DATA 

df_red <- df[1:5113, ]
sapply(df_red, function(x){sum(is.na(x))});

######## Visualize data distribution of weather stations ######

boxplot(dat[,2:11])
boxplot(dat[,12:22])
boxplot(dat[,23:33])
boxplot(dat[,34:44])
boxplot(dat[,45:55]) #dat[,46] seems with outliers
boxplot(dat[,56:66])
boxplot(dat[,67:77])
boxplot(dat[,78:88])
boxplot(dat[,89:99])

#total weather stations in one plot
boxplot(dat[,2:99])


#Visualize data distribution of PC's top 100
boxplot(dat[,100:110])
boxplot(dat[,111:121])
boxplot(dat[,122:132])
boxplot(dat[,133:143])
boxplot(dat[,144:154]) 
boxplot(dat[,155:165])
boxplot(dat[,166:176])
boxplot(dat[,177:188])
boxplot(dat[,188:198])

boxplot(dat[,100:200])



############## Finding OUTLIERS ###############

outlier <-
  function (x, opposite = FALSE, logical = FALSE) 
  {
    if (is.matrix(x)) 
      apply(x, 2, outlier, opposite = opposite, logical = logical)
    else if (is.data.frame(x)) 
      sapply(x, outlier, opposite = opposite, logical = logical)
    else {
      if (xor(((max(x,na.rm=TRUE) - mean(x,na.rm=TRUE)) < (mean(x,na.rm=TRUE) - min(x,na.rm=TRUE))),opposite)) 
      {
        if (!logical) min(x,na.rm=TRUE)
        else x == min(x,na.rm=TRUE)
      }
      else 
      {
        if (!logical) max(x,na.rm=TRUE)
        else x == max(x,na.rm=TRUE)
      }
    } 
  }

outlier(df_red[,100:102])
sapply(df_red[,2:199], outlier )


##### chisq test #######
chisq.out.test <-
  function (x,variance=var(x),opposite=FALSE) 
  {
    DNAME <- deparse(substitute(x))
    x <- sort(x[complete.cases(x)])
    n <- length(x);
    if (xor(((x[n] - mean(x)) < (mean(x) - x[1])),opposite)) {
      alt = paste("lowest value",x[1],"is an outlier");
      chi <- (mean(x)-x[1])^2/variance
    }
    else{
      alt = paste("highest value",x[n],"is an outlier");
      chi <- (x[n]-mean(x))^2/variance
    }
    pval <- 1-pchisq(chi,1)
    RVAL <- list(statistic = c("X-squared" = chi), alternative = alt, p.value = pval, method = "chi-squared test for outlier", 
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }

sapply(df_red$PC2[1:10], chisq.out.test)


##### REMOVING OUTLIERS #######
rm.outlier <-
  function (x, fill = FALSE, median = FALSE, opposite = FALSE) 
  {
    if (is.matrix(x)) 
      apply(x, 2, rm.outlier, fill = fill, median = median, opposite = opposite)
    else if (is.data.frame(x)) 
      as.data.frame(sapply(x, rm.outlier, fill = fill, median = median, opposite = opposite))
    else {
      res <- x
      if (!fill) res[-which(x == outlier(x,opposite))]
      else {
        if (median) res[which(x == outlier(x,opposite))]<-median(x[-which(x == outlier(x,opposite))])
        else res[which(x == outlier(x,opposite))]<-mean(x[-which(x == outlier(x,opposite))])
        res
      }
    }
    
  }

rm.outlier(df_red$PC2)
outlier(df_red$PC2)

sapply(df_red[,-1], rm.outlier)

#the only problem is that the function removes values once per iteration.
#would be better to try to create a loop to remove top 10 outliers for example
df_clean<- rm.outlier(df_red[,101])
outlier(df_clean[101])

boxplot(df_clean)
boxplot(df_red$PC2)

###### Plot the time series ######
y<-df_red[,2][1:5113]
nlags = 500
ts.plot(y)  
par(mfrow=c(2,1))
acf(y, nlags)  
pacf(y, nlags)

ndiffs(y, alpha=0.05, test=c("adf"))

#We can clearly notice the data is stationary, and it has seasonality

# estimating a potential forecasting time series model
s= 365 #Seasonal patterns can be shown anually in the data, since it is daily data s = 365
fit<-arima(y,order=c(0,1,3),seasonal=list(order=c(0,1,0),period=s)) 
fit  


par(mfrow=c(3,1))
ts.plot(fit$residuals)
acf(fit$residuals)
pacf(fit$residuals)    

Box.test(fit$residuals,lag=35)


# testing for normality 
shapiro.test(fit$residuals)  # 95% confidence intervals are robust for any kind of distribution

hist(fit$residuals,prob=T,ylim=c(0,1),xlim=c(mean(fit$residuals)-3*sd(fit$residuals),mean(fit$residuals)+3*sd(fit$residuals)),col="red")
lines(density(fit$residuals),lwd=2)
mu<-mean(fit$residuals)
sigma<-sd(fit$residuals)
x<-seq(mu-3*sigma,mu+3*sigma,length=100)
yy<-dnorm(x,mu,sigma)
lines(x,yy,lwd=2,col="blue") # theoretical

# Forecasting
y.pred<-predict(fit,n.ahead=1796)
y.pred$pred # point predictions
y.pred$se  # standard error for the point predictions to compute confidence intervals

# Comparing the real values with the predicted ones
real<-df_red[,2][5001:5113] 
c=1:113
plot(c,real,type="b",col="red")
lines(c,y.pred$pred,col="blue",type="b")
legend("bottomleft",c("real","forecast"),
       col = c("red","blue"),pch = c(1,1),bty ="n" )

# plotting real data with point predictions

new <- c(y,y.pred$pred) # real data + predicted values

plot.ts(new,main="Predictions",
        ylab="Solar Energy",col=3,lwd=2) # time series plot
lines(y,col=4,lwd=2) # for the second series
legend("topleft",legend=c("Predictions","Historical"),col=c(3,4),
       bty="n",lwd=2)



########### RESCALING DATA ############
##Create a function to rescale the data i.e substract
##the mean and divide by the standard deviation. 

tipify <- function(x){
  mu <- mean(x, na.rm = TRUE);
  s <- sd(x, na.rm = TRUE);
  
  # Special condition for constant variables
  s[s == 0] <- 1;
  
  # Tipify
  x <- (x - mu) / s;
}



sapply(df_red, mean, na.rm = TRUE);
sapply(df_red, sd, na.rm = TRUE);

df_tipified <- cbind(df_red[1], sapply(df_red[-1], tipify))

sapply(df_tipified, mean, na.rm = TRUE);
sapply(df_tipified, sd, na.rm = TRUE);

###changing date column from character to date type#####

lct <- Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME", "C")
df_tipified$Date <- as.Date(df_tipified[["Date"]], "%Y%m%d")
class(df_tipified$Date)


#### GGPLOT2 time series for columns in example #####
time_plot_HOBA <- ggplot(df_tipified,
                    aes(x = Date, y = HOBA)) + geom_point(na.rm = TRUE)+
                   labs(title="Weather Station data", subtitle="From HOBA weather station", 
                    y="Solar energy", x="Time");

time_plot_HOBA<- time_plot_HOBA + scale_x_date(date_breaks = "1 year", date_labels = "%Y")
time_plot_HOBA

time_plot_WILB <- ggplot(df_tipified,
                         aes(x = Date, y = WILB)) + geom_point(colour = "red") +
                         labs(title="Weather Station data", subtitle="From WILB weather station", 
                         y="Solar energy", x="Time");


time_plot_WILB <- time_plot_WILB + scale_x_date(date_breaks = "1 year", date_labels = "%Y")

time_plot_WILB




###Combining two graphs for better viisualization###

grid.arrange(time_plot_HOBA, time_plot_WILB, ncol = 2)



# Function to count Top n most frequent values
top_n_frequencies <- function(x, n = 5){
  return(names(sort(table(x), decreasing = TRUE)[1:n]));
}

#Most frequent values on ACME
sapply(dt[,2], top_n_frequencies);

#Most frequent values on ADAX
sapply(dt[,3], top_n_frequencies);

# For numerics or integer: boxplot

boxplot(dat$PC2) #Make sure to remove outliers from PC2
boxplot(dat$PC1)





############################## Correlation Analysis ######################


corr = cor(df_red[,-1]) # -1 here means we look at all columns except the first column
corr

corr2 <- rcorr(as.matrix(df_red[,-1]))
corr2

corrplot(corr, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

res = cor(df_tipified[,-1]) # -1 here means we look at all columns except the first column
res

res2 <- rcorr(as.matrix(df_tipified[,-1]))
res2
class(res2)

corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

############################## Dimension reduction & Visualization of correlations ######################

#We have identified PC1, 2, 3, 5, 7, 8, 10, 12, 13, 16, 22, 28 are highly correlated
# with weather stations (which are all intercorrelated between them)

df_dim_red <- df_tipified[, 98:128] #we know we have included uninteresting variables
head(df_dim_red)

#correlation analysis on the REDUCED dataset

res = cor(df_dim_red[,]) 
res

res2 <- rcorr(as.matrix(df_dim_red[,]))
res2

corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)



#Plots

par(mfrow=c(2,1))
plot(dat$PC1, type = 'l')
plot(dat$PC2, type = 'l')
plot(dat$ADAX, type = 'l')



##### Leaflet  WEATHER STATION MAP LOCATION #######
install.packages("leaflet")
library(leaflet)


m <- leaflet() %>%
  addTiles() %>%  # Add default OpenStreetMap map tiles
  addMarkers(lng=-98.02325, lat=34.80833, popup="ACME")%>%
  addMarkers(lng=-96.66909, lat=34.79851, popup="ADAX")%>%
  addMarkers(lng=-99.33808, lat=34.58722, popup="ALTU")%>%
  addMarkers(lng=-98.29216, lat=34.91418, popup="APAC")%>%
  addMarkers(lng=-99.90308, lat=36.07204, popup="ARNE")%>%
  addMarkers(lng=-100.53012, lat=36.80253, popup="BEAV")%>%
  addMarkers(lng=-99.05847, lat=35.40185, popup="BESS")%>%
  addMarkers(lng=-95.86621, lat=35.96305, popup="BIXB")%>%
  addMarkers(lng=-97.25452, lat=36.75443, popup="BLAC")%>%
  addMarkers(lng=-102.49713, lat=36.69256, popup="BOIS")%>%
  addMarkers(lng=-96.63121, lat=35.17156, popup="BOWL")%>%
  addMarkers(lng=-97.69394, lat=36.41201, popup="BREC")%>%
  addMarkers(lng=-96.35404, lat=35.78050, popup="BRIS")%>%
  addMarkers(lng=-99.64101, lat=36.83129, popup="BUFF")%>%
  addMarkers(lng=-96.81046, lat=36.63459, popup="BURB")%>%
  addMarkers(lng=-97.26918, lat=33.89376, popup="BURN")%>%
  addMarkers(lng=-99.27059, lat=35.59150, popup="BUTL")%>%
  addMarkers(lng=-97.00330, lat=34.84970, popup="BYAR")%>%
  addMarkers(lng=-99.34652, lat=36.02866, popup="CAMA")%>%
  addMarkers(lng=-96.33309, lat=34.60896, popup="CENT")%>%
  addMarkers(lng=-96.80407, lat=35.65282, popup="CHAN")%>%
  addMarkers(lng=-98.36274, lat=36.74813, popup="CHER")%>%
  addMarkers(lng=-99.72790, lat=35.54615, popup="CHEY")%>%
  addMarkers(lng=-97.91446, lat=35.03236, popup="CHIC")%>%
  addMarkers(lng=-95.32596, lat=34.65657, popup="CLAY")%>%
  addMarkers(lng=-95.24870, lat=34.22321, popup="CLOU")%>%
  addMarkers(lng=-94.84896, lat=35.68001, popup="COOK")%>%
  addMarkers(lng=-95.88553, lat=36.90987, popup="COPA")%>%
  addMarkers(lng=-96.32027, lat=33.92075, popup="DURA")%>%
  addMarkers(lng=-98.03654, lat=35.54848, popup="ELRE")%>%
  addMarkers(lng=-99.80344, lat=35.20494, popup="ERIC")%>%
  addMarkers(lng=-95.65707, lat=35.30324, popup="EUFA")%>%
  addMarkers(lng=-98.49766, lat=36.26353, popup="FAIR")%>%
  addMarkers(lng=-96.42777, lat=36.84053, popup="FORA")%>%
  addMarkers(lng=-99.14234, lat=36.72562, popup="FREE")%>%
  addMarkers(lng=-98.46607, lat=35.14887, popup="FTCB")%>%
  addMarkers(lng=-101.60130, lat=36.60183, popup="GOOD")%>%
  addMarkers(lng=-97.47978, lat=35.84891, popup="GUTH")%>%
  addMarkers(lng=-95.64047, lat=35.74798, popup="HASK")%>%
  addMarkers(lng=-98.48151, lat=35.48439, popup="HINT")%>%
  addMarkers(lng=-99.05283, lat=34.98971, popup="HOBA")%>%
  addMarkers(lng=-99.83331, lat=34.68550, popup="HOLL")%>%
  addMarkers(lng=-101.22547, lat=36.85518, popup="HOOK")%>%
  addMarkers(lng=-95.54011, lat=34.03084, popup="HUGO")%>%
  addMarkers(lng=-94.88030, lat=33.83013, popup="IDAB")%>%
  addMarkers(lng=-94.78287, lat=36.48210, popup="JAYX")%>%
  addMarkers(lng=-102.87820, lat=36.82937, popup="KENT")%>%
  addMarkers(lng=-97.76484, lat=34.52887, popup="KETC")%>%
  addMarkers(lng=-98.11139, lat=36.38435, popup="LAHO")%>%
  addMarkers(lng=-95.99716, lat=34.30876, popup="LANE")%>%
  addMarkers(lng=-96.94394, lat=34.03579, popup="MADI")%>%
  addMarkers(lng=-99.42398, lat=34.83592, popup="MANG")%>%
  addMarkers(lng=-97.21271, lat=36.06434, popup="MARE")%>%
  addMarkers(lng=-99.01109, lat=36.98707 , popup="MAYR")%>%
  addMarkers(lng=-95.78096, lat=34.88231, popup="MCAL")%>%
  addMarkers(lng=-97.74577, lat=36.79242, popup="MEDF")%>%
  addMarkers(lng=-98.56936, lat=34.72921, popup="MEDI")%>%
  addMarkers(lng=-94.84437, lat=36.88832, popup="MIAM")%>%
  addMarkers(lng=-97.95553, lat=35.27225, popup="MINC")%>%
  addMarkers(lng=-94.82275, lat=34.31072, popup="MTHE")%>%
  addMarkers(lng=-96.91035, lat=36.89810, popup="NEWK")%>%
  addMarkers(lng=-97.95202, lat=34.96774, popup="NINN")%>%
  addMarkers(lng=-95.60795, lat=36.74374, popup="NOWA")%>%
  addMarkers(lng=-96.49749, lat=36.03126, popup="OILT")%>%
  addMarkers(lng=-96.26265, lat=35.43172, popup="OKEM")%>%
  addMarkers(lng=-95.91473, lat=35.58211, popup="OKMU")%>%
  addMarkers(lng=-97.22924, lat=34.71550, popup="PAUL")%>%
  addMarkers(lng=-96.76986, lat=36.36114, popup="PAWN")%>%
  addMarkers(lng=-97.04831, lat=35.99865, popup="PERK")%>%
  addMarkers(lng=-95.27138, lat=36.36914, popup="PRYO")%>%
  addMarkers(lng=-98.96038, lat=35.89904, popup="PUTN")%>%
  addMarkers(lng=-97.15306, lat=36.35590, popup="REDR")%>%
  addMarkers(lng=-99.36001, lat=35.12275, popup="RETR")%>%
  addMarkers(lng=-97.58812, lat=34.19365, popup="RING")%>%
  addMarkers(lng=-94.79805, lat=35.43815, popup="SALL")%>%
  addMarkers(lng=-99.04030, lat=36.19033, popup="SEIL")%>%
  addMarkers(lng=-96.94822, lat=35.36492, popup="SHAW")%>%
  addMarkers(lng=-96.03706, lat=36.41530, popup="SKIA")%>%
  addMarkers(lng=-100.26192, lat=36.59749, popup="SLAP")%>%
  addMarkers(lng=-97.34146, lat=35.54208, popup="SPEN")%>%
  addMarkers(lng=-95.18116, lat=35.26527, popup="STIG")%>%
  addMarkers(lng=-97.09527, lat=36.12093, popup="STIL")%>%
  addMarkers(lng=-96.06982, lat=34.87642, popup="STUA")%>%
  addMarkers(lng=-96.95048, lat=34.56610, popup="SULP")%>%
  addMarkers(lng=-94.98671, lat=35.97235, popup="TAHL")%>%
  addMarkers(lng=-95.01152, lat=34.71070, popup="TALI")%>%
  addMarkers(lng=-99.13755, lat=34.43972, popup="TIPT")%>%
  addMarkers(lng=-96.67895, lat=34.33262, popup="TISH")%>%
  addMarkers(lng=-95.22094, lat=36.77536, popup="VINI")%>%
  addMarkers(lng=-97.52109, lat=34.98224, popup="WASH")%>%
  addMarkers(lng=-98.52615, lat=35.84185, popup="WATO")%>%
  addMarkers(lng=-97.98815, lat=34.16775, popup="WAUR")%>%
  addMarkers(lng=-98.77509, lat=35.50830, popup="WEAT")%>%
  addMarkers(lng=-94.64496, lat=36.01100, popup="WEST")%>%
  addMarkers(lng=-95.34805, lat=34.90092, popup="WILB")%>%
  addMarkers(lng=-94.68778, lat=34.98426, popup="WIST")%>%
  addMarkers(lng=-99.41682, lat=36.42329, popup="WOOD")%>%
  addMarkers(lng=-96.34222, lat=36.51806, popup="WYNO")
m 


# ----------- [5] Sation Zoning ------------
library(tidyverse)  # data manipulation
library(cluster)  # clustering algorithms
#install.packages("factoextra")
library(factoextra) # clustering algorithms & visualization

# ---------- [5.1] Scaling Data ----------
station_dat
new_dat <- station_dat[,2:4]
new_dat
scaled <- scale(new_dat)
scaled

# ---------- [5.2] Distance mapping between stations ----------
distance <- get_dist(scaled)
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

# ---------- [5.3] Preparing data ----------
new_dat <- station_dat
new_dat
rownames(new_dat) <- new_dat$stid
scaled
rownames(scaled) <- new_dat$stid
scaled

# ---------- [5.4] Plotting the clusters ----------
# ---------- [5.4.1] Scaled plotting ----------

k2 <- kmeans(scaled, centers = 2, nstart = 25)
k3 <- kmeans(scaled, centers = 3, nstart = 25)
k4 <- kmeans(scaled, centers = 4, nstart = 25)
k5 <- kmeans(scaled, centers = 5, nstart = 25)
k6 <- kmeans(scaled, centers = 6, nstart = 25)

# plots to compare
p1 <- fviz_cluster(k2, geom = "point", data = scaled) + ggtitle("k = 2")
p2 <- fviz_cluster(k3, geom = "point",  data = scaled) + ggtitle("k = 3")
p3 <- fviz_cluster(k4, geom = "point",  data = scaled) + ggtitle("k = 4")
p4 <- fviz_cluster(k5, geom = "point",  data = scaled) + ggtitle("k = 5")

grid.arrange(p1, p2, p3, p4, nrow = 2)

# ---------- [5.4.2] Original shape plotting ----------

# Clustering with k = 3 :
scaled %>%
  as_tibble() %>%
  mutate(cluster = k3$cluster,
         stations = row.names(scaled)) %>%
  ggplot(aes(elon, nlat, color = factor(cluster), label = stations)) +
  geom_text()

# Clustering with k = 4 :
scaled %>%
  as_tibble() %>%
  mutate(cluster = k4$cluster,
         stations = row.names(scaled)) %>%
  ggplot(aes(elon, nlat, color = factor(cluster), label = stations)) +
  geom_text()

# Clustering with k = 5 :
scaled %>%
  as_tibble() %>%
  mutate(cluster = k5$cluster,
         stations = row.names(scaled)) %>%
  ggplot(aes(elon, nlat, color = factor(cluster), label = stations)) +
  geom_text()

# Clustering with k = 6 :
scaled %>%
  as_tibble() %>%
  mutate(cluster = k6$cluster,
         stations = row.names(scaled)) %>%
  ggplot(aes(elon, nlat, color = factor(cluster), label = stations)) +
  geom_text()

# --------- [5.4.3] plotting the final model with the optimal number of clusters ----------
set.seed(123)
final <- kmeans(scaled, 5, nstart = 25)
print(final)

fviz_cluster(final, data = scaled)

# ---------- [5.5] Coordinates and Final plot ----------

# ---------- [5.5.1] Computes the coordinates of the mean of the clusters ---------- 
scaled %>%
  as_tibble() %>%
  mutate(cluster = final$cluster, stations = row.names(scaled)) %>%
  group_by(cluster) %>%
  summarise_all("mean")

#cluster   nlat   elon   elev stations
#<int>  <dbl>  <dbl>  <dbl>    <dbl>
#1       1 -0.573 -0.645  0.290       NA
#2       2  1.31  -2.53   3.48        NA
#3       3 -1.05   0.723 -0.700       NA
#4       4  0.879 -0.987  0.900       NA
#5       5  0.812  0.650 -0.490       NA

# ---------- [5.5.2] Plots the last clusters ----------
scaled %>%
  as_tibble() %>%
  mutate(cluster = final$cluster,
         stations = row.names(scaled)) %>%
  ggplot(aes(elon, nlat, color = factor(cluster), label = stations)) +
  geom_text()

# Gives the cluster to which every station belongs 
final$cluster

# ---------- [5.5.3] Coordinates as a DataFrame ---------- 

names(station_dat)
zones1 <- data.frame("zone1", -0.573, -0.645,    0.290)
names(zones1) <- c("Cluster", "nlat", "elon", "elev")
zones2 <- data.frame("zone2", 1.31, -2.53,    3.48)
names(zones2) <- c("Cluster", "nlat", "elon", "elev")
zones3 <- data.frame("zone3", -1.05, 0.723,    -0.700)
names(zones3) <- c("Cluster", "nlat", "elon", "elev")
zones4 <- data.frame("zone4", 0.879, -0.987,    0.900)
names(zones4) <- c("Cluster", "nlat", "elon", "elev")
zones5 <- data.frame("zone5", 0.812, 0.650,    -0.490)
names(zones5) <- c("Cluster", "nlat", "elon", "elev")
zone <- rbind(zones1, zones2, zones3, zones4, zones5)
zone
names(zone) <- c("stid", "nlat", "elon", "elev")
zone
# ---------- [5.6] Unscale the mean's coordinates of the clusters ---------- 

# DataFrame of the stations' coordinates :
scaled
# DataFrame of the mean of the cluster' coordinates :
zone

# The function rbind doesn't work : 
all <- rbind(scaled, zone)
all # the stations are back into their original scale

# Assuming linearity in the 3 axis x, y, z, we can find the unscanscaled values of the means of each clusters by doing : : 
cluster 1:
  FTCB unscaled : 35.14887  -98.46607  422.00
FTCB scaled: -0.47364496 -0.56683158  0.18439868

zone1 scaled: -0.573 -0.645  0.29
zone1 unscaled :  (-0.573*35.14887)/-0.47364496, (-0.645*-98.46607)/-0.56683158, (0.29*422.00)/0.18439868

cluster 2 :
  BOIS unscaled: 36.69256 -102.49713 1267.00
BOIS scaled:  1.25565352 -2.77722691  4.15064674

zone2 scaled: 1.310 -2.530  3.48
zone2 unscaled : (1.310*36.69256)/1.25565352, (-2.530*-102.49713)/-2.77722691, (3.48*1267.00)/4.15064674

cluster 3 : 
  CENT unscaled : 34.60896  -96.33309  208.00
Cent scaled : -1.07847205  0.60276873 -0.82007124

zone3: scaled: -1.050  0.723 -0.70
zone3 unscaled: (-1.050*34.60896)/-1.07847205, (0.723*-96.33309)/0.60276873, (-0.70*208.00)/-0.82007124

cluster 4:
  Wood unscaled = 36.42329,  -99.41682,  625.00
Wood scaled = 0.95400733, -1.08816675,  1.13723698

zone4 scaled = 0.879, -0.987,  0.90
zone4 unscaled = (0.879 * 36.42329)/0.95400733, (-0.987*-99.41682)/-1.08816675, (0.90*625.00)/1.13723698


cluster 5:
  SKIA unscaled :36.41530  -96.03706  282.00
SKIA scaled : 0.94505664  0.76509410 -0.47273117

zone5 scaled : 0.812  0.650 -0.49
zone5 unscaled : (0.812*36.41530)/0.94505664, (0.650*-96.03706)/0.76509410, (-0.49*282.00)/-0.47273117

# ---------- [5.6.1] Mean's coordinates of the clusters as a DataFrame ---------- 

coord1 <- data.frame("clust1", 42.52194, -112.0449,    663.6707)
names(coord1) <- c("stid", "nlat", "elon", "elev")  
coord2 <- data.frame("clust2", 38.28067, -93.3729,    1062.283)
names(coord2) <- c("stid", "nlat", "elon", "elev")
coord3 <- data.frame("clust3", 33.69527, -115.5482,    177.5456)
names(coord3) <- c("stid", "nlat", "elon", "elev")
coord4 <- data.frame("clust4", 33.55957, --90.17405,    494.6199)
names(coord4) <- c("stid", "nlat", "elon", "elev")
coord5 <- data.frame("clust5", 31.28831, -81.59008,    292.3014)
names(coord5) <- c("stid", "nlat", "elon", "elev")

coords <- rbind(coord1, coord2, coord3, coord4, coord5)
coords
names(coords) <- c("stid", "nlat", "elon", "elev")

library(data.table);

########################## [1] REGRESSION ############################

completed <- dat[1:5113]

to_predict <- dat[5114:6909]

########################## [1.1] LINEAR REGRESSION ############################

# Let's predict the variable HOBA!

########################## [1.1.1] 1 variable model using all data #############################
# build linear regression model to predict HOBA using as predictor variable PC1
model_l1_all <- lm(HOBA ~ PC1, data = dat);

# Check model info
print(model_l1_all); # HOBA = w*PC1 + (Intercept) = -5.344*PC1 + 37.285
summary(model_l1_all);

# Get model predictions
predictions <- predict(model_l1_all, newdata = dat);
manual_predictions <- dat$PC1*model_l1_all$coefficients[2] + model_l1_all$coefficients[1] ;
head(predictions);
head(manual_predictions); # High explainability (Accuracy-'explainability' tradeoff)

head(dat$HOBA);

# Get errors
errors <- predictions - dat$HOBA;
head(errors);

# Compute MAE
mae_l1_all <- round(mean(abs(errors)), 2);
sprintf("MAE = %s",mae_l1_all);

# Compute MSE
mse_l1_all <- round(mean(errors^2), 2);
sprintf("MSE = %s",mse_l1_all);


########################## [1.1.2] Train/test split #############################

# setting seed to reproduce results of random sampling
set.seed(1); 

# row indexes for training data (70%)
train_index <- sample(1:nrow(completed), 0.7*nrow(completed));  

# training data
train <- completed[train_index]; 

# test data
test  <- completed[-train_index ]

dim(completed); #rows and columns
dim(train);
dim(test);


########################## [1.1.3] 1 variable model using train/test #############################

# build linear regression model to predict HOBA using as predictor variable PC1
model_l1_train <- lm(HOBA ~ PC1, data = train);

# Get model predictions for train
predictions_train <- predict(model_l1_train, newdata = train);

# Get model predictions for test
predictions_test <- predict(model_l1_train, newdata = test);

# Get errors
errors_train <- predictions_train - train$HOBA;
errors_test <- predictions_test - test$HOBA;

# Compute Metrics
mse_train <- round(mean(errors_train^2), 2);
mae_train <- round(mean(abs(errors_train)), 2);

mse_test <- round(mean(errors_test^2), 2);
mae_test <- round(mean(abs(errors_test)), 2);

# Build comparison table
comp <- data.table(model = c("lm_1var"), 
                   mse_train = mse_train, mae_train = mae_train,
                   mse_test = mse_test, mae_test = mae_test);
comp;





########################## [1.2] More Complex Model. SVM ############################

########################## [1.2.1] Standard SVM ############################


install.packages("e1071");
library(e1071); # LIBSVM

### train standard SVM model
model <- svm(HOBA ~ PC1, data = train)  # . all columns except the target

# Check model info
print(model);


# Get model predictions
predictions_train <- predict(model, newdata = train); # There is no formula here! Less straightforward to explain
# (Accuracy-'explainability' tradeoff)
predictions_test <- predict(model, newdata = test);

# Get errors
errors_train <- predictions_train - train$HOBA;
errors_test <- predictions_test - test$HOBA;

# Compute Metrics
mse_train <- round(mean(errors_train^2), 2);
mae_train <- round(mean(abs(errors_train)), 2);

mse_test <- round(mean(errors_test^2), 2);
mae_test <- round(mean(abs(errors_test)), 2);

# Build comparison table
comp <- rbind(comp,
              data.table(model = c("standard svm"), 
                         mse_train = mse_train, mae_train = mae_train,
                         mse_test = mse_test, mae_test = mae_test));
comp; # Worse than linear regression. What is wrong





########################## [1.2.3] hyperparameter optimization ############################
 # Best model in terms of error metrics!! More accuracy (Accuracy-'explainability' tradeoff)
