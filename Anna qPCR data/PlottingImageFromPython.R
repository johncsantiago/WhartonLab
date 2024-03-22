library(dplyr)
temp = read.csv("/Users/johncsantiago/Documents/Image.csv", header = F)
temp = as.matrix(temp)

##temp = t(temp)
row.names(temp) = c(1:nrow(temp))
##x = c(1:nrow(temp))
##temp = cbind(x,temp)
##colnames(temp) = c("x", 1:(ncol(temp)-1))

colnames(temp)= c(1:ncol(temp))

##temp = data.frame(temp)
df = temp %>%
  reshape2::melt() %>%
  rename(y=Var1,x=Var2)

##df$value[df$value >-50] = 50
##df$value[df$value < (-25)] = -50
##df$color = df$value
##df$color[df$color ==(-50)] = 'orange'
##df$color[df$color ==50] = 'navy'

df$color = ceiling(df$value + abs(min(df$value))+1)

color.pallette = colorRampPalette(c('orange','blue3', "navy"))(max(df$color))

df$color = color.pallette[df$color]

plot(x = df$x,
     y = df$y,
     col = df$color,
     ylim = c(max(df$y),1),
     ylab = NA,
     xlab = NA)
