pUAST.seq = read.csv("/Users/johncsantiago/Documents/GitHub/WhartonLab/GeneralDataFiles/pUAST.Sequence.csv", header = F)

pUAST.seq = cbind(pUAST.seq, pUAST.seq)
colnames(pUAST.seq) = c("top", "bottom")

pUAST.seq$bottom[pUAST.seq$top == "c"] = "g"
pUAST.seq$bottom[pUAST.seq$top == "g"] = "c"
pUAST.seq$bottom[pUAST.seq$top == "a"] = "t"
pUAST.seq$bottom[pUAST.seq$top == "t"] = "a"

three.prime = c(6061,6293)
five.prime = c(1,587)

plot(x = NA,
     y = NA,
     type = 'n',
     ylim = c(0,3),
     xlim = c(0, 9100),
     #xaxt = "n",
     ylab="",
     yaxt = "n",
     xlab = NA)
lines(x = c(1,9053),
      y = c(2,2),
      lwd = 2)
lines(x = c(1,9053),
      y = c(1,1),
      lwd = 2)

lines(x = three.prime,
      y = c(2,2),
      lwd = 4,
      col = "red")

lines(x = three.prime,
      y = c(1,1),
      lwd = 4,
      col = "red")

lines(x = five.prime,
      y = c(2,2),
      lwd = 4,
      col = "green")

lines(x = five.prime,
      y = c(1,1),
      lwd = 4,
      col = "green")
