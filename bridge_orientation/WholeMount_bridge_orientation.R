 # Christoph Jüschke
 # IMBA - Institute of Molecular Biotechnology. All rights reserved.
 # Edited by Sabrina Syamsul
 #
 # Save compiled csv/excel as tab delimited file (tab delimited .txt)
 #
 # executing the script in powershell:
 # Rscript .\calculate_angles_rVers450.R .\fileName.txt
 #
 # e. g.:
 # Rscript .\calculate_angles_rVers450.R .\Raw_angles.txt
 #
 # output file: vers450angles_fileName.txt

 set.seed(1)
 cmd_args = commandArgs(T)
 file = cmd_args[1]
 base = basename(file)
 cat("analysing ", file, "...\n", sep="")
 Rx.stats <- function(x) {
  n = nrow(x)
  s = colSums(abs(x))
  r = sqrt(crossprod(s)) # length of the resultant vector
  m = s / r # mean resultant vector
  rbar = round(r / n, 15) # length of the mean resultant vector
  var = 1 - rbar # variance
  aSD = sqrt(2 * var) # angular std
  v = rep(0, ncol(x))
  v[ncol(x)] = 1
  mean = abs(asin(crossprod(m, v)))
  data.frame(n, r, rbar, var, aSD, mean)
 }
 rawdata = read.table(file, sep="\t", as.is=T, header=T)
 n = dim(rawdata)[1]
 attach(rawdata)
 Cen = array(data=c(C1_x, C1_y, C1_z, C2_x, C2_y, C2_z), 
            dim=c(n, 3, 2), 
            dimnames=list(cell=1:n, coord=c("x", "y", "z"), points=c("C1", "C2")))
 Pla = array(data=c(P1_x, P1_y, P1_z, P2_x, P2_y, P2_z, P3_x, P3_y, P3_z, P4_x, P4_y, P4_z, P5_x, P5_y, P5_z), 
            dim=c(n, 3, 5),
            dimnames=list(cell=1:n, coord=c("x", "y", "z"), points=paste("P", 1:5, sep="")))
 detach(rawdata)
 C_v = Cen[,,"C2"] - Cen[,,"C1"]
 C_n = C_v / sqrt(apply(C_v, 1, crossprod))
 centroid = apply(Pla, c(1,2), mean)
 Pl = sweep(Pla, c(1,2), centroid, "-")
 P_n = apply(Pl, 1, function(x) { s=svd(x); s$u[,match(min(s$d), s$d)] })
 dimnames(P_n) = list(coord=c("X","Y","Z"), cell=1:n)
 alp5 = abs(asin(colSums(P_n * t(C_n))))
 al4 = c()
 P4_n = c()
 for (i in 1:5) {
  Pla4 = Pla[,,-i]
  centroid = apply(Pla4, c(1,2), mean)
  Pl = sweep(Pla4, c(1,2), centroid, "-")
  P_n = apply(Pl, 1, function(x) { s=svd(x); s$u[,match(min(s$d), s$d)] })
  dimnames(P_n) = list(coord=c("X","Y","Z"), cell=1:n)
  alp4 = abs(asin(colSums(P_n * t(C_n))))
  P4_n = c(P4_n, P_n)
  al4 = cbind(al4, alp4)
 }
 P4_n = array(P4_n, dim=c(3, n, 5), 
  dimnames=list(coord=c("x","y","z"), dat=1:n, points=paste("dP", 1:5, sep="")))
 normSD = c()
 phiSD = c()
 for (i in 1:n) {
  normSD = c(normSD, Rx.stats(t(P4_n[,i,]))$aSD)
  phiSD = c(phiSD, Rx.stats(matrix(cbind(cos(al4[i,]),sin(al4[i,])), nrow=5))$aSD)
 }
 results = data.frame(rawdata[,c("treatment","individual","bridge")], 
                     alpha = round(alp5 * 180/pi, 2), 
                     alphaSD = round(phiSD * 180/pi, 2), 
                     normSD = round(normSD * 180/pi, 2))
 write.table(results, paste("vers450angles_", base, sep=""), row.names=F, quote=F, sep="\t")
