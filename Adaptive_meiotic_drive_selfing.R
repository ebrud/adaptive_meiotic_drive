Adaptive.drive.self <- function(ks2, ke2, wAA, waa, S){
  
  #### ks2 is the mutant ratio in males, and is referred to as "km" in the article. ke2 is "kf".
  
  #### Initializing values assume:
  #### (1) mutant modifiers at initial frequency of 1e-4, 
  #### (2) Heterozygote fitness equals 1.
  #### (3) Resident Mendelian segregation schemes
  #### These can be altered to explore other scenarios.
  
  b0 <- 1e-4
  
  wAa <- 1
  r <- 1/2
  
  ks1 <- 0.5
  ke1 <- 0.5

  #### For a mixed-selfing population, the values for initial equilibrium are calculated below.
  if(0 < S && S < 1){
  qhat = -((2*waa - 4*S*waa + 3*S*(waa^2) + 6*wAA - 4*S*wAA - 12*waa*wAA + 12*S*waa*wAA + 2*(waa^2)*wAA - 6*S*(waa^2)*wAA - 4*(wAA^2) + S*(wAA^2) + 6*waa*(wAA^2) - 2*S*waa*(wAA^2) + waa*sqrt(8*S*(-1 + waa)*(-1 + wAA)*(-wAA + waa*(-1 + 2*wAA)) + (2 - S*wAA - 2*waa*wAA + S*waa*(-1 + 2*wAA))^2) - wAA*sqrt(8*S*(-1 + waa)*(-1 + wAA)*(-wAA + waa*(-1 + 2*wAA)) + (2 - S*wAA - 2*waa*wAA + S*waa*(-1 + 2*wAA))^2))/(4*(-1 + S)*(-(-2 + wAA)*wAA + (waa^2)*(-1 + 2*wAA) + 2*waa*(1 - 3*wAA + (wAA^2)))))
  phat = 1-qhat
  
  Fhat = -(-2 + S*waa + S*wAA + 2*waa*wAA - 2*S*waa*wAA + sqrt(8*S*(-1 + waa)*(-1 + wAA)*(-wAA + waa*(-1 + 2*wAA)) + (2 - S*wAA - 2*waa*wAA + S*waa*(-1 + 2*wAA))^2))/(4*(-1 + waa)*(-1 + wAA))
  }
  
  
  #### The equilibrium under obligate selfing:
  if(S == 1){ 
    qhat = ((-1 + waa)*(-1 + 2*wAA))/(2 - 3*wAA + waa*(-3 + 4*wAA))
    phat = 1-qhat
    
    Fhat = (waa + wAA - 2*waa*wAA)/(2 - 2*waa - 2*wAA + 2*waa*wAA)
    
  }

  
  #### The genotypic frequencies at the B locus, assuming initial Hardy-Weinberg proportions at B
  BB.0 <-(1- b0)^2
  bb.0 <- (b0)^2
  Bb.0 <- 1-(BB.0 + bb.0)
  
  #### Two-locus diploid genotypic frequencies, assuming initial probabilistic independence
  
  ABAB <- (phat^2 + phat*qhat*Fhat)*(BB.0)
  ABAb <- (phat^2 + phat*qhat*Fhat)*(Bb.0)
  AbAb <- (phat^2 + phat*qhat*Fhat)*(bb.0)
  ABaB <- 2*phat*qhat*(1-Fhat)*(BB.0)
  Abab <- 2*phat*qhat*(1-Fhat)*(bb.0)
  ABab <- phat*qhat*(1-Fhat)*(Bb.0)
  AbaB <- phat*qhat*(1-Fhat)*(Bb.0)
  aBaB <- (qhat^2 + phat*qhat*Fhat)*(BB.0)
  aBab <- (qhat^2 + phat*qhat*Fhat)*(Bb.0)
  abab <- (qhat^2 + phat*qhat*Fhat)*(bb.0)
  
  #### Pre-allocation of vectors recording the evolution of allele frequencies at the B and A loci
  p.b <- b0
  p.b.history <- vector(length = 1e5)
  p.b.history[1] <- b0
  gen <- 1
  
  p.a <- qhat
  p.a.history <- vector(length = 1e5)
  p.a.history[1] <- qhat
  
  
  #### Numerical iteration of the model for 100,000 generations
  while(gen < 1e5){
  
    Wmean = ABAB*(wAA) + ABAb*(wAA) + AbAb*(wAA) + ABaB*(wAa) + Abab*(wAa) + ABab*(wAa) + AbaB*(wAa) + aBaB*(waa) + aBab*(waa) + abab*(waa)

eggAB = (1/Wmean)*(ABAB*(wAA)*(1) + ABAb*(wAA)*(1/2) + AbAb*(wAA)*(0) + ABaB*(wAa)*(1 - ke1) + Abab*(wAa)*(0) + ABab*(wAa)*(1 - r)*(1 - ke2) + AbaB*(wAa)*(r)*(1 - ke2) + aBaB*(waa)*(0) + aBab*(waa)*(0) + abab*(waa)*(0))

eggAb = (1/Wmean)*(ABAB*(wAA)*(0) + ABAb*(wAA)*(1/2) + AbAb*(wAA)*(1) + ABaB*(wAa)*(0) + Abab*(wAa)*(1 - ke2) + ABab*(wAa)*(r)*(1 - ke2) + AbaB*(wAa)*(1 - r)*(1 - ke2) + aBaB*(waa)*(0) + aBab*(waa)*(0) + abab*(waa)*(0))

eggaB = (1/Wmean)*(ABAB*(wAA)*(0) + ABAb*(wAA)*(0) + AbAb*(wAA)*(0) + ABaB*(wAa)*(ke1) + Abab*(wAa)*(0) + ABab*(wAa)*(r)*(ke2) + AbaB*(wAa)*(1 - r)*(ke2) + aBaB*(waa)*(1) + aBab*(waa)*(1/2) + abab*(waa)*(0))

eggab = (1/Wmean)*(ABAB*(wAA)*(0) + ABAb*(wAA)*(0) + AbAb*(wAA)*(0) + ABaB*(wAa)*(0) + Abab*(wAa)*(ke2) + ABab*(wAa)*(1 - r)*(ke2) + AbaB*(wAa)*(r)*(ke2) + aBaB*(waa)*(0) + aBab*(waa)*(1/2) + abab*(waa)*(1))


spermAB = (1/Wmean)*(ABAB*(wAA)*(1) + ABAb*(wAA)*(1/2) + AbAb*(wAA)*(0) + ABaB*(wAa)*(1 - ks1) + Abab*(wAa)*(0) + ABab*(wAa)*(1 - r)*(1 - ks2) + AbaB*(wAa)*(r)*(1 - ks2) + aBaB*(waa)*(0) + aBab*(waa)*(0) + abab*(waa)*(0))

spermAb = (1/Wmean)*(ABAB*(wAA)*(0) + ABAb*(wAA)*(1/2) + AbAb*(wAA)*(1) + ABaB*(wAa)*(0) + Abab*(wAa)*(1 - ks2) + ABab*(wAa)*(r)*(1 - ks2) + AbaB*(wAa)*(1 - r)*(1 - ks2) + aBaB*(waa)*(0) + aBab*(waa)*(0) + abab*(waa)*(0))

spermaB = (1/Wmean)*(ABAB*(wAA)*(0) + ABAb*(wAA)*(0) + AbAb*(wAA)*(0) + ABaB*(wAa)*(ks1) + Abab*(wAa)*(0) + ABab*(wAa)*(r)*(ks2) + AbaB*(wAa)*(1 - r)*(ks2) + aBaB*(waa)*(1) + aBab*(waa)*(1/2) + abab*(waa)*(0))

spermab = (1/Wmean)*(ABAB*(wAA)*(0) + ABAb*(wAA)*(0) + AbAb*(wAA)*(0) + ABaB*(wAa)*(0) + Abab*(wAa)*(ks2) + ABab*(wAa)*(1 - r)*(ks2) + AbaB*(wAa)*(r)*(ks2) + aBaB*(waa)*(0) + aBab*(waa)*(1/2) + abab*(waa)*(1))


ABAB.p = (1 - S)*(spermAB*eggAB) + (1/Wmean)*S*(ABAB*(wAA)*(1) + ABAb*(wAA)*(1/4) + AbAb*(wAA)*(0) + ABaB*(wAa)*((1 - ks1)*(1 - ke1)) + Abab*(wAa)*(0) + ABab*(wAa)*((1 - r)^2)*((1 - ks2)*(1 - ke2)) + AbaB*(wAa)*(r^2)*((1 - ks2)*(1 - ke2)) + aBaB*(waa)*(0) + aBab*(waa)*(0) + abab*(waa)*(0))

ABAb.p = (1 - S)*((spermAB*eggAb) + (spermAb*eggAB)) + (1/Wmean)*S*(ABAB*(wAA)*(0) + ABAb*(wAA)*(1/2) + AbAb*(wAA)*(0) + ABaB*(wAa)*(0) + Abab*(wAa)*(0) + ABab*(wAa)*(2*r*(1 - r)*((1 - ks2)*(1 - ke2))) + AbaB*(wAa)*(2*r*(1 - r)*((1 - ks2)*(1 - ke2))) + aBaB*(waa)*(0) + aBab*(waa)*(0) + abab*(waa)*(0))

AbAb.p = (1 - S)*(spermAb*eggAb) + (1/Wmean)*S*(ABAB*(wAA)*(0) + ABAb*(wAA)*(1/4) + AbAb*(wAA)*(1) + ABaB*(wAa)*(0) + Abab*(wAa)*((1 - ks2)*(1 - ke2)) + ABab*(wAa)*((r^2)*((1 - ks2)*(1 - ke2))) + AbaB*(wAa)*(((1 - r)^2)*((1 - ks2)*(1 - ke2))) + aBaB*(waa)*(0) + aBab*(waa)*(0) + abab*(waa)*(0))

ABaB.p = (1 - S)*((spermAB*eggaB) + (spermaB*eggAB)) + (1/Wmean)*S*(ABAB*(wAA)*(0) + ABAb*(wAA)*(0) + AbAb*(wAA)*(0) + ABaB*(wAa)*(ks1*(1 - ke1) + (1 - ks1)*ke1) + Abab*(wAa)*(0) + ABab*(wAa)*(r*(1 - r)*(ks2*(1 - ke2) + (1 - ks2)*ke2)) + AbaB*(wAa)*(r*(1 - r)*(ks2*(1 - ke2) + (1 - ks2)*ke2)) + aBaB*(waa)*(0) + aBab*(waa)*(0) + abab*(waa)*(0))

Abab.p = (1 - S)*((spermAb*eggab) + (spermab*eggAb)) + (1/Wmean)*S*(ABAB*(wAA)*(0) + ABAb*(wAA)*(0) + AbAb*(wAA)*(0) + ABaB*(wAa)*(0) + Abab*(wAa)*(ks2*(1 - ke2) + (1 - ks2)*ke2) + ABab*(wAa)*(r*(1 - r)*(ks2*(1 - ke2) + (1 - ks2)*ke2)) + AbaB*(wAa)*(r*(1 - r)*(ks2*(1 - ke2) + (1 - ks2)*ke2)) + aBaB*(waa)*(0) + aBab*(waa)*(0) + abab*(waa)*(0))

ABab.p = (1 - S)*((spermAB*eggab) + (spermab*eggAB)) + (1/Wmean)*S*(ABAB*(wAA)*(0) + ABAb*(wAA)*(0) + AbAb*(wAA)*(0) + ABaB*(wAa)*(0) + Abab*(wAa)*(0) + ABab*(wAa)*(((1 - r)^2)*(ks2*(1 - ke2) + (1 - ks2)*ke2)) + AbaB*(wAa)*((r^2)*(ks2*(1 - ke2) + (1 - ks2)*ke2)) + aBaB*(waa)*(0) + aBab*(waa)*(0) + abab*(waa)*(0))

AbaB.p = (1 - S)*((spermAb*eggaB) + (spermaB*eggAb)) + (1/Wmean)*S*(ABAB*(wAA)*(0) + ABAb*(wAA)*(0) + AbAb*(wAA)*(0) + ABaB*(wAa)*(0) + Abab*(wAa)*(0) + ABab*(wAa)*((r^2)*(ks2*(1 - ke2) + (1 - ks2)*ke2)) + AbaB*(wAa)*(((1 - r)^2)*(ks2*(1 - ke2) + (1 - ks2)*ke2)) + aBaB*(waa)*(0) + aBab*(waa)*(0) + abab*(waa)*(0))

aBaB.p = (1 - S)*(spermaB*eggaB) + (1/Wmean)*S*(ABAB*(wAA)*(0) + ABAb*(wAA)*(0) + AbAb*(wAA)*(0) + ABaB*(wAa)*(ks1*ke1) + Abab*(wAa)*(0) + ABab*(wAa)*((r^2)*(ks2*ke2)) + AbaB*(wAa)*(((1 - r)^2)*(ks2*ke2)) + aBaB*(waa)*(1) + aBab*(waa)*(1/4) + abab*(waa)*(0))

aBab.p = (1 - S)*((spermaB*eggab) + (spermab*eggaB)) + (1/Wmean)*S*(ABAB*(wAA)*(0) + ABAb*(wAA)*(0) + AbAb*(wAA)*(0) + ABaB*(wAa)*(0) + Abab*(wAa)*(0) + ABab*(wAa)*(2*r*(1 - r)*(ks2*ke2)) + AbaB*(wAa)*(2*r*(1 - r)*(ks2*ke2)) + aBaB*(waa)*(0) + aBab*(waa)*(1/2) + abab*(waa)*(0))

abab.p = (1 - S)*(spermab*eggab) + (1/Wmean)*S*(ABAB*(wAA)*(0) + ABAb*(wAA)*(0) + AbAb*(wAA)*(0) + ABaB*(wAa)*(0) + Abab*(wAa)*(ks2*ke2) + ABab*(wAa)*(((1 - r)^2)*(ks2*ke2)) + AbaB*(wAa)*((r^2)*(ks2*ke2)) + aBaB*(waa)*(0) + aBab*(waa)*(1/4) + abab*(waa)*(1))

#### Updating genotypic frequencies
ABAB <- ABAB.p
ABAb <- ABAb.p
AbAb <- AbAb.p
ABaB <- ABaB.p
Abab <- Abab.p
ABab <- ABab.p
AbaB <- AbaB.p
aBaB <- aBaB.p
aBab <- aBab.p
abab <- abab.p

gen <- gen+1

p.b <-(sum(AbAb, Abab, abab))+0.5*sum(ABAb, ABab, AbaB,aBab)
p.a <-(sum(aBaB, aBab, abab))+0.5*sum(ABaB, Abab, ABab, AbaB)

p.b.history[gen] <-p.b
p.a.history[gen] <- p.a

  }

  #### Generate trajectory plots...

plot(p.b.history, type = "l", lwd = 4, ylim = c(0,1))
lines(p.a.history, col = "dark green", lty = 2, lwd = 4)

}
