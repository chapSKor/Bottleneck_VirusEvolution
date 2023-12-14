library(pracma)
library(deSolve)

####Note: The full forloop takes many hours to complete. 
##### pchip sometimes returns an error, however, it does not stop the code and the system of ODE's will still 
# solve. 

y0 = 1.33e5
BottleNeckparams_k3 = c(
  
  alpha = 7.9e-5,
  D = 0.32, 
  E = 5, 
  B = 15.7,
  C = 15, 
  k = 3
)
ComputeR0 = function(alpha, D, E, B, C){
  
  R0 = B*alpha*y0*((3*E)^(3))/(D*(alpha*y0 + C)*((D + 3*E)^(3)) )
  
  return(R0)
}

R0_start = ComputeR0(BottleNeckparams_k3[1], BottleNeckparams_k3[2],
                     BottleNeckparams_k3[3], BottleNeckparams_k3[4], 
                     BottleNeckparams_k3[5])
Total = 50/0.5
Parameters_selection = data.frame(matrix(ncol = 5, nrow = Total))
colnames(Parameters_selection) = c("alpha", "D", "E", "B", "C")

paramCounter = 1
alpha_FINAL = numeric(1)
D_FINAL = numeric(1)
E_FINAL = numeric(1)
B_FINAL = numeric(1)
C_FINAL = numeric(1)

alpha_start = BottleNeckparams_k3[1]
D_start = BottleNeckparams_k3[2]
E_start = BottleNeckparams_k3[3]
B_start = BottleNeckparams_k3[4]
C_start = BottleNeckparams_k3[5]

for(i in 1:Total){#Begin outer loop
  if(i == 1){
    alpha_FINAL = alpha_start
    D_FINAL = D_start
    E_FINAL = E_start
    B_FINAL = B_start
    C_FINAL = C_start
  }
  #Desired R0 is:
  #Calculate percent increase
  R0_final = R0_start + ((i*0.5)/100)*R0_start
  
  #E, B, and alpha you add. C and D you substract (to increase R0). 
  for(j in 1:1000000){#Go for as long as you need to find the parameter change required then save and break
    
    alpha_FINAL = alpha_FINAL + alpha_FINAL/1000
    #We only want to examine alpha here (not the other params)
    R0 = ComputeR0(alpha_FINAL, D_start, E_start, B_start, C_start)
    if(R0 >= R0_final){
      break
    }
  }
  for(k in 1:1000000){#D
    D_FINAL = D_FINAL - D_FINAL/1000
    R0 = ComputeR0(alpha_start, D_FINAL, E_start, B_start, C_start)
    if(R0 >= R0_final){
      break
    }
  }
  for(q in 1:1000000){#E
    
    E_FINAL = E_FINAL + E_FINAL/1000
    
    R0 = ComputeR0(alpha_start, D_start, E_FINAL, B_start, C_start)
    if(R0 >= R0_final){
      break
    }
  }
  for(w in 1:1000000){#B
    B_FINAL = B_FINAL + B_FINAL/1000
    R0 = ComputeR0(alpha_start, D_start, E_start, B_FINAL, C_start)
    if(R0 >= R0_final){
      break
    }
  }
  for(j in 1:1000000){#C
    C_FINAL = C_FINAL - C_FINAL/1000
    R0 = ComputeR0(alpha_start, D_start, E_start, B_start, C_FINAL)
    if(R0 >= R0_final){
      break
    }
  }
  
  #Save final param sets
  Parameters_selection$alpha[paramCounter] = alpha_FINAL
  Parameters_selection$D[paramCounter] = D_FINAL
  Parameters_selection$E[paramCounter] = E_FINAL #Only use 1:12, maybe 13 max
  Parameters_selection$B[paramCounter] = B_FINAL
  Parameters_selection$C[paramCounter] = C_FINAL
  
  #Reset dummy params 
  #alpha_FINAL = numeric(1)
  #D_FINAL = numeric(1)
  #E_FINAL = numeric(1)
  #B_FINAL = numeric(1)
  #C_FINAL = numeric(1)
  
  paramCounter = paramCounter + 1
  print(i)
}#End outer loop

###Plot parameter percent difference as a function of selective coffecient perfect increase. 

Parameters_PercentDifference = data.frame(matrix(ncol = 5, nrow = Total))
colnames(Parameters_PercentDifference) = c("alpha", "D", "E", "B", "C")

counter = 1

for(i in 1:100){
  
  Parameters_PercentDifference$alpha[counter] = abs((Parameters_selection$alpha[counter] - alpha_start)/alpha_start)*100
  Parameters_PercentDifference$D[counter] = abs((Parameters_selection$D[counter] - D_start)/D_start)*100
  Parameters_PercentDifference$E[counter] = abs((Parameters_selection$E[counter] - E_start)/E_start)*100
  Parameters_PercentDifference$B[counter] = abs((Parameters_selection$B[counter] - B_start)/B_start)*100
  Parameters_PercentDifference$C[counter] = abs((Parameters_selection$C[counter] - C_start)/C_start)*100
  
  counter = counter + 1 
  
}

colfunc2 <- colorRampPalette(c("#674A40","#50A3A4", "#FCAF38", "#F95335","#B4C540", "#575A6C", "#3686C9", "#E80566", "#BA53B3", "#FAE442"))
colours2 = colfunc2(10)

# Will produce Fig 2B
Select_change = seq(0.5, 50, 0.5) #X-axis
length(Select_change)
quartz()
maxRange = max(Parameters_PercentDifference)
minRange = min(Parameters_PercentDifference)
plotseq = seq(1, length(Select_change), 4)
plot(Select_change[plotseq], Parameters_PercentDifference$alpha[plotseq], col = "#F95335", pch = 8,
     cex = 1.5, ylim = c(0,175), yaxt = "n", xaxt = "n", xlab = "R0 increase (%)", ylab = "Parameter percent difference (|%|)")
axis(1,las= 1, tck=0.02,mgp=c(3, .5, 0), cex.axis = 1.2)
axis(2,las= 1, tck=0.02,mgp=c(3, .5, 0), cex.axis = 1.2)
points(Select_change[plotseq], Parameters_PercentDifference$D[plotseq], col = "#FAE442", pch = 19,
       cex = 1.5)
points(Select_change, Parameters_PercentDifference$E, col = "#BA53B3", pch = 18,
       cex = 2.3)
points(Select_change[plotseq], Parameters_PercentDifference$B[plotseq], col = "#B4C540", pch = 15,
       cex = 1.5)
points(Select_change[plotseq], Parameters_PercentDifference$C[plotseq], col = "#3686C9", pch = 17,
       cex = 1.5)

########################################################################################################
########### Set up the preamble and solve the stochastic bottleneck system for all the bottleneck sizes
#########################################################################################################

#Convert all parameters to hours
## !!!
Parameters_selection = Parameters_selection/24 #Only run this once!!! 
## !!!

times_k3 = seq(0,500,1) #In units of hours currently

#### Just infectious viral load #####
#### Jones WT params convert to hours for higher resolution#####
alpha_start = 7.9e-5/24
B_start = 15.7/24
C_start = 15/24
D_start = 0.32/24
E_start = 5/24
#Set up and solve the TEIV model for k=3

BottleNeckparams_k3 = c(
  
  #alpha = Parameters_selection$alpha[q],
  alpha = alpha_start,
  B = B_start,
  C = C_start, 
  D = D_start, 
  E = E_start, 
  k = 3
)
initConds_k3 = c(yT = 1.33e5, y1 = 0, y2 = 0, y3 = 0,yB = 1/30, v = 0)
bottleNeck_k3 = function(t,x,params){
  
  yT = x[1]
  y1 = x[2]
  y2 = x[3]
  y3 = x[4]
  yB = x[5]
  v = x[6]
  
  with(as.list(params),{
    dyT = -alpha*yT*v
    dy1 = alpha*yT*v - (D + k*E)*y1 #Eclipse stage 1
    dy2 = k*E*y1 - (D + k*E)*y2     #Eclispie stage 2
    dy3 = k*E*y2 - (D + k*E)*y3 
    dyB = k*E*y3 - D*yB
    dv = -C*v + B*yB - alpha*yT*v   
    
    list(c(dyT, dy1, dy2, dy3, dyB, dv))
  })  
}

out_k3_WT = as.data.frame(lsoda(initConds_k3,times_k3,bottleNeck_k3,BottleNeckparams_k3))

TAU = which.max(out_k3_WT$v)
#quartz()
#plot(out_k3_WT$v)
#Fcon = 2/(max(out_k3_WT$v))
BottleneckSize = c(5,50,5)
for(p in BottleneckSize){
  Fcon = p/(30*out_k3_WT$v[TAU])
  
  #Parameters_selection = Parameters_selection/24
  #Parameters_selection is calculated above
  Total = 50/0.5
  P_selection = data.frame(matrix(ncol = 5, nrow = Total))
  S_selection = data.frame(matrix(ncol = 5, nrow = Total))
  colnames(P_selection) = c("alpha", "D", "E", "B", "C")
  colnames(S_selection) = c("alpha", "D", "E", "B", "C")
  
  #These are (1-X(tm)):
  AllExtinctionProbs_alpha = data.frame(matrix(ncol = 100, nrow = TAU-3))
  AllExtinctionProbs_B = data.frame(matrix(ncol = 100, nrow = TAU-3))
  AllExtinctionProbs_C = data.frame(matrix(ncol = 100, nrow = TAU-3))
  AllExtinctionProbs_D = data.frame(matrix(ncol = 100, nrow = TAU-3))
  AllExtinctionProbs_E = data.frame(matrix(ncol = 13, nrow = TAU-3))
  ## So I need to calculate P for each of the 100 set of params in Parameters_selection
  mu = 1.3e-6 #Mutation rate
  initCondsIAVPGF = c(x1 = 1-Fcon, x2 = 1, x3 = 1,x4 = 1, x5 = 1, tau = TAU)
  tmLIM = TAU-3
  yB = out_k3_WT$yB[1:tmLIM]
  counterP = 1
  counter = 1
  #saveMAX = numeric(100)
  
  #Try chaning TAU back when starting again. 
  
  for(q in 1:100){#Begin outer loop
    
    alpha = alpha_start
    alpha_mutant = Parameters_selection$alpha[q]
    
    BottleNeckparams_PGF_COVID_C = c(
      C = Parameters_selection$C[q], 
      B = B_start,
      D = D_start, 
      E = E_start
    )
    BottleNeckparams_PGF_COVID_B = c(
      C = C_start, 
      B = Parameters_selection$B[q],
      D = D_start, 
      E = E_start
    )
    BottleNeckparams_PGF_COVID_D = c(
      C = C_start, 
      B = B_start,
      D = Parameters_selection$D[q], 
      E = E_start
    )
    BottleNeckparams_PGF_COVID_E = c(
      C = C_start, 
      B = B_start,
      D = D_start, 
      E = Parameters_selection$E[q]
    )
    BottleNeckparams_PGF_COVID_alpha = c(
      C = C_start, 
      B = B_start,
      D = D_start, 
      E = E_start
    )
    
    ExtinctionProbabilities = data.frame(matrix(ncol = 5, nrow = TAU-3))
    colnames(ExtinctionProbabilities) = c("B", "C", "D", "E", "alpha")
    for(i in 1:(TAU)){
      if(i >= (TAU-2)){
        break
      }
      t0 = i #Time mutation appears
      timesIAV = seq(t0,(TAU),1) 
      A = alpha*out_k3_WT$yT[t0:(TAU)]
      A_mutant = alpha_mutant*out_k3_WT$yT[t0:(TAU)]
      xi = seq(t0,(TAU),1)
      #Update TEIV past tm and cmopare ratios. 
      bottleNeckIAVPGF = function(t,x,params){
        
        x1 = x[1]
        x2 = x[2]
        x3 = x[3]
        x4 = x[4]
        x5 = x[5]
        tau = x[6]
        
        with(as.list(params),{
          
          dx1 = pchip(xi,A,tau)*x2 + C - (pchip(xi,A,tau) + C)*x1
          dx2 = -(3*E + D)*x2 + 3*E*x3 + D
          dx3 = -(3*E + D)*x3 + 3*E*x4 + D
          dx4 = -(3*E + D)*x4 + 3*E*x5 + D
          dx5 = B*x1*x5 + D - (D + B)*x5
          dtau = -1 
          #So the time in pchip should be x[4].
          #tau is the time that needs to be put into the attachment rate. 
          list(c(dx1, dx2, dx3,dx4,dx5, dtau))
        })  
      }
      
      #Only call this function for the alpha mutant
      bottleNeckIAVPGF_alphaMutant = function(t,x,params){
        
        x1 = x[1]
        x2 = x[2]
        x3 = x[3]
        x4 = x[4]
        x5 = x[5]
        tau = x[6]
        
        with(as.list(params),{
          
          dx1 = pchip(xi,A_mutant,tau)*x2 + C - (pchip(xi,A_mutant,tau) + C)*x1
          dx2 = -(3*E + D)*x2 + 3*E*x3 + D
          dx3 = -(3*E + D)*x3 + 3*E*x4 + D
          dx4 = -(3*E + D)*x4 + 3*E*x5 + D
          dx5 = B*x1*x5 + D - (D + B)*x5
          dtau = -1 
          #So the time in pchip should be x[4].
          #tau is the time that needs to be put into the attachment rate. 
          list(c(dx1, dx2, dx3,dx4,dx5, dtau))
        })  
      }
      
      outPGF_B = as.data.frame(lsoda(initCondsIAVPGF,timesIAV,bottleNeckIAVPGF,BottleNeckparams_PGF_COVID_B))
      outPGF_C = as.data.frame(lsoda(initCondsIAVPGF,timesIAV,bottleNeckIAVPGF,BottleNeckparams_PGF_COVID_C))
      outPGF_D = as.data.frame(lsoda(initCondsIAVPGF,timesIAV,bottleNeckIAVPGF,BottleNeckparams_PGF_COVID_D))
      outPGF_alpha = as.data.frame(lsoda(initCondsIAVPGF,timesIAV,bottleNeckIAVPGF_alphaMutant,BottleNeckparams_PGF_COVID_alpha))
      
      #Can only do a subset for E:
      if(q <= 13){
        outPGF_E = as.data.frame(lsoda(initCondsIAVPGF,timesIAV,bottleNeckIAVPGF,BottleNeckparams_PGF_COVID_E))
      }
      
      end_B = length(outPGF_B$x1)
      end_C = length(outPGF_C$x1)
      end_D = length(outPGF_D$x1)
      end_alpha = length(outPGF_alpha$x1)
      ExtinctionProbabilities$B[i] = outPGF_B$x1[end_B]
      ExtinctionProbabilities$C[i] = outPGF_C$x1[end_C]
      ExtinctionProbabilities$D[i] = outPGF_D$x1[end_D]
      ExtinctionProbabilities$alpha[i] = outPGF_alpha$x1[end_alpha]
      
      #Can only do a subset for E:
      if(q <= 13){
        end_E = length(outPGF_E$x1)
        ExtinctionProbabilities$E[i] = outPGF_E$x1[end_E]
      }
    }
    
    COVIDBottleneck_extinction_k3_B = 1 - ExtinctionProbabilities$B
    COVIDBottleneck_extinction_k3_C = 1 - ExtinctionProbabilities$C
    COVIDBottleneck_extinction_k3_D = 1 - ExtinctionProbabilities$D
    COVIDBottleneck_extinction_k3_alpha = 1 - ExtinctionProbabilities$alpha
  
    #These will be of length tau-3
    #Save Extinction probabilities
    AllExtinctionProbs_B[1:length(COVIDBottleneck_extinction_k3_B),counter] = COVIDBottleneck_extinction_k3_B
    AllExtinctionProbs_C[1:length(COVIDBottleneck_extinction_k3_C),counter] = COVIDBottleneck_extinction_k3_C
    AllExtinctionProbs_D[1:length(COVIDBottleneck_extinction_k3_D),counter] = COVIDBottleneck_extinction_k3_D
    AllExtinctionProbs_alpha[1:length(COVIDBottleneck_extinction_k3_alpha),counter] = COVIDBottleneck_extinction_k3_alpha
    
    #1 - X(tm) = survivalP_mean, below from line 402
    #survivalP = AllExtinctionProbs_alpha[,counter]
    survivalP_B = COVIDBottleneck_extinction_k3_B
    survivalP_C = COVIDBottleneck_extinction_k3_C
    survivalP_D = COVIDBottleneck_extinction_k3_D
    survivalP_alpha = COVIDBottleneck_extinction_k3_alpha
    
    #Budding cells, mean response from line 323
    
    transRate_nu_B = B_start*yB*mu*survivalP_B
    transRate_nu_C = B_start*yB*mu*survivalP_C
    transRate_nu_D = B_start*yB*mu*survivalP_D
    transRate_nu_alpha = B_start*yB*mu*survivalP_alpha
    
    #S; the expected number of times that the mutation of interest occurs de novo, 
    #over the course of the infection, and survives to be transmitted to the next host:
    S_B = sum(transRate_nu_B)
    S_C = sum(transRate_nu_C)
    S_D = sum(transRate_nu_D)
    S_alpha = sum(transRate_nu_alpha)
    
    S_selection$B[counterP] = S_B
    S_selection$C[counterP] = S_C
    S_selection$D[counterP] = S_D
    S_selection$alpha[counterP] = S_alpha
    #P; the probability that at least one copy of the mutation of interest arises 
    #de novo during the course of the infection, survives the bottleneck, and is 
    #transmitted to the new host:
    P_B = 1 - exp(-S_B)
    P_C = 1 - exp(-S_C)
    P_D = 1 - exp(-S_D)
    P_alpha = 1 - exp(-S_alpha)
    #P_selection$alpha[counterP] = P
    P_selection$B[counterP] = P_B
    P_selection$C[counterP] = P_C
    P_selection$D[counterP] = P_D
    P_selection$alpha[counterP] = P_alpha 
    
    #Again do E separately for a subset of iterations
    if(q <= 13){
      COVIDBottleneck_extinction_k3_E = 1 - ExtinctionProbabilities$E
      AllExtinctionProbs_E[1:length(COVIDBottleneck_extinction_k3_E),counter] = COVIDBottleneck_extinction_k3_E
      survivalP_E = COVIDBottleneck_extinction_k3_E
      transRate_nu_E = B_start*yB*mu*survivalP_E
      S_E = sum(transRate_nu_E)
      S_selection$E[counterP] = S_E
      P_E = 1 - exp(-S_E)
      P_selection$E[counterP] = P_E
    }
    
    counter = counter + 1
    counterP = counterP + 1
    
    #Update progress to the user
    if(q %% 5 == 0){
      String1_q = "Bottleneck size "
      comma = ","
      paramNum = "parameter set "
      Pchar = as.character(p)
      qchar = as.character(q)
      complete = "(complete)"
      String_update_q = paste(String1_q,Pchar,comma,paramNum,qchar,complete, sep = "")
      print(String_update_q) 
    }
    
  }#End outer loop and dont save (dont compile saving code below)
  
  Pchar = as.character(p)
  
  # String1_B = "AllExtinctionProb_B_bottleneck"
  # String1_D = "AllExtinctionProb_D_bottleneck"
  # String1_C = "AllExtinctionProb_C_bottleneck"
  # String1_E = "AllExtinctionProb_E_bottleneck"
  # String1_alpha = "AllExtinctionProb_alpha_bottleneck"
  # String2 = "_WTJonesk3Params_R019_TAUmax.csv"
  # 
  # dataFile_B = paste(String1_B,Pchar,String2,sep = "")
  # dataFile_D = paste(String1_D,Pchar,String2,sep = "")
  # dataFile_C = paste(String1_C,Pchar,String2,sep = "")
  # dataFile_E = paste(String1_E,Pchar,String2,sep = "")
  # dataFile_alpha = paste(String1_alpha,Pchar,String2,sep = "")
  # 
  # 
  # write.table(AllExtinctionProbs_B, dataFile_B)
  # write.table(AllExtinctionProbs_D, dataFile_D)
  # write.table(AllExtinctionProbs_C, dataFile_C)
  # write.table(AllExtinctionProbs_E, dataFile_E)
  # write.table(AllExtinctionProbs_alpha, dataFile_alpha)
  # 
  # String1_P = "AllP_bottleneck"
  # String2_P = "_WTJonesk3Params_R019_TAUmax.csv"
  # 
  # dataFile_P = paste(String1_P,Pchar, String2_P, sep = "")
  # write.table(P_selection, dataFile_P)
  
  String_update = "Bottleneck size complete: "
  String_update_2 = paste(String_update,Pchar, sep = "")
  print(String_update_2)
}#uncomment this loop end if you want to save everything
##################################
##################################



