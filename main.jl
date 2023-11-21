push!(LOAD_PATH, pwd())
import reading
import diagonalize
import dynamics

println("Non Linear Potentials 1D ")
println("Diagonalization in the Fock Space")

# Reading data from input file
#------------------------------------
open("input.dat") do f
 K1=readline(f)
 K2=readline(f)
 K3=readline(f)
 K4=readline(f)
 K5=readline(f)
 K6=readline(f)
 K7=readline(f)
 K8=readline(f)
 K9=readline(f)
 K10=readline(f)
 K11=readline(f)
 K12=readline(f)
 K13=readline(f)
 K14=readline(f)
 K15=readline(f)
 K16=readline(f)
 K17=readline(f)
 K18=readline(f)
 K19=readline(f)
 K20=readline(f)
 K21=readline(f)
 K22=readline(f)
 K23=readline(f)
 K24=readline(f)
 K25=readline(f)
 K26=readline(f)


#------ converting string to a list
 xcl = reading.stringtofloatlist(K6)
 pcl = reading.stringtofloatlist(K8)
 xps = reading.stringtofloatlist(K10)
 pps = reading.stringtofloatlist(K12)
 ts  = reading.stringtofloatlist(K14)
 xpi = reading.stringtofloatlist(K20)
 xcl0= reading.stringtofloatlist(K22)
 pcl0= reading.stringtofloatlist(K24)

#------ converting string to float or Integer
 Nmax=parse(Int64,K2)   
 hbar=parse(Float64,K4)
 tmax=parse(Float64,K16)
 flagi=parse(Int64,K18)
 ki=parse(Int64,K26)
 
 

#------ calling method which perform diagonalization
diag=diagonalize.diagonalization(xcl,pcl,Nmax,hbar,1)

println("Ground state Energy:",diag[1])

#----- calling methods for the dynamics of a coherent states
if flagi==1
  ics = dynamics.initialcoherent(xpi[1],xpi[2],hbar,Nmax)
  phasespace = dynamics.timeevolution(xcl,xps,pcl,pps,ics,ts,hbar,Nmax)
  survival = dynamics.survivalp(xcl,pcl,ics,tmax,hbar,Nmax)
end

#----- calling methods for the dynamics of a quench
if flagi==2
   eigvecs0 = diagonalize.diagonalization(xcl0,pcl0,Nmax,hbar,2)
   iqs=[eigvecs0[ki,i] for i in 1:Nmax+1]
   phasespace = dynamics.timeevolution(xcl,xps,pcl,pps,iqs,ts,hbar,Nmax)
   survival = dynamics.survivalp(xcl,pcl,iqs,tmax,hbar,Nmax)
end

end # Final