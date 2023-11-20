module dynamics
push!(LOAD_PATH, pwd())
using LinearAlgebra
export initialcoherent
export timeevolution
#export fotoc
#export SurvivalP

function initialcoherent(xi::Float64,pi::Float64,hbar::Float64,Nmax::Int64)
	 alabs=((xi^2+pi^2)*(1/(2*hbar)))^(1/2)
	 al=(1/(2*hbar)^(1/2))*(xi+im*pi)
	 cs=[(al^n)/((factorial(big(n)))^(1/2)) for n in 0:Nmax]
	 cs = exp(-alabs^2/2)*cs
	 cs=[round(cs[i],digits=20) for i in 1:length(cs)]
	 csf=[convert(Complex{Float64},cs[i]) for i in 1:length(cs)]
	 return csf
	 end


function timeevolution(cx::Vector{Float64},xlim::Vector{Float64},cp::Vector{Float64},plim::Vector{Float64},psi0::Vector{Complex{Float64}},Tl::Vector{Float64},hbar::Float64,Nmax::Int64)
	 dar=[0.0 for i in 1:Nmax]
	 d=[0.0 for i in 1:Nmax+1]
	 dab=[(i)^(1/2) for i in 1:Nmax]
	 a=Array(Tridiagonal(dar,d,dab))
	 ad=a'
	 xo=((10000/2.)^(1/2))*(a + ad);
	 po=-im*((10000/2.)^(1/2))*(a - ad)
	 H=sum(cx[i]*xo^(i) for i in 1:length(cx))+sum(cp[i]*po^(i) for i in 1:length(cp))
	 pmax=plim[2]
	 pmin=plim[1]
	 xmax=xlim[2]
	 xmin=xlim[1]
	 NN=50
	 intp=(pmax-pmin)/NN
	 intx=(xmax-xmin)/NN
	 open("tevolution.dat","w") do io
	 for t in Tl
	   evol=exp(-im*H*t/10000)
	   x=xmin
	   p=pmin
	   for i in 1:NN+1
 	     for j in 1:NN+1
	     	alabs=((x^2+p^2)*(1/(2*hbar)))^(1/2)
	        al=(1/(2*hbar)^(1/2))*(x-im*p)
		psi0a=Array{Complex{Float64}}(undef,1,Nmax+1) 
	        for k in 1:Nmax+1
		   fac=(al^(k-1))/((factorial(big((k-1))))^(1/2))
		   facf=convert(Complex{Float64},fac)
		   psi0a[1,k]=facf
		end
	        psi0a = exp(-alabs^2/2)*psi0a
		psi0t=evol*psi0
		wvec=(psi0a*psi0t)
		w=abs2(wvec[1])
   	 	println(io,round(x,digits=16)," ",round(p,digits=16)," ",w)
		p=p+intp
	     end
	   p=pmin
	   x=x+intx
	   end
	 end
	 end
	 println("-------------Go to file tevolution.dat for dynamics----------------")
	 println("The file contains ",length(Tl)," screenshots of the dynamics    ")
	 println("Each screenshot corresponds to 2601 rows of the data file          ")
	 println("-------------------------------------------------------------------")
	 return "done"
	 end





#csi=initialcoherent(0.5,0.0,0.1,50)
#L1=[0,0.5]
#L2=[0,0.5]
#tlist=[0.0]
#xlim=[-4.0,4.0]
#plim=[-4.0,4.0]
#Tlist=[0.0,10]
#tt=timeevolution(L1,xlim,L2,plim,csi,Tlist,0.1,50)
#println(tt)

end