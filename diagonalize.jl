module diagonalize
push!(LOAD_PATH, pwd())
using LinearAlgebra
export diagonalization
export expectationvalueH

function diagonalization(cx::Vector{Float64},cp::Vector{Float64},Nmax::Int64,hbar::Float64,flag::Int64)
	 dar=[0.0 for i in 1:Nmax]
	 d=[0.0 for i in 1:Nmax+1]
	 dab=[(i)^(1/2) for i in 1:Nmax]
	 a=Array(Tridiagonal(dar,d,dab))
	 ad=a'
	 xo=((hbar/2.)^(1/2))*(a + ad);
	 po=-im*((hbar/2.)^(1/2))*(a - ad)
	 H=sum(cx[i]*xo^(i) for i in 1:length(cx))+sum(cp[i]*po^(i) for i in 1:length(cp))
	 if flag==1
	    energies=eigvals(H)
	    return energies
	 end
	 if flag==2
	    states=eigvecs(H)
	    return states
	 end
	 end
	 
function expectationvalueH(Nmax::Int64,hbar::Float64,psi::Vector{ComplexF64},cx::Vector{Float64},cp::Vector{Float64})
	 dar=[0.0 for i in 1:Nmax]
	 d=[0.0 for i in 1:Nmax+1]
	 dab=[(i)^(1/2) for i in 1:Nmax]
	 a=Array(Tridiagonal(dar,d,dab))
	 ad=a'
	 xo=((hbar/2.)^(1/2))*(a + ad);
	 po=-im*((hbar/2.)^(1/2))*(a - ad)
	 H=sum(cx[i]*xo^(i) for i in 1:length(cx))+sum(cp[i]*po^(i) for i in 1:length(cp))
	 psia=Array{Complex{Float64}}(undef,1,Nmax+1) 
	 for k in 1:Nmax+1
             psia[1,k]=psi[k]
         end
	 psih=H*psi
	 val=psia*psih
	 return val
         end

function probabilitiesE(Nmax::Int64,hbar::Float64,ent::Float64,psi::Vector{ComplexF64},cx::Vector{Float64},cp::Vector{Float64})
	 dar=[0.0 for i in 1:Nmax]
	 d=[0.0 for i in 1:Nmax+1]
	 dab=[(i)^(1/2) for i in 1:Nmax]
	 a=Array(Tridiagonal(dar,d,dab))
	 ad=a'
	 xo=((hbar/2.)^(1/2))*(a + ad);
	 po=-im*((hbar/2.)^(1/2))*(a - ad)
	 H=sum(cx[i]*xo^(i) for i in 1:length(cx))+sum(cp[i]*po^(i) for i in 1:length(cp))
	 psia=Array{Complex{Float64}}(undef,1,Nmax+1) 
	 for k in 1:Nmax+1
             psia[1,k]=psi[k]
         end
	 eigsorted=eigvals(H)
         statessorted=eigvecs(H)
	 ik=1
	 suma=0
	 eps=0.00000000001
	 while (real(eigsorted[ik])<=ent+eps)
	    state=[statessorted[ik,i] for i in 1:Nmax+1]
	    coef=psia*state
	    suma=suma + abs2(coef[1])
	    println("---state below"," ",ik," ",real(eigsorted[ik]))
	    ik=ik+1
	 end
	 return 1-suma
	 end
end