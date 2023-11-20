module diagonalize
push!(LOAD_PATH, pwd())
using LinearAlgebra
export diagonalization

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
	 



end