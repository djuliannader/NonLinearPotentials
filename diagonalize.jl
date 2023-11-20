module diagonalize
push!(LOAD_PATH, pwd())
using LinearAlgebra
export diagonalization

function diagonalization(cx::Vector{Float64},cp::Vector{Float64},Nmax::Int64,hbar::Float64)
	 dar=[0.0 for i in 1:Nmax]
	 d=[0.0 for i in 1:Nmax+1]
	 dab=[(i)^(1/2) for i in 1:Nmax]
	 a=Array(Tridiagonal(dar,d,dab))
	 ad=a'
	 xo=((hbar/2.)^(1/2))*(a + ad);
	 po=-im*((hbar/2.)^(1/2))*(a - ad)
	 H=sum(cx[i]*xo^(i) for i in 1:length(cx))+sum(cp[i]*po^(i) for i in 1:length(cp))
	 return eigvals(H)
	 end
	 
#L1=[0,-10.0,0,1]
#L2=[0,0.5]
#ex=diagonalization(L1,L2,100,0.5)
#println(length(ex))



end