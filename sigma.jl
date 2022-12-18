function fsigma(x::Float64,w::Float64,epsilon::Float64,gamma::Float64,s::Float64,T::Float64)
	gamma=(gamma^2)*(x^(2*s))
	epsilon=epsilon^2
	x2=x^2
	w2=(w*T)^2
    return  4*(2)*exp(w)*((exp(w)+1)^(-2))*(x+1)*gamma*((w2+gamma+epsilon+x2)^2)/((((w2+gamma+epsilon+x2)^2)-4*w2*(epsilon+x2))^2)
end
function selffsigmaw(x::Float64,a::Float64,b::Float64,epsilon::Float64,gamma::Float64,s::Float64,T::Float64,S::Float64,errors::Float64=10^-10)

h=(b-a)/2
s1=h/6*(fsigma(x,a,epsilon,gamma,s,T)+4*fsigma(x,a+h/2,epsilon,gamma,s,T)+fsigma(x,a+h,epsilon,gamma,s,T))
s2=h/6*(fsigma(x,a+h,epsilon,gamma,s,T)+4*fsigma(x,a+3*h/2,epsilon,gamma,s,T)+fsigma(x,b,epsilon,gamma,s,T))
if abs(S-s1-s2) <= errors
	I= (s1+s2)+(s1+s2-S)/15
else
	I=selffsigmaw(x,a,(a+b)/2,epsilon,gamma,s,T,s1,errors)+selffsigmaw(x,(a+b)/2,b,epsilon,gamma,s,T,s2,errors)
end
          return I
end


function selffsigmax(a::Float64,b::Float64,w::Float64,epsilon::Float64,gamma::Float64,s::Float64,T::Float64,S::Float64,errors::Float64=10^-10)

h=(b-a)/2
s1=h/6*(fsigma(a,w,epsilon,gamma,s,T)+4*fsigma(a+h/2,w,epsilon,gamma,s,T)+fsigma(a+h,w,epsilon,gamma,s,T))
s2=h/6*(fsigma(a+h,w,epsilon,gamma,s,T)+4*fsigma(a+3*h/2,w,epsilon,gamma,s,T)+fsigma(b,w,epsilon,gamma,s,T))
if abs(S-s1-s2) <= errors
	I= (s1+s2)+(s1+s2-S)/15
else
	I=selffsigmax(a,(a+b)/2,w,epsilon,gamma,s,T,s1,errors)+selffsigmax((a+b)/2,b,w,epsilon,gamma,s,T,s2,errors)
end
          return I
end

function selffsigmae(x::Float64,w::Float64,a::Float64,b::Float64,gamma::Float64,s::Float64,T::Float64,S::Float64,errors::Float64=10^-10,i::Float64=0.0)
i=i+1.0
h=(b-a)/2
s1=h/6*(fsigma(x,w,a,gamma,s,T)+4*fsigma(x,w,a+h/2,gamma,s,T)+fsigma(x,w,a+h,gamma,s,T))
s2=h/6*(fsigma(x,w,a+h,gamma,s,T)+4*fsigma(x,w,a+3*h/2,gamma,s,T)+fsigma(x,w,b,gamma,s,T))
if abs(S-s1-s2) <= errors
	I= (s1+s2)+(s1+s2-S)/15
else
	I=selffsigmae(x,w,a,(a+b)/2,gamma,s,T,s1,errors,i)+selffsigmae(x,w,(a+b)/2,b,gamma,s,T,s2,errors,i)
end
          #println("i:",i," ",errors,"  ",abs(S-s1-s2))
          return I
end

function selffsigmaestart(x::Float64,w::Float64,a::Float64,b::Float64,gamma::Float64,s::Float64,T::Float64)
	h=(b-a)
	#println(h/6*fsigma(a,w,epsilon,gamma,s,T)," ",4*fsigma((a+b)/2,w,epsilon,gamma,s,T)," ",fsigma(b,w,epsilon,gamma,s,T))
	return (h/6*(fsigma(x,w,a,gamma,s,T)+4*fsigma(x,w,(a+b)/2,gamma,s,T)+fsigma(x,w,b,gamma,s,T)))
end

function selffsigmaxstart(a::Float64,b::Float64,w::Float64,epsilon::Float64,gamma::Float64,s::Float64,T::Float64)
	h=(b-a)
	#println(h/6*fsigma(a,w,epsilon,gamma,s,T)," ",4*fsigma((a+b)/2,w,epsilon,gamma,s,T)," ",fsigma(b,w,epsilon,gamma,s,T))
	return (h/6*(fsigma(a,w,epsilon,gamma,s,T)+4*fsigma((a+b)/2,w,epsilon,gamma,s,T)+fsigma(b,w,epsilon,gamma,s,T)))
end

function selffsigmawstart(x::Float64,a::Float64,b::Float64,epsilon::Float64,gamma::Float64,s::Float64,T::Float64)
	h=(b-a)
	#println(fsigma(x,a,epsilon,gamma,s,T)," ",(a+b)/2," ",fsigma(x,(a+b)/2,epsilon,gamma,s,T)," ",b," ",fsigma(x,b,epsilon,gamma,s,T))
	return (h/6*(fsigma(x,a,epsilon,gamma,s,T)+4*fsigma(x,(a+b)/2,epsilon,gamma,s,T)+fsigma(x,b,epsilon,gamma,s,T)))
end


# println(selffsigmawstart(0.1,0.0,5.0,0.1,0.1,1.0,0.001))
# for  i in 5.0:1.0:100.0
#  	println(i," ",selffsigmawstart(0.1,0.0,i,0.1,0.1,1.0,0.001)," ",selffsigmaw(0.1,0.0,i,0.1,0.1,1.0,0.001,selffsigmawstart(0.1,0.0,i,0.1,0.1,1.0,0.001)))
# end

function fselffsigmaw(x::Float64,epsilon::Float64,gamma::Float64,s::Float64,T::Float64)

	return selffsigmaw(x,0.0,10.0,epsilon,gamma,s,T,selffsigmawstart(x,0.0,10.0,epsilon,gamma,s,T))
end

function selffsigmawestart(a::Float64,b::Float64,x::Float64,gamma::Float64,s::Float64,T::Float64)
	h=b-a
	return (h/6*(fselffsigmaw(x,a,gamma,s,T)+4*fselffsigmaw(x,(a+b)/2,gamma,s,T)+fselffsigmaw(x,b,gamma,s,T)))
end

#println(fselffsigmaw(0.1,0.1,0.1,1.0,0.001))

function selffsigmawe(x::Float64,a::Float64,b::Float64,gamma::Float64,s::Float64,T::Float64,S::Float64,errors::Float64=10^-10,i::Float64=0.0)
  h=(b-a)/2
s1=h/6*(fselffsigmaw(x,a,gamma,s,T)+4*fselffsigmaw(x,a+h/2,gamma,s,T)+fselffsigmaw(x,a+h,gamma,s,T))
s2=h/6*(fselffsigmaw(x,a+h,gamma,s,T)+4*fselffsigmaw(x,a+3*h/2,gamma,s,T)+fselffsigmaw(x,b,gamma,s,T))
if abs(S-s1-s2) <= errors
	I= (s1+s2)+(s1+s2-S)/15
else
	I=selffsigmawe(x,a,(a+b)/2,gamma,s,T,s1,errors)+selffsigmawe(x,(a+b)/2,b,gamma,s,T,s2,errors)
end
          return I

end
# for i in 1.0:1.0:200.0
# println(selffsigmawe(0.1,0.0,i,0.1,1.0,0.01,selffsigmawestart(0.0,i,0.1,0.1,1.0,0.01)))
# end

function fselffsigmawe(x::Float64,gamma::Float64,s::Float64,T::Float64)

	return selffsigmawe(x,0.0,10.0,gamma,s,T,selffsigmawestart(x,0.0,10.0,gamma,s,T))
end
#println(fselffsigmawe(0.1,0.1,1.0,0.01))


function selffsigmawexstart(a::Float64,b::Float64,gamma::Float64,s::Float64,T::Float64)
	h=b-a

	println(fselffsigmawe(a,gamma,s,T))
	println(4*fselffsigmawe((a+b)/2,gamma,s,T))
	println(fselffsigmawe(b,gamma,s,T))
	return (h/6*(fselffsigmawe(a,gamma,s,T)+4*fselffsigmawe((a+b)/2,gamma,s,T)+fselffsigmawe(b,gamma,s,T)))
end


function selffsigmawex(a::Float64,b::Float64,gamma::Float64,s::Float64,T::Float64,S::Float64,errors::Float64=10^-10,i::Float64=0.0)
  h=(b-a)/2
s1=h/6*(fselffsigmawe(a,gamma,s,T)+4*fselffsigmawe(a+h/2,gamma,s,T)+fselffsigmawe(a+h,gamma,s,T))
s2=h/6*(fselffsigmawe(a+h,gamma,s,T)+4*fselffsigmawe(a+3*h/2,gamma,s,T)fselffsigmawe(b,gamma,s,T))
if abs(S-s1-s2) <= errors
	I= (s1+s2)+(s1+s2-S)/15
else
	I= selffsigmawex(a,(a+b)/2,gamma,s,T,s1,errors)+ selffsigmawex((a+b)/2,b,gamma,s,T,s2,errors)
end
          return I

end

function logsigmawex(T::Float64,gamma::Float64,s::Float64,n::Int64=50000)
	b=((pi/4)/(10^-10))^(1/n)
	I=0
	for i in 0:n-1
		x=(0.5*10^-10)*(b^(i+1)+b^(i))
		
        I=I+(10^-10)*(b^(i+1)-b^(i))*fselffsigmawe(x,gamma,s,T)
    end
    return I
end
for i in 1.:1.:100.
println(i," ",fselffsigmawe(10^(-i),0.1,2.0,0.00001))
end
#println(logsigmawex(0.001,0.1,2.0))
#println(selffsigmawex(0.1,pi/4.,0.1,2.0,0.001,selffsigmawexstart(0.1,pi/4.,0.1,2.0,0.001)))
#println(selffsigmawexstart(0.1,pi/4.,0.1,2.0,0.01))
#println(selffsigmawex(0.1,0.0,i,0.1,1.0,0.01,selffsigmawexstart(0.0,i,0.1,0.1,1.0,0.01)))

#println(selffsigmaxstart(10^-10,pi/4,0.1,0.1,0.1,2.0,0.1))
# for i in 1.0:1.0:100.0
# println(i," ",selffsigmax(10^(-i),pi/4.,0.1,0.1,0.1,2.0,0.1,selffsigmaxstart(10^(-i),pi/4.,0.1,0.1,0.1,2.0,0.1)))
# end
#x::Float64,w::Float64,a::Float64,b::Float64,gamma::Float64,s::Float64,T::Float64,S::Float64,errors::Float64=10^-10,i::Float64=0.0)
# println(selffsigmaestart(0.1,0.1,0.0,1.0,0.1,2.0,0.1))
# println(selffsigmae(0.1,0.1,0.0,1.0,0.1,2.0,0.1,selffsigmaestart(0.1,0.1,0.0,1.0,0.1,2.0,0.1)))
# for  i in 5000.0:1.0:10000.0
# println(i,"  ",selffsigmae(0.1,0.1,0.0,i,0.1,0.1,2.0,0.1,selffsigmaestart(0.1,0.1,0.0,i,0.1,0.1,2.0,0.1)))
# end
# ss=0
# for  i in 0.0:0.0001:100
# 	global ss
# ss=ss+0.0001*fsigma(0.1,0.1,i,0.1,2.0,0.1)
# end
# println(ss)
# ss=0
# for  i in 0.0:0.000001:100
# 	global ss
# ss=ss+0.000001*fsigma(0.1,0.1,i,0.1,2.0,0.1)
# end
# println(ss)
