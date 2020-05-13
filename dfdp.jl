function dfdp(block::Block,x,f,p,dp,func)
# numerical partial derivatives [Jacobian] df/dp for use with leasqr
# --------INPUT VARIABLES---------
# x=vec or matrix of indep var(used as arg to func) x=[x0 x1 ...]
# f=func[x,p] vector initialsed by user before each call to dfdp
# p= vec of current parameter values
# dp= fractional increment of p for numerical derivatives
#      dp[j]>0 central differences calculated
#      dp[j]<0 one sided differences calculated
#      dp[j]=0 sets corresponding partials to zero; i.e. holds p[j] fixed
# func=string naming the function (.m) file
#      e.g. to calc Jacobian for function expsum prt=dfdp[x,f,p,dp,"expsum"]
#----------OUTPUT VARIABLES-------
# prt= Jacobian Matrix prt[i,j]=df[i]/dp[j]
#---------------------------------

m=size(x,1)
if m==1
    m=size(x,2)
end  ## PAK: in case() #cols > #rows
n=length(p)      #dimensions
ps=p; prt=zeros(m,n); del=zeros(n,1);       # initialise Jacobian to Zero
for j=1:n
      del[j]=dp[j] .*p[j];    #cal delx=fract[dp]*param value[p]
      if p[j]==0
           del[j]=dp[j];     #if param=0 delx=fraction
      end
      p[j]=ps[j] + del[j]
      if del[j]!=0
           #f1=feval(func,x,p)
           f1=func(block,x,p)  #for Julia
           if dp[j] < 0 prt[:,j]=(f1-f)./del[j]
           else
                p[j]=ps[j]- del[j]
                #prt[:,j]=(f1-feval(func,x,p))./(2 .*del[j])
                prt[:,j]=(f1-func(block,x,p))./(2 .*del[j])  #for Julia
           end
      end
      p[j]=ps[j];     #restore p[j]
end
    #println("prt", size(prt))  #for debugging
return prt  #for Julia
    #return prt'  #for test 9/9/2018
end


#******************************************************************
