"""
function solve_quadratic(a,b,c)
#Solve a quadratic equation given the coefficients
#(allows complex solutions).
"""
function solve_quadratic(a,b,c)
    if a==0 #not quadratic--only one solution
        (x1, x2) = (-c/b, -c/b)
    #elseif c==0 #solve by factoring
    #    (x1, x2) = (-b/a, 0)
    else
        discr = b^2 - 4a*c
        if discr>=0
            s = sqrt(discr)
        else
            s = sqrt(complex(discr))
        end
        (x1, x2) = ((-b + s)/2a, (-b - s)/2a)
    end
    #check
    #println(s)
    #println(x1, x2)
    #println(abs(a*x1^2 + b*x1 +c))
    maxerr = 1.e-8
    if abs(a*x1^2 + b*x1 +c) > maxerr ||
      abs(a*x2^2 + b*x2 +c) > maxerr
        #println("x1 error: ", abs(a*x1^2 + b*x1 +c))
        #println("x2 error: ", abs(a*x2^2 + b*x2 +c))
        error("solve_quadratic failed")
    end
    return (x1, x2)
end
