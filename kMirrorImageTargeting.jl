module kMirrorImageTargeting
using jBlock
using jGeom3D: Point
using LinearAlgebra: norm
using jConstants: MAX_ERR_MIRROR_IMAGE_TARGETING

import kBackwardModel, kForwardModel  #for testing test3_mit()  #2020/1/27


export mit, test1_mit, test2_mit, test_mit_all, test3_mit

#function mit(forw_model, bakw_model, targetpt, block=missing, maxiter=20, maxerr=MAX_ERR_MIRROR_IMAGE_TARGETING(); block=missing) #2020/1/25
function mit(forw_model, bakw_model, targetpt, block, maxiter=20, maxerr=MAX_ERR_MIRROR_IMAGE_TARGETING())
#generic function for mirror image targeting
  aimpt = targetpt
  xin = bakw_model(aimpt, block)   #creation of global variable
  k = 1
  while k <= maxiter
    #println("xin: ", xin)
    hitpt = forw_model(xin, block, targetpt) #targetpt added 2020/1/25
    #println("hitpt: ", hitpt)
    errvec = hitpt.vec .- targetpt.vec  #need a vector here instead of Direc  #2020/1/26
  println("errvec: ", errvec)  #temporary for testing 2020/1/27
    if norm(errvec) <= maxerr
      break
    end
    aimpt = Point(aimpt.vec .+ errvec)  #mirror image of the hit point
    xin = bakw_model(aimpt, block)
    #println("aimpt: ", aimpt)
    k += 1
  end
  if k > maxiter
    error("mirror-image-targeting failed to converge in " * string(k-1) * " iterations")
  end
  return xin
end

function mit_all(forw_model, bakw_model, targetpts, block=missing, maxiter=20, maxerr=MAX_ERR_MIRROR_IMAGE_TARGETING())

  xin = []
  for targetpt in targetpts
    push!(xin, mit(forw_model, bakw_model, targetpt, block, maxiter, maxerr))
  end

  nrows = size(targetpts,1)

  hitpts = forw_model(xin, block)
  hitpts = reshape(hitpts, nrows, :)
  println("hitpts: ", hitpts)

  xin = reshape(xin, nrows, :)
  return xin
end

function test_mit_all(targetpts, maxiter=20)
  forw_model(x) = x + 0.1 ./x
  bakw_model(x) = x
  xin = mit_all(forw_model, bakw_model, targetpts, maxiter)
end

function test3_mit(dirin, sensor)
  targetpt = Point(-99.9994,   -231.9994,  -1029.9951)  #S/N 628578 Run123, point c000 #2020/1/27
  block = Block(dirin, sensor)
  xin = mit(kForwardModel.forw_model, kBackwardModel.bakw_model, targetpt, block)
  hitpt = kForwardModel.forw_model(xin, block, targetpt)
  println("targetpt: ", targetpt)
  println("xin, hitpt: ", xin, hitpt)
  return hitpt
end

function test2_mit(maxiter=20)
  forw_model(x) = [ x[1] + .1/x[1], x[2] - .1/x[2] ]
  bakw_model(x) = [ x[1], x[2] ]
  xin = mit(forw_model, bakw_model, [2., 3.], maxiter)
  hitpt = forw_model(xin)
  println("xin, hitpt: ", xin, hitpt)
  return hitpt
end

function test1_mit()
  forw_model(x)::Float64 = x + .1/x
  bakw_model(x)::Float64 = x - .1
  xin = mit(forw_model, bakw_model, 5., 20)
  hitpt = forw_model(xin)
  println("xin: ", xin, "   hitpt: ", hitpt)
end

end  #endmodule
