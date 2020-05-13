

function mirrorimagetargeting(forw_model, bakw_model, targetpt, niter)

  aimpt = targetpt
  xin = bakw_model(aimpt)
  for i in 1:niter
    hitpt = forw_model(xin)
    errvec = hitpt - targetpt
    aimpt = aimpt - errvec
    xin = bakw_model(aimpt)
  end
  return xin
end
