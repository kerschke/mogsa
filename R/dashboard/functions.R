f1_1 = function(x) (x[1]**4 - 2*x[1]**2 + x[2]**2 + 1)
f2_1 = function(x) ((x[1] + 0.5)**2 + (x[2]-2)**2)
f3_1 = function(x) (x[1] ** 4 + 3 * (x[2] - 1) ** 2)
f_2d2d = function(x) c(f1_1(x), f2_1(x))
f_2d3d = function(x) c(f1_1(x), f2_1(x), f3_1(x))

test.2d.2d = smoof::makeMultiObjectiveFunction(name = "2D/2D Test Function", id = "", description = "", fn = f_2d2d,
                                               par.set = ParamHelpers::makeNumericParamSet(len = 2, lower = c(-2,-1), upper = c(2,3)))
test.2d.3d = smoof::makeMultiObjectiveFunction(name = "2D/3D Test Function", id = "", description = "", fn = f_2d3d,
                                               par.set = ParamHelpers::makeNumericParamSet(len = 2, lower = c(-2,-1), upper = c(2,3)))

f1_2 = function(x) (x[1]**4 - 2*x[1]**2 + x[2]**2 + 1 + x[3] ** 2)
f2_2 = function(x) ((x[1] + 0.5)**2 + (x[2]-2)**2 + (x[3] - 1) ** 4)
f3_2 = function(x) (x[1] ** 4 + 3 * (x[2] - 1) ** 2 + (x[3] - 1) ** 2)
f_3d2d = function(x) c(f1_2(x), f2_2(x))
f_3d3d = function(x) c(f1_2(x), f2_2(x), f3_2(x))
test.3d.2d = smoof::makeMultiObjectiveFunction(name = "3D/2D Test Function", id = "", description = "", fn = f_3d2d,
                                               par.set = ParamHelpers::makeNumericParamSet(len = 3, lower = c(-2,-1,-2), upper = c(2,3,2)))
test.3d.3d = smoof::makeMultiObjectiveFunction(name = "3D/3D Test Function", id = "", description = "", fn = f_3d3d,
                                               par.set = ParamHelpers::makeNumericParamSet(len = 3, lower = c(-2,-1,-2), upper = c(2,3,2)))

f1_3 = smoof::makeMPM2Function(3, 3, "random", seed=4)
f2_3 = smoof::makeMPM2Function(3, 3, "random", seed=8)
f3_3 = smoof::makeMPM2Function(1, 3, "random", seed=12)
f_3 = function(x) c(f1_3(x), f2_3(x), f3_3(x))
mpm2.fn = smoof::makeMultiObjectiveFunction(name = "MPM(3,3,random,{4,8,12})", id = "", description = "", fn = f_3,
                                       par.set = ParamHelpers::makeNumericParamSet(len = 3, lower = rep(0,3), upper = rep(1,3)))


test.functions = list(
  test.2d.2d,
  test.2d.3d,
  test.3d.2d,
  test.3d.3d,
  smoof::makeMMF10Function(),
  smoof::makeMMF6Function(),
  mpm2.fn,
  smoof::makeViennetFunction()
)

test.function.ids = as.list(1:length(test.functions))

names(test.function.ids) = lapply(test.function.ids, function(id) {
  smoof::getName(test.functions[[id]])
})
