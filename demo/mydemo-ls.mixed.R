#Play with zelig
library(Zelig)
library(ZeligMultilevel)
data(voteincome)
l5 <- lme4::lmer(income ~ education + age + female + (1 +education| state),data = voteincome, REML = T)
n5 <- nlme::nlme(income ~ B0+B1*education + B2*age + B3*female, fixed = B0+B1+B2+B3~1, random = B0+B1~1 | state, start=c(B0=8,B1=2,B2=0,B3=-.5), data = voteincome)#, method = "REML", control = list(maxIter = 50))

nlFormExp <- ~B0 + B1*(1 - exp((B2)*timeVar)) 
nlFuncExp <- deriv(nlFormExp,namevec=c("B0","B1","B2"), function.arg=c("timeVar","B0","B1","B2"))

creatExpDat <- function(nlev1=50,nlev2=10){
  rancoef <- mvrnorm(n = nlev2, c(8,1), Sigma = diag(c(1,.5),2))
  expDatList <- list()
  for(i in 1:nrow(rancoef)){
    expDatList[[i]] <- data.frame(y=nlFuncExp(1:nlev1,B0=rancoef[i,1],B1=rancoef[i,2],B2=.05)+rnorm(nlev1,sd=1),x=1:nlev1,lev2=i)
  }
  return(do.call("rbind",expDatList))
}
myExpDat <- creatExpDat()

library(ggplot2)
myCurveFunc <- function(x) nlFuncExp(x,B0=8,B1=1,B2=.05) #SSasymp(input = -x,Asym=2,R0=-1,lrc=-1)
ggplot(myExpDat, aes(x, y)) + geom_point(size = 2, colour = "darkblue", shape = 19) + theme_bw() + stat_function(fun = myCurveFunc,color="red",lwd=2) + stat_smooth(aes(group=lev2,color=lev2),size=2) 

myNlOb <- lme4::nlmer(y ~ nlFuncExp(x,B0,B1,B2) ~ B0+B1 | lev2, start=c(B0=8,B1=1,B2=.05), data = myExpDat)

nz5 <- znlsmixed$new()
nz5$zelig(formula= y ~ nlFuncExp(x,B0,B1,B2) ~ B0+B1 | lev2, start=c(B0=8,B1=1,B2=.05), data = myExpDat)
nz5$setx(x = quantile(myExpDat$x, 0.8)) 
nz5$setx1(x = quantile(myExpDat$x, 0.2))
nz5$sim()
plot(nz5)
summary(nz5)

lform <- ~B0 + B1*x
lFunc <- deriv(lform, namevec = c("B0","B1"), function.arg=c("x","B0","B1")) 
nlmer(y ~ lFunc(x,B0,B1) ~ B0 + B1 | lev2, data = myExpDat, start = c(B0=1,B1=-1))
lmer(y ~ x + (1+x | lev2), data = myExpDat, REML = FALSE)

set.seed(1234)
ls1 <- zlsmixed$new()
ls1$zelig(formula= y ~ x + (1+x | lev2), data = myExpDat, REML = FALSE)
ls1$zelig.out$z.out
ls1$setx(x = quantile(myExpDat$x, 0.8))
ls1$setx1(x = quantile(myExpDat$x, 0.2))
ls1$sim()
plot(ls1)
summary(ls1)

set.seed(1234)
ls2 <- znlsmixed$new()
ls2$zelig(formula= y ~ lFunc(x,B0,B1) ~ B0 + B1 | lev2, data = myExpDat, start = c(B0=1,B1=-1))
ls2$zelig.out$z.out
ls2$setx(x = quantile(myExpDat$x, 0.8))
ls2$setx1(x = quantile(myExpDat$x, 0.2))
ls2$sim()
plot(ls2)
summary(ls2)

# startvec <- c(Asym = 200, xmid = 725, scal = 350)
# (nm1 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
#               Orange, start = startvec))


# SSlogis(nm1@frame$age[1],fixef(nm1)[1]+ranef(nm1)$Tree[row.names(ranef(nm1)$Tree) %in% nm1@frame$Tree[1],],fixef(nm1)[2],fixef(nm1)[3])
# formula(nm1)

z5 <- zlsmixed$new()
z5$zelig(formula= income ~ education + age + female + (1 +education| state),data = voteincome)
z5$zelig.out$z.out
z5$setx(education = quantile(voteincome$education, 0.8))
z5$setx1(education = quantile(voteincome$education, 0.2))
z5$sim()
plot(z5)
summary(z5)

length(z5$getcoef())

sample(1:10,size = 10, replace = T)

z5$zelig.out$z.out


myzeligList <- list()
#myzeligList$zelig.call <- as.call(myzeligList$zelig(formula = income ~ education + age + female + (1 + education | state), data = voteincome))
myzeligList$model.call <- l5@call
myzeligList$formula <- as.formula(income ~ education + age + female + (1 +education| state))
myzeligList$originaldata <- voteincome
myzeligList$bootstrap <- FALSE
myzeligList$model.call$bootstrap <- NULL
myzeligList$matched  <- FALSE
myzeligList$mi <- FALSE
myzeligList$data <- voteincome
myzeligList$weights <- NULL # No valid weights
myzeligList$model.call$weights <- NULL

myzeligList$model.call[[1]] <- quote(lme4::lmer) #.self$fn
myzeligList$model.call$by <- NULL

myzeligList$data <- cbind(1, myzeligList$data)
names(myzeligList$data)[1] <- "by"
myzeligList$by <- "by"
myzeligList$data <- tbl_df(myzeligList$data)
fn2 <- function(fc, data) {
  fc$data <- data
  return(fc)
}
myzeligList$zelig.out <- myzeligList$data %>% group_by_(myzeligList$by) %>% do(z.out = eval(fn2(myzeligList$model.call, quote(as.data.frame(.)))))
myzeligList$zelig.out$z.out

myzeligList$formula.full <- myzeligList$formula # fixed and random effects
myzeligList$formula <- formula(myzeligList$zelig.out$z.out[[1]], fixed.only = TRUE)

myzeligList$bsetx <- TRUE
#myzeligList$setx.out$x  <- myzeligList$set(..., fn = fn)

f2 <- myzeligList$formula
f <- update(myzeligList$formula, 1 ~ .) #this won't work for nlme
formula(n5)

myzeligList$setforeveryby <- TRUE
update <- myzeligList$data %>%
  group_by_(myzeligList$by) %>%
  do(mm = model.matrix(f, ldata))
update$mm

myzeligList$simparam <- myzeligList$zelig.out %>%
  do(simparam = myzeligList$param(.$z.out))
myzeligList$zelig.out$z.out

mvrnorm(5, coef(myzeligList$zelig.out$z.out[[1]]), vcov(myzeligList$zelig.out$z.out[[1]]))

vcov(n5)

myzeligList$zelig.out %>% do(simparam = mvrnorm(5,coef(.$z.out),vcov(.$z.out))) 
chkMe <- myzeligList$zelig.out %>% do(simparam = arm::sim(.$z.out,5)) 
chkMe$simparam
arm::sim(myzeligList$zelig.out$z.out[[1]], 5)
Lz <- getME(myzeligList$zelig.out$z.out[[1]], "L")
Rzx <- getME(regression, "RZX")
Rx <- getME(regression, "RX")


list(simparam = arm::sim(z.out, .self$num), simalpha = z.out)
arm::sim(n5, 5)
myArmSim(myNlOb, 5)

RLRsim::extract.lmeDesign()
class (n5)
showMethods("sim")
getMethod("sim","merMod")
#reduce(dataset = "MEANINGLESS ARGUMENT", s, formula = f2, data = ., avg = .self$avg))

reduce(dataset = "MEANINGLESS ARGUMENT", s=list(education = quantile(voteincome$education, 0.8)), formula = f2, data = ., avg = .self$avg)

z5$setx.out$x$mm



pred <- try(terms(fit <- lm(f2, myzeligList$data), "predvars"), silent = TRUE) #probably won't work for nlme
dataset <- model.frame(fit) #won't work for nlme
ldata <- lapply(dataset, myavg)
n <- union(as.character(attr(pred, "predvars"))[-1], names(dataset))
if (is.list(s[[1]])) s <- s[[1]]
m <- match(names(s), n)
ma <- m[!is.na(m)]
for (i in seq(n[ma])) ldata[n[ma]][i][[1]] <- setval(dataset[n[ma]][i][[1]], s[n[ma]][i][[1]])


myavg <- function(val, fn = list(numeric = mean, ordered = Median, other = Mode)) {
  if (is.numeric(val))
    ifelse(is.null(fn$numeric), mean(val), fn$numeric(val))
  else if (is.ordered(val))
    ifelse(is.null(fn$ordered), Median(val), fn$ordered(val))
  else
    ifelse(is.null(fn$other), Mode(val), fn$other(val))
}

z5$setx(education = quantile(voteincome$education, 0.8))

z5$model.call
class(voteincome)
z5$bootstrap
is.logical(z5$bootstrap)
as.character(z5$fn)[2]

is.character(z5$mcformula)
nlme::getGroups()
chkMe <- getME(z5$zelig.out$z.out[[1]], "devcomp")
chkMe$dims
nlme::gsummary(n5)
n5$dims$Q
dims <- list(n=n5$dims$N,q=n5b$dims$ngrps[1]*n5b$dims$qvec[1],reTrms=n5b$dims)
?nlmeObject
n5b$dims$ngrps[1]*n5b$dims$qvec[1]

l5b <- lme4::lmer(income ~ education + age + female + (1+education | mystate) + (1 | state),data = myvote)
getME(l5b, "devcomp")$dims
head(getME(l5b, "X"))

incidence <-
  matrix(c(1,1,2,2,2,1,3,2,0), ncol = 3, byrow = TRUE,
         dimnames = list(NULL, c("Par", "From", "To")))
incidence
nlme::LDEsysMat(c(1.2, 0.3, 0.4), incidence)



z5$setx.out$x1$mm
head(z5$setx.out$x$by)
X <- getME(l5, "X");
head(X)
X <- matrix(rep(mm, length(X)), nrow(X), ncol(X), byrow = TRUE)

head(getME(l5, "Zt"))

model.matrix(n5$modelStruct$reStruct)

n5$groups$state
model.matrix(formula(n5$modelStruct$reStr)[[1]],data=nlme::getData(n5))

head(nlme::getData(n5))

## Linear predictor
x.beta <- as.matrix(tcrossprod(as.matrix(X), sims@fixef))



myvote <- voteincome
1500/4
myvote$mystate <- rep(c("AR","SC","SC2","SC3"),each=375)
n5b <- nlme::nlme(income ~ B0+B1*education + B2*age + B3*female, fixed = B0+B1+B2+B3~1, random = B0~1 | mystate/vote, start=c(B0=0,B1=0,B2=0,B3=0), data = myvote)

str(l5)
s <- list()
intersect(names(s), group)
# data(coalition)
# z.out <- zelig(duration ~ fract + numst2, model = "gamma", data = coalition)
# summary(z.out)

names(ranef(z5$zelig.out$z.out[[1]]))
lme4::ranef(l5)
group <- names(n5$coefficients$random)
head(n5$coefficients$random)
names(n5$groups)
n5$call 
table(voteincome$state)
head(nlme::random.effects(n5))


formula(z5$zelig.out$z.out[[1]], fixed.only = TRUE)
formula(l5, fixed.only = TRUE)
formula(n5, fixed.only = TRUE)
nlme::getGroups(n5)
nlme::getCovariateFormula(n5)
nlme::getGroupsFormula(n5)
nlme::getResponseFormula(n5)
nlme::getResponse(n5)
nlme::getCovariate(n5)
z5$formula
n5$formula

head(model.matrix(z5$formula, data=voteincome))
head(model.matrix(formula(l5, fixed.only = TRUE), data=voteincome))
head(model.matrix(n5, data=voteincome))

n5$modelStruct$reStruct$state
formula(n5$modelStruct$reStruct$state)
n5$plist
model.matrix(formula(n5$modelStruct$reStruct),data=voteincome)
formula(n5)

f <- update(z5$formula, 1 ~ .)

pd1 <- nlme::pdSymm(~Sex*age)
formula(pd1)

pd1 <- nlme::pdMat(diag(1:4), pdClass = "pdDiag")
pd1

formula(n5$call$fixed)
update(n5$call$fixed, 1 ~ .)
n5$dims
head(n5$fitted)

nlme:::nlme
getMethods("nlme","nlme")
methods("lme")
nlme:::nlme.formula

asOneSidedFormula(n5$call$model)
nlme::asOneFormula(y ~ x + z | g, list(~ w, ~ t * sin(2 * pi)))
all.vars(expression(sin(x+y)))
all.vars(n5$call$model)
asOneSidedFormula(~ age)

#TODO:: not sure how or if it's necessary to solve the formula issue

