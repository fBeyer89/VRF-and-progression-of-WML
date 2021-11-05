library("DHARMa")


##Three different ways to fit lesion load for baseline
usable=subset(bl,(bl$Age_all>50&(!is.na(bl$ADULT_BP_SBP.1))&(!is.na(bl$waist2hip))&(!is.na(bl$icv))))

fit1 <- lm(log10(usable$lesionload/1000) ~ Age_all, data = usable)
fit3 <- glm(lesionload/1000 ~  Age_all, data = usable, family = gaussian(link = "log"))
fit5 <- glm(lesionload/1000 ~ Age_all, data = usable, family = Gamma(link = "log"))

#residuals are not normally distributed for fit1
diagnostics.plot(fit1)
plot(fitted(fit1),residuals(fit1))
hist(residuals(fit1))
#residuals are not without pattern for fit3
scatter.smooth(predict(fit3, type='response'), rstandard(fit3, type='deviance'), col='gray')
plot(fitted(fit3),residuals(fit3))
#residuals are largely without pattern for fit5
scatter.smooth(predict(fit5, type='response'), rstandard(fit5, type='deviance'), col='gray')
plot(fitted(fit5),residuals(fit5))
check_gamma_model <- simulateResiduals(fittedModel = fit5, n = 500)
plot(check_gamma_model)

plot(lesionload/1000 ~ Age_all, data = usable)
curve(10**(predict(fit1,newdata = data.frame(Age_all = x))), col = "green", add = TRUE)
curve(predict(fit3, newdata = data.frame(Age_all = x), type = "response"), col = "red", add = TRUE, lty = 2)
curve(predict(fit5, newdata = data.frame(Age_all = x), type = "response"), col = "cyan", add = TRUE)
      
coefficients(fit1)
coefficients(fit3)
coefficients(fit5)


qqplot(residuals(fit5))

plot(log(bl$lesionload/1000) ~ Age_all, data = bl)
curve(predict(fit1,newdata = data.frame(Age_all = x)), col = "green", add = TRUE)




summary(fit3)

######
usable=subset(bl,(bl$Age_all>50&(!is.na(bl$ADULT_BP_SBP.1))&(!is.na(bl$waist2hip))&(!is.na(bl$icv))))
usable$age_sc=usable$Age_all-mean(usable$Age_all)
#+ scale(ADULT_BP_SBP.1, scale = F) +
#  scale(waist2hip, scale = F)
test=glm(formula = lesionload/1000 ~ age_sc , 
    family = Gamma(link = "log"), data = usable) # + sex + scale(icv,scale = F)

test2=lm(log(lesionload/1000) ~ age_sc , data = usable) # + sex + scale(icv,scale = F)

plot(usable$age_sc, usable$lesionload/1000)

curve(predict(test, newdata=data.frame(age_sc=x)),
                                       col = "green", add = TRUE)
curve(exp(predict(test2, newdata=data.frame(age_sc=x))),
      col = "red", add = TRUE)



